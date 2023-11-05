#%%
import shutil
import joblib
import numpy as np
import yaml
from ship_plasmid.utils.motifs import BulkMotifFinder
from ship_plasmid.utils.files import get_amr
from ship_plasmid.utils.genome import FastBreakDistance
from ship_plasmid.utils.phylogeny import PlasmidDistance
from ship_plasmid.utils.amrfinder import make_amrfinder_df
from ship_plasmid.utils.cdhit import process_cdhit
import matplotlib.pyplot as plt
import os
from ship_plasmid.utils.clustering import AdaptativeAgglomerativeClustering
import datetime
import argparse
import logging

def get_gff_filenames(annotations_dir: str):
    '''
    Retrieves the plasmid sequence names through the annotations directory.
    Should be the name of the folder containing each individual Prokka output. 
    '''
    return list(os.walk(annotations_dir))[0][1]

def main():

    parser = argparse.ArgumentParser(description = 'SHIP')
    parser.add_argument(
        '--annotations', '-a', nargs = 1, required = True,
        help = 'Directory containing the output from Prokka.'
    )
    parser.add_argument(
        '--out', '-o', nargs = 1, required = True,
        help = 'Output directory.'
    )
    parser.add_argument(
        '--cdhit', '-c', nargs = 1, required = True,
        help = 'Path to the .clstr file resulting from CD-HIT. Must contain information about all annotated plasmids.'
    )
    parser.add_argument(
        '--amr', '-r', nargs = 1, required = True,
        help = 'Path to the directory containing the AMRFinder+ output files. All .txt files in this directory will be considered as the main AMRFinder+ output files.'
    )
    parser.add_argument(
        '--plot-dendrogram', '-d', action='store_true',
        help = 'If set, plots a dendrogram of the plasmid phylogeny as estimated by SHIP.'
    )
    parser.add_argument('--min_dist', '-m', default = 0.1, nargs = 1,
        help = 'Minimum average plasmid distance between plasmids containing a region for it to be considered as having evidence for HGT. Default is 0.1.'
    )
    parser.add_argument('--min_len', '-l', default = 5, nargs = 1,
        help = 'Minimum number of CDS in a fragment with evidence for HGT for it to be included. Default is 5.'
    )
    parser.add_argument('--max_len', '-L', default = 9, nargs = 1,
        help = 'Maximum searched number of CDS in fragments. Longer fragments with evidence for HGT will be split in fragments of --max_len. Default is 9.'
    )
    parser.add_argument('--min_n', '-n', default = 3, nargs = 1,
        help = 'Minimum number of plasmids containing a region with evidence for HGT for it to be included. Default is 3'
    )
    parser.add_argument('--keep-intermediate', '-i', action='store_true',
        help = 'If set, keeps all intermediate files in path_to_SHIP_output/tmp.')

    args = parser.parse_args()
    args.out = args.out[0]
    args.cdhit = args.cdhit[0]
    args.amr = args.amr[0]
    args.annotations = args.annotations[0]

    # Prepare logger
    logging.basicConfig(filename=os.path.join(args.out, 'ship.log'), filemode='w', format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=50)
    logging.info(f'Starting SHIP.')
    logging.info(f'Output directory: {args.out}\nReading annotations from {args.annotations}\nUsing CH-HIT clusters at {args.cdhit}')
    with open('configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)
    logging.info(f"Phylo config file (at configs/phylo_config):") 
    logging.info('\n'.join([k+': '+str(v) for k, v in phylo_config.items()]))


    timestamp = datetime.datetime.now().strftime('%d-%b-%Y__%H-%M-%S')

    #########################################
    #            Prepare inputs             #
    #########################################

    tmp_path = os.path.join(args.out, 'tmp')
    logging.debug(f'Creating tmp directory at {tmp_path}.')
    if not os.path.exists(tmp_path): os.mkdir(tmp_path)

    # Parse AMRFinder+ output
    logging.info('Processing the AMRFinder+ output files.')
    path_to_amr_tsv = os.path.join(tmp_path, 'amrfinder_output.tsv')
    make_amrfinder_df(args.amr, path_to_amr_tsv)

    # Parse CD-HIT output
    logging.info('Processing the CD-HIT output file.')
    process_cdhit(args.cdhit, args.annotations, tmp_path)

    #########################################
    #  Build a plasmid dissimilarity matrix #
    #########################################

    if phylo_config['agglomerative-clustering']['adaptative']: ratio_ = phylo_config['agglomerative-clustering']['distance-threshold']
    else: ratio_ = None

    clustering_method = AdaptativeAgglomerativeClustering(
        n_clusters = None,
        affinity = 'precomputed',
        linkage = phylo_config['agglomerative-clustering']['linkage'],
        distance_threshold = phylo_config['agglomerative-clustering']['distance-threshold'],
        distance_ratio = ratio_,
        compute_distances = True
    )

    distance_function = FastBreakDistance(
        phylo_config['break-distance']['k'],
        phylo_config['break-distance']['indel-size-penalty']
    )

    # Get the name of the plasmid sequences by checking the subdirectories in the annotations directory
    logging.debug('Retrieving plasmid sequence names through the annotations directory.')
    plasmid_names = get_gff_filenames(args.annotations)

    # Get plasmid phylogeny
    logging.debug('Computing plasmid distances.')
    phylo = PlasmidDistance(plasmid_names, args.annotations, os.path.join(tmp_path, 'protein_cluster_membership.tsv'),
        path_to_amr_tsv, clustering_method = clustering_method, learn_weights = phylo_config['learn-weights'],
        learning_method = phylo_config['learning-method'], weights = phylo_config['weights'])

    clusters = phylo.fit_predict(distance_function)
    clusters = phylo.get_clusters_as_series()
    affinity = phylo.get_affinity_as_frame()

    if args.plot_dendrogram:
        logging.info('Plotting the plasmid phylogeny dendrogram.')
        phylo.plot_dendogram()
        plt.show()

    # Store plasmid distances
    plasmid_distance_path = os.path.join(args.out, 'plasmid_dissimilarity.pkl')
    logging.info(f'Storing plasmid distance/dissimilarity (PlasmidDistance object) in {plasmid_distance_path}')
    joblib.dump(phylo, plasmid_distance_path, compress=3)

    #########################################
    #  Find conserved regions in dissimilar #
    #              plasmids                 #
    #########################################

    # Convert args to kwargs-like dict
    regions_config = {'min-distance': float(args.min_dist[0]), 'min-length': int(args.min_len[0]), 
        'max-length': int(args.max_len[0]), 'min-n-plasmids': int(args.min_n[0])}

    # Load information about AMR genes
    logging.debug('Building AMR DataFrame with AMRFinder+ output.')
    amr_df = get_amr(plasmid_names, path_to_amr_tsv, args.annotations)
    logging.debug('Getting the IDs of resistant plasmids.')
    resistant_ids = np.unique(amr_df[amr_df['Gene symbol']!= 'Susceptible'].index)

    logging.debug('Instanciating a BulkMotifFinder object.')
    bmf = BulkMotifFinder(resistant_ids, phylo, tmp_path, args.annotations)
    logging.info('Searching for AMR regions with evidence for HGT.')
    result = bmf.search(regions_config['min-distance'], regions_config['min-length'],
        regions_config['max-length'], regions_config['min-n-plasmids'])
        
    if result is not None:
        logging.info('Saving the search results.')
        bmf.save_report(os.path.join(args.out, 'HGT_AMR_regions.tsv'))
        joblib.dump(bmf, os.path.join(args.out, 'HGT_AMR_finder.pkl'), compress = 3)

    if not args.keep_intermediate: 
        logging.info('Removing temporary files.')
        shutil.rmtree(tmp_path)

    logging.info('Done!')

if __name__ == '__main__':
    main()