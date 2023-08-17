#%%
import joblib
import numpy as np
import yaml
from utils.motifs import BulkMotifFinder
from utils.files import get_amr
from utils.genome import FastBreakDistance
from utils.phylogeny import PlasmidDistance
import matplotlib.pyplot as plt
import os
from utils.clustering import AdaptativeAgglomerativeClustering
import datetime
import argparse
import logging

parser = argparse.ArgumentParser(
    description = '''
Creates plasmid similarity networks based on the proposed method. Stores a SimpleDendrogram object
that can be used for further analysis.
'''
)
parser.add_argument(
    '--annotations', '-a', nargs = 1, required = True, const = True,
    help = 'Directory containing the output from Prokka.'
)
parser.add_argument(
    '--out', '-o', nargs = 1, required = True, const = True,
    help = 'Output directory.'
)
parser.add_argument(
    '--cdhit-clusters', '-c', nargs = 1, required = True, const = True,
    help = 'Path to the TSV file containing the homolog cluster assigned to each protein by CD-HIT. Can be obtained with replace_with_homologs.py.'
)
parser.add_argument(
    '--amr', '-r', nargs = 1, required = True, const = True,
    help = 'Path to the AMRFinder+ output as a TSV file. This file can be obtained from the original AMRFinder+ output with make_amrfinder_df.py.'
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

def get_gff_filenames(annotations_dir: str):
    '''
    Retrieves the plasmid sequence names through the annotations directory.
    Should be the name of the folder containing each individual Prokka output. 
    '''
    return list(os.walk(annotations_dir))[0][1]

if __name__ == '__main__':

    args = parser.parse_args()

    # Prepare logger
    logging.basicConfig(filename=os.path.join(args.out, 'ship.log'), filemode='w', format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logging.info(f'Starting SHIP.')
    logging.info(f'Output directory: {args.out}\nReading annotations from {args.annotations}\nUsing CH-HIT clusters at {args.cdhit_clusters}')
    with open('configs/clustering_config.yaml', 'r') as config_file:
        clustering_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)
    logging.info(f"Phylo config file (at configs/phylo_config):") 
    logging.info('\n'.join([k+': '+str(v) for k, v in phylo_config.items()]))


    timestamp = datetime.datetime.now().strftime('%d-%b-%Y__%H-%M-%S')

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
    phylo = PlasmidDistance(plasmid_names, args.annotations, args.cdhit_clusters,
        args.amr, clustering_method = clustering_method, learn_weights = phylo_config['learn-weights'],
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
        'max-length': int(args.max_len[0]), 'min-n-plasmids': int(args.min_n[0]), 
        'motif-finder-output-dir': args.out}

    # Load information about AMR genes
    logging.debug('Building AMR DataFrame with AMRFinder+ output.')
    amr_df = get_amr(plasmid_names, args.amr)
    logging.debug('Getting the IDs of resistant plasmids.')
    resistant_ids = np.unique(amr_df[amr_df['Gene symbol']!= 'Susceptible'].index)

    #TODO: Check the files BMF is using through the paths in the configs
    logging.debug('Instanciating a BulkMotifFinder object.')
    bmf = BulkMotifFinder(resistant_ids, phylo, data_config, phylo_config)
    logging.info('Searching for AMR regions with evidence for HGT.')
    result = bmf.search(regions_config['min-distance'], regions_config['min-length'],
        regions_config['max-length'], regions_config['min-n-plasmids'])
        
    if result is not None:
        logging.info('Saving the search results.')
        bmf.save_report(os.path.join(args.out, 'HGT_AMR_regions.tsv'))
        joblib.dump(bmf, os.path.join(args.out, 'HGT_AMR_finder.pkl'), compress = 3)

    logging.info('Done!')
