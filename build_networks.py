#%%
from copy import deepcopy
import warnings
import joblib
import yaml
from utils.genome import GraphGenome, FastBreakDistance
from utils.phylogeny import Pangenome, SimpleDendogram
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from utils.clustering import AdaptativeAgglomerativeClustering
import datetime
from utils.jaccard import jaccard_distance
import argparse

parser = argparse.ArgumentParser(
    description = '''
Creates plasmid similarity networks based on the proposed method. Stores a SimpleDendrogram object
that can be used for further analysis.
'''
)
parser.add_argument(
    '--clusters',
    nargs = '+',
    required = False,
    help = 'Jaccard cluster numbers for which to compute plasmid networks'
)
parser.add_argument(
    '--paper',
    action = 'store_const',
    const = True, 
    default = False, 
    help = 'Loads the Jaccard-based clusters used in the original paper.'
)

if __name__ == '__main__':
    args = parser.parse_args()
    with open('configs/clustering_config.yaml', 'r') as config_file:
        clustering_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)

    timestamp = datetime.datetime.now().strftime('%d-%b-%Y__%H-%M-%S')
    home_dir = os.path.join(
        phylo_config['output']
    )
    for k in phylo_config['output-paths']:
        phylo_config['output-paths'][k] = os.path.join(
            home_dir,
            phylo_config['output-paths'][k]
        )

    if phylo_config['agglomerative-clustering']['adaptative']:
        ratio_ = phylo_config['agglomerative-clustering']['distance-threshold']
    else:
        ratio_ = None

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

    all_clusters = pd.DataFrame(
        [],
        columns = ['Cluster', 'Megacluster']
    )
    all_phylos = {}

    if args.paper:
        clusters_input = phylo_config['clusters']
        phylo_config['paths']['plasmid-clusters'] = os.path.join(
            *(phylo_config['paths']['plasmid-clusters'].split('/')[:-1]+['Paper', 'Clusters.pkl'])
        )
    else:
        clusters_input = [int(x) for x in args.clusters]

    for cluster in clusters_input:

        plasmid_accessions = joblib.load(phylo_config['paths']['plasmid-clusters']).loc[
            cluster
        ]['Plasmids']

        phylo = SimpleDendogram(
            plasmid_accessions,
            data_config['paths']['annotations'],
            phylo_config['paths']['representative-proteins'],
            data_config['paths']['amr_hits'],
            clustering_method = clustering_method,
            learn_weights = phylo_config['learn-weights'],
            learning_method = phylo_config['learning-method'],
            weights = phylo_config['weights']
        )

        clusters = phylo.fit_predict(distance_function)
        clusters = phylo.get_clusters_as_series()
        affinity = phylo.get_affinity_as_frame()
        phylo.plot_dendogram()
        plt.show()

        all_clusters = pd.concat(
            [
                all_clusters,
                pd.concat(
                    [
                        clusters,
                        pd.Series(
                            [cluster]*len(clusters),
                            index = clusters.index,
                            name = 'Megacluster'
                        )
                    ],
                    axis = 1
                )
            ]
        )
        all_phylos[cluster] = deepcopy(phylo)

        pangenome = Pangenome(
            [
                GraphGenome(
                    x,
                    data_config['paths']['annotations'],
                    phylo_config['paths']['representative-proteins'],
                    data_config['paths']['amr_hits']
                )
                for x in phylo.accessions
            ],
            path_to_protein_names = phylo_config['paths']['protein-names'],
            path = phylo_config['output-paths']['pangenomes']+str(cluster)+'.html'
        )
        pangenome.show()

        joblib.dump(
            [deepcopy(phylo.linkage_matrix), deepcopy(phylo.accessions)],
            phylo_config['output-paths']['linkage']+str(cluster)+'.pkl'
        )

    all_clusters.to_csv(
        phylo_config['output-paths']['sub-clusters'],
        sep = '\t'
    )
    joblib.dump(
        all_phylos,
        phylo_config['output-paths']['phylogenies']
    )

    with open(phylo_config['output-paths']['logfile'], 'x') as outstream:
        outstream.write(
            '\n'.join(
                [k+': '+str(v) for k, v in phylo_config.items()]
            )
        )
