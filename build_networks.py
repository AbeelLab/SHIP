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
    '--description',
    default = 'Plasmid_Similarity',
    nargs = 1,
    help = 'Name of the output files. The current date and time will be added after the description. Default is "Plasmid_Similarity".'
)
args = parser.parse_args() 

if __name__ == '__main__':
    with open('./configs/clustering_config.yaml', 'r') as config_file:
        clustering_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('./configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('./configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)

    timestamp = datetime.datetime.now().strftime('%d-%b-%Y__%H-%M-%S')
    home_dir = os.path.join(
        phylo_config['output'],
        args.description[0]+'_'+timestamp
    )
    for k in phylo_config['output-paths']:
        phylo_config['output-paths'][k] = os.path.join(
            home_dir,
            phylo_config['output-paths'][k]
        )
    os.mkdir(home_dir)

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

    for cluster in phylo_config['clusters']:

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
