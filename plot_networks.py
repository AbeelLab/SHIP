'''
Plots plasmid networks obtained with the proposed method. Shows the estimated frequencies
for evolutionary events. Optionally, shows networks colored on plasmid types.
'''

from sklearn.cluster import AgglomerativeClustering
from utils.phylogeny import Pangenome
from utils.genome import GraphGenome
from utils.plasmid_typing import LabeledNetwork
import yaml
import pandas as pd
import numpy as np
from utils.cluster_images import PlasmidNetwork, SubClusterImages
import joblib
import os
import warnings
import matplotlib.cm as cm
import argparse

parser = argparse.ArgumentParser(
    description = '''
Plots plasmid similarity networks obtained with the proposed method. Shows the estimated frequencies
for evolutionary events. Optionally, shows networks colored on plasmid types.
'''
)
parser.add_argument(
    '--types',
    action = 'store_true',
    help = 'If set, plots plasmid networks colored on plasmid types from MOBSuite.'
)
parser.add_argument(
    '--paper',
    action = 'store_const',
    const = True, 
    default = False, 
    help = 'Loads the Jaccard-based clusters used in the original paper.'
)
parser.add_argument(
    '--panplasmidomes',
    action = 'store_const',
    const = True, 
    default = False, 
    help = 'Shows the panplasmidome graphs for each subcluster in the networks.'
)
args = parser.parse_args() 

if __name__ == '__main__':
    with open('configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)

    if args.paper:
        phylo_config['results-dir'] = os.path.join(
            phylo_config['results-dir'],
            'Paper'
        )
        data_config['paths']['mob-types'] = os.path.join(
            data_config['paths']['mob-types'],
            'Paper'
        )

    for k in phylo_config['output-paths']:
        phylo_config['output-paths'][k] = os.path.join(
            phylo_config['results-dir'],
            phylo_config['output-paths'][k]
        )

    clusters = pd.read_csv(
        phylo_config['output-paths']['sub-clusters'],
        index_col = 0,
        sep = '\t'
    )
    SHOW_PANGENOME = args.panplasmidomes
    #%%
    for cluster in np.unique(clusters['Megacluster']):

        print(
    f'''
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                            CLUSTER {cluster}
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    '''
        )

        subclusters = clusters.loc[clusters['Megacluster']==cluster]
        phylo = joblib.load(phylo_config['output-paths']['phylogenies'])[cluster]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            genomes = {
                x: GraphGenome(
                    x,
                    data_config['paths']['annotations'],
                    phylo_config['paths']['representative-proteins'],
                    data_config['paths']['amr_hits']
                )
                for x in phylo.accessions
            }

        plots = SubClusterImages(
            phylo,
            genomes,
            data_config['paths']['info_df'],
            data_config['paths']['amr_hits'],
            phylo_config['output-paths']['linkage']+str(cluster)+'.pkl',
            figsize=(15,15),
            distance_threshold = phylo.clustering_method.distance_threshold
        )
        plots.show(subclusters, ratio=True, jaccard_distance_correlation=True,size_correlation=False, sizes=False)

        if SHOW_PANGENOME:
            for subcluster in np.unique(subclusters['Cluster']):

                ids = subclusters[subclusters['Cluster']==subcluster].index.to_list()
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    genomes = {
                        x: GraphGenome(
                            x,
                            data_config['paths']['annotations'],
                            phylo_config['paths']['representative-proteins'],
                            data_config['paths']['amr_hits']
                        )
                        for x in ids
                    }
                
                pangenome = Pangenome(
                    list(genomes.values()),
                    path_to_protein_names = phylo_config['paths']['protein-names'],
                    path = None
                )
                pangenome.show(with_duplicates = True)
                print(f'\nShowing cluster {cluster}, subcluster {subcluster}.\nAccessions: {ids}\n')
                input('Press any key to continue.')

        graph = PlasmidNetwork(
            data_config['paths']['info_df'],
            edge_threshold = .5,
            hide_threshold = .5
        )

        graph.show(
            phylo.get_clusters_as_series(),
            phylo.get_affinity_as_frame(),
            path = None
        )
        AgglomerativeClustering()

    # %%
    if args.types:
        for m, (method, name_) in enumerate(zip(
            [
                'Mate-Pair Formation',
                'Origin of Transfer',
                'Replicon',
                'Relaxase',
                'Mobility'
            ],
            [
                'mpfs',
                'orit',
                'rep',
                'relaxases',
                'mobility'
            ]
        )):
            print(f'''
        ::::::::::: {method} Typing :::::::::::
            ''')

            for n, cluster in enumerate(np.unique(clusters['Megacluster'])):

                phylo = joblib.load(phylo_config['output-paths']['phylogenies'])[cluster]

                print(
        f'''
        ::::::::::: CLUSTER {cluster} :::::::::::                     
        '''
                )
                typings = pd.read_csv(
                    os.path.join(
                        data_config['paths']['mob-types'],
                        f'MOBSuite_{name_}.tsv'
                    ),
                    sep = '\t',
                    index_col = 0
                )
                network = LabeledNetwork(
                    distance_threshold=.5,
                    colors = np.vstack([[
                        np.array([33/255, 111/255, 221/255, 1]),
                        np.array([219/255, 132/255, 11/255, 1]),
                        np.array([73/255, 190/255, 227/255, 1]),
                        np.array([100/255, 219/255, 22/255, 1]),
                        np.array([192/255, 0/255, 0/255, 1]),
                        np.array([127/255, 127/255, 127/255, 1]),
                        np.array([250/255, 250/255, 0/255, 1])
                    ], cm.rainbow(np.linspace(0,1,100))])
                )
                network.from_distance(
                    phylo.get_affinity_as_frame(),
                    typings
                )
                network.show()
                input()
                metrics = network.get_metrics(phylo.get_clusters_as_series())

                if n == 0:
                    silhouette = metrics['Silhouette Score']
                    ami = metrics['Adjusted Mutual Information']
                else:
                    silhouette = pd.concat(
                        [silhouette, metrics['Silhouette Score']],
                        axis = 'columns'
                    )
                    ami = pd.concat(
                        [ami, metrics['Adjusted Mutual Information']],
                        axis='columns'
                    )
            silhouette.columns = np.unique(clusters['Megacluster'])
            ami.columns = np.unique(clusters['Megacluster'])

            if m == 0:
                mean_silhouette = silhouette.mean(skipna = True)
                mean_ami = ami.mean(skipna = True)
            else:
                mean_silhouette = pd.concat(
                    [
                        mean_silhouette,
                        silhouette.mean(skipna = True)
                    ],
                    axis = 'columns'
                )
                mean_ami = pd.concat(
                    [
                        mean_ami,
                        ami.mean(skipna = True)
                    ],
                    axis = 'columns'
                )
                input()

        for x in [mean_silhouette, mean_ami]:
            x.columns = [
                    'Mate-Pair Formation',
                    'Origin of Transfer',
                    'Replicon',
                    'Relaxase',
                    'Mobility'
                ]
    # %% Frequency estimates for each cluster
    frequencies = {}
    for n, cluster in enumerate(np.unique(clusters['Megacluster'])):

            phylo = joblib.load(phylo_config['output-paths']['phylogenies'])[cluster]
            frequencies[cluster+1] = {
                k.capitalize(): (1/v)/(sum(
                    [
                        1/phylo.weights[kk] for kk in [
                            'reversal', 
                            'duplicate', 
                            'mobile element', 
                            'transposition', 
                            'mutation'
                        ]
                    ]
                ))*100 
                for k,v in phylo.weights.items()
                if k != 'recombination'
            }

    frequencies = pd.DataFrame.from_dict(
        frequencies
    ).transpose().round(1)
    print('Estimated frequencies of evolutionary events:')
    print(frequencies)

# %%
