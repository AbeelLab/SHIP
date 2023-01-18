'''
Clusters plasmids according to their Jaccard distance on gene content.
'''
import os
import joblib
from matplotlib import pyplot as plt
import pandas as pd
import yaml
from utils.clustering import JaccardClusters
import numpy as np
import datetime
from utils.cluster_images import GraphPlot
import argparse

parser = argparse.ArgumentParser(
    description = 'Clusters plasmids according to their Jaccard distance on gene content.'
)
parser.add_argument(
    '--load', 
    action = 'store_const',
    const = True, 
    default = False, 
    help='Loads precomputed clusters in the file given by the phylo_config file instead of clustering.'
)
args = parser.parse_args()

if __name__ == '__main__':
    LOAD = args.load
    with open('./configs/clustering_config.yaml', 'r') as config_file:
        config = yaml.load(config_file, Loader=yaml.Loader)
    with open('./configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('./configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)

    #%%
    if not LOAD:
        timestamp = datetime.datetime.now().strftime('%d-%b-%Y__%H-%M-%S')
        write_dir = os.path.join(
                config['paths'][config['similarity']],
                'Clustering_Results_' + timestamp
            )
        os.mkdir(write_dir)

        with open(
            os.path.join(
                write_dir,
                'clustering_parameters.txt'
            ),
            'x'
        ) as outstream:
            outstream.write(
                str(config)
            )
        
        clustering = JaccardClusters(
            config['paths']['plasmids'],
            inflation = config['inflation'],
            reg_factor = config['reg_factor'],
            aggregate = config['aggregate']['apply'],
            min_cluster_size = config['aggregate']['min_cluster_size'],
            min_n_core = config['aggregate']['min_n_core'],
            drop_isolates = True,
            score = config['score'],
            core_threshold = config['core_threshold'],
            resolution = config['resolution'],
            distance_threshold = config['distance_threshold'],
            max_size = config['max_size'],
            linkage = config['linkage']
        )

        affinity = clustering.build_affinity_matrix()
        clusters = clustering.grid_search(
            config['inflation'],
            affinity
        )
        clustering.save(write_dir)
        affinity = clustering.affinity.loc[clusters['Plasmids'].sum()][clusters['Plasmids'].sum()]
        clustering.score()
    else:
        clusters = joblib.load(
            phylo_config['paths']['plasmid-clusters']
        )
        affinity = pd.read_csv(
            '/'+os.path.join(
                *(phylo_config['paths']['plasmid-clusters'].split('/')[:-1]),
                'Affinity.tsv'
            ),
            sep = '\t',
            index_col = 0
        )
        timestamp = 'loaded'
    counts = clusters['Plasmids'].apply(lambda x: len(x))

    print(clusters)
    with plt.style.context('ggplot'):
        f, axs = plt.subplots(
            1,2,
            figsize = (14, 7),
            constrained_layout = True
        )
        color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

        axs[0].boxplot(
            counts.values,
            patch_artist = True,
            boxprops={
                'fill': True,
                'facecolor': color_cycle[0],
                'edgecolor': color_cycle[0]
                },
            medianprops={
                'color': color_cycle[4],
                'linewidth': 3
                },
            capprops={'color': color_cycle[0]},
            flierprops={'markeredgecolor': color_cycle[0]},
            whiskerprops={'color': color_cycle[0]}
        )
        axs[0].set_ylabel('Number of Plasmids')
        axs[0].set_ylim((0, max(counts)+0.05*max(counts)))
        axs[1].hist(
            counts.values,
            bins = np.arange(max(counts)+1)
        )
        axs[1].set_xlabel('Number of Plasmids')
        axs[1].set_ylabel('Number of Clusters')
        axs[1].set_xlim((0, max(counts)))

    plot = GraphPlot(
        joblib.load(config['paths']['plasmids']+'/Plasmids with Clustered Proteins_s9_k5.pkl'),
        clusters,
        affinity,
        data_config['paths']['info_df'],
        data_config['paths']['amr_hits']
    )

    with plt.style.context('ggplot'):
        
        plot.draw_graph(
            config['paths']['plots']+f'/Jaccard-Plasmid-Clustering-{timestamp}.html',
            min_size = 0,
            edge_threshold = .1
        )
        
        plot.plot_stats(figsize = (12, 12), heatmap=True)
