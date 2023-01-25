'''
Computes statitics for Jaccard-based clusters.
'''
#%%
import os
import joblib
import pandas as pd
import yaml
import numpy as np
from scipy.stats import fisher_exact
from utils.plasmid_typing import LabeledNetwork
#%%

if __name__ == '__main__':
    with open('configs/clustering_config.yaml', 'r') as config_file:
        config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)

    clusters = joblib.load(
            phylo_config['paths']['plasmid-clusters']
        )
    counts = clusters['Plasmids'].apply(lambda x: len(x))
    info = pd.read_csv(
        data_config['paths']['info_df'],
        sep = '\t',
        index_col = 0
    )
    info.index = [x.split('.')[0] for x in info.index]
    clusters = clusters['Plasmids']
    # %% Fisher's exact test for species inside a cluster vs species outside cluster
    '''
    Contingency table for each species X and cluster Y:

    No. species X     | No. species not X
    in cluster Y      | in cluster Y
    ------------------|-------------------
    No. species X     | No. species not X
    outside cluster Y | outside cluster Y

    '''

    for species in np.unique(info['Organism']):
        for cluster in clusters.index:

            # Info slice for the query cluster
            cluster_info = info.loc[
                clusters[cluster]
            ]
            n_species_in = len(
                cluster_info[cluster_info['Organism']==species]
            )
            if n_species_in < 20:
                continue
            n_not_species_in = counts[cluster] - n_species_in

            n_species_out = len(
                [
                    x
                    for x in info.index
                    if x not in clusters[cluster]
                    and info.loc[x]['Organism'] == species
                ]
            )

            n_not_species_out = len(
                [
                    x
                    for x in info.index
                    if x not in clusters[cluster]
                    and info.loc[x]['Organism'] != species
                ]
            )

            table = np.array([
                [n_species_in, n_not_species_in],
                [n_species_out, n_not_species_out]
            ])

            _, p = fisher_exact(table)

            if p < 0.95:
                print(table)

                print(f'Cluster {cluster} has a unique distribution for species {species} (p={p:.5}).')
                print(f'{np.round(table[0,0]/counts[cluster]*100)}% of the plasmids in cluster {cluster} are from {species}.')
                print(f'{np.round(table[0,0]/(table[1,0]+table[0,0])*100, 0)}% of the {species} plasmids are in cluster {cluster}.\n')
    sequences = joblib.load(
        os.path.join(
            config['paths']['plasmids'],
            'Plasmids with Clustered Proteins_s9_k5.pkl'
        )
    )['Proteins']
    for cluster in phylo_config['clusters']:
        cluster_seqs = sequences.loc[clusters[cluster]]
        all_genes = np.unique(cluster_seqs.sum())
        n_plasmids = len(cluster_seqs)
        
        ratio_presence = np.array(
            [
                len([x for x in cluster_seqs if gene in x])/n_plasmids
                for gene in all_genes
            ]
        )

        n_core = sum(ratio_presence > .5)

        print(f'Number of core genes in cluster {cluster}: {n_core}')
    # %% Plot Network with Species Labels
    net = LabeledNetwork(
        distance_threshold = .5
    )
    affinity = pd.read_csv(
        os.path.join(
            config['paths']['jaccard'],
            'Affinity.tsv'
        ),
        sep = '\t',
        index_col = 0
    )
    affinity = affinity.loc[np.sum(clusters)][np.sum(clusters)]
    # Get common indices
    common_idx = list(set(affinity.index).intersection(info.index))
    net.from_distance(
        1-affinity.loc[common_idx][common_idx],
        info['Organism'].loc[common_idx]
    )
    net.show()
# %%
