'''
Obtains the plasmid gene sequence with homolog cluster representatives.
'''
#%%
import os
import pandas as pd
import yaml
from utils.proteins import ProteinClusters
import joblib

if __name__ == '__main__':

    with open('configs/clustering_config.yaml', 'r') as config_file:
        config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)

    prot_clust = ProteinClusters(
        similarity = 0.9,
        word_length = 5,
        path_clustering = config['paths']['cdhit']
    )
    prot_clust.to_tsv(
        os.path.join(
            config['paths']['membership'],
            'protein_cluster_membership'
        )
    )
    
    joblib.dump(
        prot_clust.clusters,
        os.path.join(
            config['paths']['membership'],
            'protein_cluster_membership_s9_k5.pkl'
        ),
        compress = 3
    )

    prot_clust.build_protein_db(
        os.path.join(
            config['paths']['membership'],
            'Protein Clusters.tsv'
        )
    )
    # Find all annotated plasmids
    accessions = []
    for x in list(os.walk(data_config['paths']['annotations']))[1:]:
        for subfile in x[-1]:
            if subfile.endswith('.faa') and subfile.startswith('PROKKA'):
                accessions.append(
                    os.path.join(
                        x[0], subfile
                    )
                )

    #prot_clust.clusters = prot_clust.clusters.drop_duplicates(keep='first')
    raw_proteins = prot_clust.get_raw_proteins(accessions)
    clustered_plasmids = prot_clust.replace_with_cluster(
        raw_proteins
    )
    prot_clust.clustered_proteins_to_tsv(
        os.path.join(
            config['paths']['plasmids'],
            'Plasmids with Clustered Proteins'
        )
    )

# %%
