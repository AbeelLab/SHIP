#%%
'''
Obtains the plasmid gene sequence with homolog cluster representatives.
'''
import os
import pandas as pd
import yaml
from utils.proteins import ProteinClusters
from utils.files import find_annotation_paths

if __name__ == '__main__':

    with open('configs/clustering_config.yaml', 'r') as config_file:
        config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)

    prot_clust = ProteinClusters(
        similarity = 0.9,
        word_length = 5,
        path_clustering = config['paths']['cdhit']
    )
    prot_clust.build_protein_db(
        os.path.join(
            config['paths']['membership'],
            'Protein Clusters.pkl'
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

    prot_clust.replace_with_cluster(
        prot_clust.get_raw_proteins(accessions)
    )
    prot_clust.clustered_proteins_to_tsv(
        os.path.join(
            config['paths']['plasmids'],
            'Plasmids with Clustered Proteins_s9_k5'
        )
    )


# %%
