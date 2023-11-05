'''
Processes the CD-HIT output.
Obtains the plasmid gene sequence with homolog cluster representatives.
'''
#%%
import os
import pandas as pd
import yaml
from ship_plasmid.utils.proteins import ProteinClusters
import joblib

def process_cdhit(
    cdhit_out: str,
    annotations_path: str,
    out: str
):

    prot_clust = ProteinClusters(path_clustering = cdhit_out)
    
    prot_clust.to_tsv(os.path.join(out, 'protein_cluster_membership'))
    joblib.dump(prot_clust.clusters, os.path.join(out, 'protein_cluster_membership.pkl'),compress = 3)
    prot_clust.build_protein_db(os.path.join(out, 'Protein Clusters.tsv'))

    # Find all annotated plasmids
    accessions = []
    for x in list(os.walk(annotations_path))[1:]:
        for subfile in x[-1]:
            if subfile.endswith('.faa') and subfile.startswith('PROKKA'):
                accessions.append(os.path.join(x[0], subfile))

    raw_proteins = prot_clust.get_raw_proteins(accessions)
    clustered_plasmids = prot_clust.replace_with_cluster(raw_proteins)
    prot_clust.clustered_proteins_to_tsv(os.path.join(out, 'Plasmids with Clustered Proteins'))

# %%
