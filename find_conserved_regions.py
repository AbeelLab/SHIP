'''
This script enables the manual search of conserved AMR regions in plasmids.
'''
#%%
import os
import yaml
import joblib
import numpy as np
from utils.files import get_amr
from utils.motifs import MotifFinder
import argparse

parser = argparse.ArgumentParser(
    description = '''
Manual search for AMR regions in plasmids. Requires a pre-computed phylogeny.
'''
)
parser.add_argument(
    '--cluster',
    default = None,
    nargs = 1,
    help = 'Jaccard cluster in which to find AMR genes. If not specified, assumes that there is only one cluster.'
)
args = parser.parse_args() 

if __name__ == '__main__':
    with open('./configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('./configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)

    for k in phylo_config['output-paths']:
        phylo_config['output-paths'][k] = os.path.join(
            phylo_config['results-dir'],
            phylo_config['output-paths'][k]
        )

    cluster = args.cluster[0]
    phylos = joblib.load(phylo_config['output-paths']['phylogenies'])

    if cluster is not None:
        phylo = phylos[int(cluster)]
    else:
        phylo = phylos

    plasmid_ids = phylo.accessions
    amr_df = get_amr(
        plasmid_ids,
        data_config['paths']['amr_hits']
    )

    resistant_ids = np.unique(amr_df[amr_df['Gene symbol']!= 'Susceptible'].index)

    motifs = MotifFinder(
        phylo,
        data_config,
        phylo_config,
        min_n_edges = 2,
        graph_max_distance = .5
    )
    motifs.start(resistant_ids)

# %%
