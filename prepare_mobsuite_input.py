import os
import yaml
import pandas as pd
from utils.files import find_annotation_paths
import shutil
import argparse

parser = argparse.ArgumentParser(
    description = '''
Copies the annotated plasmid files inside one or more clusters into one folder, 
ready to be used as input for MOBFinder. Uses the files in data/Plasmid Networks
and includes all plasmids in the Jaccard-Subclusters.tsv file.
'''
)

parser.add_argument(
    '--paper',
    action = 'store_const',
    const = True, 
    default = False, 
    help = 'Loads the Jaccard-based clusters used in the original paper.'
)

args = parser.parse_args() 

if __name__ == '__main__':
    with open('configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/clustering_config.yaml', 'r') as config_file:
        clustering_config = yaml.load(config_file, Loader=yaml.Loader)

    # Retrieve the plasmid IDs of the plasmids in the mega clusters
    # Create a new folder with the FASTA files for those plasmids
    if args.paper:
        phylo_config['results-dir'] = os.path.join(
            phylo_config['results-dir'],
            'Paper'
        )
    clusters = pd.read_csv(
        os.path.join(
            phylo_config['results-dir'],
            'Jaccard-Subclusters.tsv'
        ),
        sep = '\t',
        index_col = 0
    )

    ids = clusters.index.to_list()
    paths = find_annotation_paths(ids, data_config['paths']['annotations'], format='.fna')

    for path in paths:
        accession = path.split('/')[-2]+'.fa'
        shutil.copyfile(
            path, 
            os.path.join(
                data_config['paths']['mob-suite'],
                accession
            )
        )

# %%
