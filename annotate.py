'''
Runs Prokka on plasmid FASTA files to transfer the original annotations. 
Outputs the annotations to an Annotations folder. Requires a FASTA file
with all products concatenated.
'''
from annotations.prokka import batch_annotate
from annotations.databases import cat_genbank
import yaml
import argparse
import os
import subprocess

parser = argparse.ArgumentParser(
    description = 'Annotate plasmid FASTA nucleotide sequences in the path specified in data_config/paths/plasmids, using Prokka.'
)
parser.add_argument(
    '--build', 
    action = 'store_const',
    const = True, 
    default = False, 
    help='Builds a protein FASTA database with all plasmid features. Must be set at least once.'
)
args = parser.parse_args()

if __name__ == "__main__":
    with open('./configs/data_config.yaml', 'r') as config_file:
        config = yaml.load(config_file, Loader=yaml.Loader)

    if args.build:
        cat_genbank(
            config['paths']['plasmids_gb'],
            os.path.join(
                config['paths']['plasmids_gb'],
                'concat_features.gb'
            ),
            format='gb',
            bash = False
        )
        process = subprocess.run(
            [
                'makeblastdb',
                '-dbtype',
                'prot',
                '-in',
                os.path.join(
                    config['paths']['plasmids_gb'],
                    'concat_features.gb'
                ),
                '-out',
                config['paths']['proteins']
            ]
        )
        assert process.returncode != 0, 'Failed to execute makeblastdb.'


    batch_annotate(
        config['paths']['plasmids'],
        config['paths']['proteins'],
        config['paths']['annotations']
    )
