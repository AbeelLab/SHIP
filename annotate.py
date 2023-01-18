'''
Runs Prokka on plasmid FASTA files to transfer the original annotations. 
Outputs the annotations to an Annotations folder. Requires a FASTA file
with all products concatenated.
'''
from annotations.prokka import batch_annotate
import yaml
import sys
import getopt

if __name__ == "__main__":
    with open('./configs/data_config.yaml', 'r') as config_file:
        config = yaml.load(config_file, Loader=yaml.Loader)

    batch_annotate(
        config['paths']['plasmids'],
        config['paths']['proteins'],
        config['paths']['annotations']
    )
