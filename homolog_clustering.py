'''
Clusters homologs with CD-HIT.
'''
import os
import yaml
import subprocess
from annotations.databases import cat_genbank

if __name__ == '__main__':
    with open('configs/clustering_config.yaml', 'r') as config_file:
        config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)

    all_fasta = []
    for full_dir, _, files in os.walk(data_config['paths']['annotations']):
        for file in files:
            if file.endswith('.faa'):
                with open(os.path.join(full_dir, file)) as instream:
                    all_fasta.append(
                        instream.read()
                    )
    with open(
        os.path.join(
            config['paths']['cdhit'],
            'tmp-concat-fasta.fa'
        ),
        'w'
    ) as outstream:
        outstream.write(
            '\n'.join(all_fasta)
        )
    subp = subprocess.run(
        [
            'cdhit',
            '-i',
            os.path.join(
                config['paths']['cdhit'],
                'tmp-concat-fasta.fa'
            ),
            '-c',
            '0.9',
            '-n',
            '5',
            '-o',
            os.path.join(
                config['paths']['cdhit'],
                'protein-clustering-s9-k5'
            )
        ]
    )

    os.remove(
        os.path.join(
            config['paths']['cdhit'],
            'tmp-concat-fasta.fa'
        )
    )
