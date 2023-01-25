
'''
Calls AMRFinder+ to find AMR genes in all plasmids in 
data_config/paths/amrfinder_output
'''

import subprocess
import os
import yaml

if __name__ == "__main__":
    with open('configs/data_config.yaml', 'r') as config_file:
        config = yaml.load(config_file, Loader=yaml.Loader)

    files = []
    for x in list(
        os.walk(config['paths']['concat-fasta'])
    )[0][-1]:
        print(x)
        if x.endswith('.fa'):
            subprocess.run(
                [
                    'amrfinder',
                    '-n',
                    os.path.join(config['paths']['concat-fasta'], x),
                    '-o',
                    os.path.join(
                        config['paths']['amrfinder_output'],
                        x.split('.')[0]+'.txt'
                    ),
                    '--log',
                    os.path.join(
                        config['paths']['amrfinder_output'],
                        'Logs',
                        'Logfile_'+x.split('.')[0]+'.err'
                    )
                ]
            )
