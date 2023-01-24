
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
        os.walk(config['paths']['plasmids'])
    )[1:]:
        for xx in x[-1]:
            if xx.endswith('.fa') or xx.endswith('.fasta') or xx.endswith('.faa'):
                subprocess.run(
                    [
                        'amrfinder',
                        '-n',
                        os.path.join(x[0], xx),
                        '-o',
                        os.path.join(
                            config['paths']['amrfinder_output'],
                            xx.split('.')[0]+'.txt'
                        ),
                        '--log',
                        os.path.join(
                            config['paths']['amrfinder_output'],
                            'Logs',
                            'Logfile_'+xx.split('.')[0]+'.err'
                        )
                    ]
                )
    
