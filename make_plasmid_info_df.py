'''
Writes plasmid information in a TSV file.
'''
import yaml
import pandas as pd
from Bio.SeqIO import parse
import numpy as np
import os

if __name__ == '__main__':
    with open('configs/data_config.yaml', 'r') as config_file:
        config = yaml.load(config_file, Loader=yaml.Loader)

    with open(
        os.path.join(
            config['paths']['plasmids_gb'],
            'concat_features.gb'
        )
    ) as instream:
        for n, plasmid in enumerate(parse(instream, 'genbank')):
            species = ' '.join(plasmid.description.split(' ')[:2])
            strain = plasmid.description.split('plasmid')[0].split(' ')[-2]
            plasmid_name = plasmid.description.split('plasmid')[-1].split(',')[0]
            biosample = plasmid.dbxrefs[-1].split(':')[-1]
            id_ = plasmid.id
            size = int(np.round(len(plasmid)/1000))

            if n == 0:
                info_df = pd.DataFrame(
                    [[strain, plasmid_name, biosample, species, size]],
                    columns = ['Strain', 'Plasmid', 'BioSample', 'Organism', 'Size'],
                    index = [id_]
                )
            else: 
                info_df = pd.concat(
                    [
                        info_df,
                        pd.DataFrame(
                            [[strain, plasmid_name, biosample, species, size]],
                            columns = ['Strain', 'Plasmid', 'BioSample', 'Organism', 'Size'],
                            index = [id_]
                        )
                    ]
                )

    info_df.to_csv(
        config['paths']['info_df'],
        sep = '\t'
    )

# %%
