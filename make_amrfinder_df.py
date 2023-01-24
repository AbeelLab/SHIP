'''
Builds a DataFrame containing the AMRFinder+ gene hits
and the information for each plasmid, contained in 
./Data/GenBank/plasmid_info.csv
'''
import pandas as pd
import os
import yaml

if __name__ == 'main':
    with open('configs/data_config.yaml', 'r') as config_file:
        config = yaml.load(config_file, Loader=yaml.Loader)

    # Load info file
    info = pd.read_csv(
        config['paths']['info_df'],
        sep = '\t',
        index_col=0
    )
    files = [x for x in os.walk(config['paths']['amrfinder_output'])][0][2]

    amr_df = pd.concat(
        [
            pd.read_csv(
                os.path.join(
                    config['paths']['amrfinder_output'],
                    filename
                ),
                sep = '\t'
            )
            for filename in files
        ]
    ).drop(
        [
            'Protein identifier',
            'HMM id',
            'HMM description'
        ],
        axis = 1
    ).set_index(
        [
            'Contig id',
            'Gene symbol'
        ]
    )

    # Drop indices not in the info DataFrame
    amr_df.drop(
        set(
            amr_df.index.get_level_values('Contig id')
        ).difference(info.index),
        inplace = True
    )

    # Join the two dataframes. This is weird because amr_df has a pair of
    # primary keys
    for column in info:
        amr_df[column] = info[column].loc[
            amr_df.index.get_level_values('Contig id')
        ].values

    # Save to tsv
    amr_df.to_csv(
        config['paths']['amr_hits'],
        sep = '\t'
    )
