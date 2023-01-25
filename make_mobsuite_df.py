import pandas as pd
import os
import numpy as np
import yaml


if __name__ == '__main__':
    with open('configs/data_config.yaml', 'r') as config_file:
        config = yaml.load(config_file, Loader=yaml.Loader)
    files = list(os.walk(
        config['paths']['mob-types']
    ))[-1][-1]

    files = iter(files)
    file = next(files)
    mob = pd.read_csv(
        file,
        sep = '\t',
        index_col = 0
    )

    for file in files:
        if file.endswith('.fa.txt'):
            mob = pd.concat(
                [
                    mob,
                    pd.read_csv(
                        os.path.join(config['paths']['mob-types'], file),
                        sep = '\t',
                        index_col = 0
                    )
                ]
            )

    mob.index = [x.split('.')[0] for x in mob.index]
    mob.to_csv(
        os.path.join(
            config['paths']['mob-types'],
            'MOBSuite_output_concat.tsv'
        ), 
        sep='\t'
    )

    def concatenate(column: str, mob: pd.DataFrame, name: str = None):
        series = pd.Series(name = name)
        for id_ in np.unique(mob.index):
            types = np.unique(str(mob.loc[id_][column]).split(','))
            series = pd.concat(
                [
                    series,
                    pd.Series(
                        types,
                        name = name,
                        index = [id_] * len(types)
                    )
                ]
            )

        return series

    relaxase = concatenate(
        'relaxase_type(s)',
        mob,
        'Relaxases'
    )
    relaxase.to_csv(
        os.path.join(
            config['paths']['mob-types'],
            'MOBSuite_relaxases.tsv'
        ),
        sep = '\t'
    )

    reps = concatenate(
        'rep_type(s)',
        mob,
        'Rep Types'
    )
    reps.to_csv(
        os.path.join(
            config['paths']['mob-types'],
            'MOBSuite_rep.tsv'
        ),
        sep = '\t'
    )

    mpfs = concatenate(
        'mpf_type',
        mob,
        'MPF Types'
    )
    mpfs.to_csv(
        os.path.join(
            config['paths']['mob-types'],
            'MOBSuite_mpfs.tsv'
        ),
        sep = '\t'
    )

    orit = concatenate(
        'orit_type(s)',
        mob,
        'Orit Types'
    )
    orit.to_csv(
        os.path.join(
            config['paths']['mob-types'],
            'MOBSuite_orit.tsv'
        ),
        sep = '\t'
    )

    orit = concatenate(
        'orit_type(s)',
        mob,
        'Orit Types'
    )
    orit.to_csv(
        os.path.join(
            config['paths']['mob-types'],
            'MOBSuite_orit.tsv'
        ),
        sep = '\t'
    )

    mobility = concatenate(
        'predicted_mobility',
        mob,
        'Predicted Mobility'
    )
    mobility.to_csv(
        os.path.join(
            config['paths']['mob-types'],
            'MOBSuite_mobility.tsv'
        ),
        sep = '\t'
    )
    # %%
