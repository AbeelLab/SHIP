# Functions to transform the AMRFinder+ output into files accepted by SHIP

import pandas as pd
import os

def make_amrfinder_df(
    amrfinder_out: str,
    out: str
) -> pd.DataFrame:
    '''
    Builds a DataFrame containing the AMRFinder+ gene hits.
    '''

    files = list(os.walk(amrfinder_out))[0][2]

    individual_df = []
    for filename in files:
        if filename.endswith('.txt'):
            try:
                individual_df.append(pd.read_csv(os.path.join(amrfinder_out, filename), sep = '\t'))
            except pd.errors.EmptyDataError: 
                # No AMR genes found
                pass  
    amr_df = pd.concat(individual_df).drop(['Protein identifier', 'HMM id', 'HMM description'], axis = 1).set_index(
        ['Contig id', 'Gene symbol'])

    # Save to tsv
    amr_df.to_csv(out, sep = '\t')

    return amr_df