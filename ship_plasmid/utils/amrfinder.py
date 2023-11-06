# Functions to transform the AMRFinder+ output into files accepted by SHIP

import pandas as pd
import os
from typing import Iterable
import glob

def make_amrfinder_df(
    amrfinder_out: str,
    out: str
) -> pd.DataFrame:
    '''
    Builds a DataFrame containing the AMRFinder+ gene hits.
    '''
    print(amrfinder_out)
    print(list(os.walk(amrfinder_out)))

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

def filenames_from_contig_ids(
    annotations_path: str
) -> pd.Series:
    
    fna_files = glob.glob(f'{annotations_path}/*/*.fna')
    plasmid_names, contig_ids = [], []
    for path in fna_files:
        plasmid_names.append(path.split('/')[-2])
        with open(path) as instream:
            contig_ids.append(instream.readline()[1:-1])

    return pd.Series(plasmid_names, index = contig_ids, name = 'contig_ids')
