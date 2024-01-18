# Functions to transform the AMRFinder+ output into files accepted by SHIP

import pandas as pd
import os
from typing import Iterable
import glob
import logging

def make_amrfinder_df(
    amrfinder_out: str,
    out: str
) -> pd.DataFrame:
    '''
    Builds a DataFrame containing the AMRFinder+ gene hits.
    '''
    logging.info(f"Reading AMRFinderPlus output in {amrfinder_out}")

    amrfinder_df = pd.read_table(amrfinder_out)
    amrfinder_df = amrfinder_df.drop(['Protein identifier', 'HMM id', 'HMM description'], axis=1).set_index(['Contig id', 'Gene symbol'])

    # Save to tsv
    logging.debug(f"Saving a CSV of AMRFinderPlus results in {out}")
    amrfinder_df.to_csv(out, sep = '\t')

    return amrfinder_df

def filenames_from_contig_ids(annotations_path: str) -> pd.Series:
    
    fna_files = glob.glob(f'{annotations_path}/*/*.fna')
    plasmid_names, contig_ids = [], []
    for path in fna_files:
        plasmid_names.append(path.split('/')[-2])
        with open(path) as instream: contig_ids.append(instream.readline().split(" ")[0][1:].split("\n")[0])

    return pd.Series(plasmid_names, index = contig_ids, name = 'contig_ids')
