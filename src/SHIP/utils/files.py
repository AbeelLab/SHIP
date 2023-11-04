'''
Useful classes for file manipulation and navigating the folder structure.
'''

import os
from typing import Iterable
import numpy as np
import warnings
import pandas as pd
from utils.amrfinder import filenames_from_contig_ids

def find_annotation_paths(
    accessions: Iterable,
    annotations_dir: str,
    format: str = '.fa',
):
    '''
    Returns the paths to the annotations of plasmids with accession number 
    specified in the parameter "accessions", in the format given by "format".
    Looks for files in "annotations_dir", a folder containing the Prokka
    output files.

    Parameters
    ---------

    accessions: Iterable of str
        Accession numbers of the plasmids to search for.

    annotations_dir: str
        Path to the folder containing the Prokka annotations.
    
    format: str. Default is .fa
        File extension of the files to retrieve.

    Returns
    -------

    paths: List of str
        List of paths to each plasmid annotation.
    
    '''
    folders = list(os.walk(annotations_dir))[0][1]
    
    paths = []
    for accession in accessions:
        if accession in folders:
            for file in list(os.walk(os.path.join(annotations_dir, accession)))[0][2]:
                if file.endswith(format) and file.startswith('PROKKA'):
                    paths.append(os.path.join(annotations_dir, accession, file))
                    break
        else:
            warnings.warn(f'Could not find {format} annotation file for accession {accession}. Skipping this plasmid...')

    return paths

def concat_files(
    paths: str,
    output: str
) -> str:

    if not os.path.exists(output):
        open_mode = 'x'
    else:
        open_mode = 'w'

    with open(output, open_mode) as out_stream:
        out_stream.write('')

    for path in paths:
        print(f'Writting annotations in {path}')
        with open(path, 'r') as in_stream:
            with open(output, 'a') as out_stream:
                out_stream.write(in_stream.read()+'\n')
    
    return output

def annotation_qc(
    paths: str
):
    '''
    Deletes all annotations with empy FASTA files.
    '''

    keep = []
    delete = []
    for path in paths:
        with open(path, 'r') as instream:
            fasta = instream.read()

        if len(fasta)>0:
            keep.append(path)
        else:
            delete.append(path)

    print(f'Discarding files in {delete}.')
    return keep

def get_amr(
    accessions: Iterable,
    path_to_amr: str,
    annotations_path: str,
    no_amr_key: str = 'Susceptible'
) -> pd.DataFrame:
    '''
    Returns a DataFrame containing the AMR gene information in the AMRFinder
    output table for a set of plasmids given their accession numbers.

    All fields for susceptible plasmids (no AMR genes found) are set to the
    string passed to no_amr_key. A plasmid is assumed to be susceptible if
    it is not in the index of the AMRFinder table.
    '''

    amr_df = pd.read_csv(path_to_amr, sep = '\t', index_col = 0)
    contig_ids = filenames_from_contig_ids(annotations_path)
    amr_df.index = [contig_ids[x] for x in amr_df.index]


    plasmids_in_amr_df = amr_df.index.to_numpy()
    common_plasmids = np.array(
        set(plasmids_in_amr_df).intersection(accessions)
    )
    susceptible_plasmids = np.array(
        list(set(accessions).difference(plasmids_in_amr_df))
    )

    amr_df = amr_df.loc[common_plasmids]
    susceptible_df = pd.DataFrame(
        [[no_amr_key]*amr_df.shape[1]]*len(susceptible_plasmids),
        index = susceptible_plasmids,
        columns = amr_df.columns
    )

    return pd.concat(
        [amr_df, susceptible_df]
    )
