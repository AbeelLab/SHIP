'''
Functions to display annotations and CDS product names.
'''

from typing import Iterable
import numpy as np
import pandas as pd

def join_gene_names(
    before,
    anchor,
    after
):
    return '-'.join(
        [
            np.sort([before, after])[0],
            anchor,
            np.sort([before, after])[1]
            
        ]
    )

def make_3mers(x: list):

    # Defines the direction according to alphabetical
    # order of the genes in the extremeties, without
    # actually looking at strand information
    return [
        join_gene_names(before, anchor, after)
        for after, anchor, before in zip(
            x[-2:]+x[:-2], x, x[2:]+x[:2]
        )
    ]

def get_protein_function(
    proteins: Iterable,
    path_to_representative_names: str
):
    '''
    Returns the protein function of representative proteins (CD-Hit output), given an iterable
    of protein accessions.

    Parameters
    ----------

    proteins: Iterable of str
        Protein accession numbers.
    
    path_to_representative_names: str
        Path to the file containing the representative protein names. Usually in
        DBL/InitialClustering/Representative Proteins/protein-cluster-names-s9-k5.tsv.
    '''
    
    with open(path_to_representative_names, 'r') as instream:
        content = instream.read().split('\n')[2:-1]

    protein_names = [
        ' '.join([first_line, second_line])
        for first_line, second_line in zip(
            content[1::3],
            content[2::3]
        )
    ]
    protein_ids = content[::3]

    names = pd.Series(
        protein_names,
        protein_ids
    )

    common_proteins = set(proteins).intersection(names.index.to_list())
    missing_proteins = set(proteins).difference(names.index.to_list())
    if len(missing_proteins): 
        names = pd.concat(
            [
                names,
                pd.Series(
                    list(missing_proteins),
                    index = list(missing_proteins)
                )
            ]
        )

    return names.loc[list(proteins)]
   
def linebreak(
    s: str,
    max_char: int = 15,
    hyphenate: bool = True
):
    '''
    Adds linebreaks to a string.

    Attributes
    ----------

    max_char: int. Default is 15
        Maximum number of characters allowed by line.
    
    hyphenate: bool. Default is True
        Specifies the behaviour if a word has more than max_char characters. If True,
        the work is broken into two using a hyphen. If False, the word is preserved
        as is, running beyond the line character limit.
    '''
    final_string = ''
    line = ''
    words = s.split()
    for word in words:
        line_len = len(line)
        word_len = len(word)
        if line_len + word_len > max_char:
            # Is the word longer than the character limit?
            if word_len > max_char and hyphenate:
                while word_len > max_char:
                    margin = max_char - line_len - 2
                    line += (' ' + word[:margin] + '-')
                    final_string += ('\n' + line[1:])
                    word = word[margin:]
                    word_len = len(word)
                    line = ''
                    line_len = 0
                line += (' ' + word)
            else:
                # Trigger line break
                final_string += ('\n' + line[1:])
                line = ' ' + word
        else:
            line += (' ' + word)
    final_string += ('\n' + line[1:])
    return final_string[1:]
# %%
