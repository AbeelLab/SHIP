'''
Extracts the regions described in the paper (recombination in E. faecalis and complex class 1 integron
in E. coli and K. pneumoniae) into FASTA files. Performs alignment with Biopython's pairwise2 module.
'''
from utils.files import find_annotation_paths
import yaml
import os
import pandas as pd
import numpy as np
from BCBio.GFF import parse
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SeqRecord import SeqRecord
import argparse

parser = argparse.ArgumentParser(
    description = '''
Extracts the regions described in the paper (recombination in E. faecalis and complex class 1 integron
in E. coli and K. pneumoniae) into FASTA files. Performs alignment with Biopython's pairwise2 module.
'''
)
parser.add_argument(
    '--region', 
    default = 'integron', 
    nargs = 1,
    help = 'Conserved region to analyse. Either "integron" (default) or "recombination".'
)
parser.add_argument(
    '--out', 
    default = './', 
    nargs = 1,
    help = 'Output directory.'
)
args = parser.parse_args() 

if __name__ == '__main__':
    with open('./configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('./configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)

    for k in phylo_config['output-paths']:
        phylo_config['output-paths'][k] = os.path.join(
            phylo_config['results-dir'],
            phylo_config['output-paths'][k]
        )

    def get_alignment_score(
        seq1, seq2
    ):
        alignment = pairwise2.align.globalms(
            seq1,
            seq2,
            match = 1,
            mismatch = 0,
            open = 0,
            extend = 0,
            score_only = True
        )
        return alignment/(
            np.minimum(
                len(seq1), len(seq2)
            )
        )

    assert args.region[0] in ['integron', 'recombination'], f'--region must be either "integron" or "recombination". Got {args.region}.'

    description = args.region[0] #'integron' or 'recombination'

    if description == 'integron':
        ids = ['NZ_CP101516', 'NZ_CP102878', 'NZ_CP102884', 'NZ_CP103504']
        start, end = 'KEBDPHEG_00078', 'HFMAGJGA_00569'
    elif description == 'recombination':
        ids = ['NZ_CP053182', 'NZ_CP068250', 'NZ_CP098027', 'NZ_CP098420']
        start, end = 'CBGDDEJM_00007', 'CBGDDEJM_00013' 

    representatives = pd.read_csv(phylo_config['paths']['representative-proteins'], sep='\t', index_col = 0)['Representative']
    start = representatives[representatives==start].index.to_numpy()
    end = representatives[representatives==end].index.to_numpy()

    paths = find_annotation_paths(
        ids,
        data_config['paths']['annotations'],
        format = '.gff'
    )

    limits = {k:None for k in ids}
    sequences = {k:None for k in ids}
    for path, id in zip(paths, ids):
        found_start, found_end = False, False
        limits_it = []
        reverse_limits = []
        with open(path, 'r') as instream:
            features = next(parse(instream)).features
        start_path = False
        end_path = False
        for feature in features:
            append = False

            if np.any(start == feature.id) and not found_start:
                append = True
                found_start = True
            if np.any(end == feature.id) and not found_end:
                append = True
                found_end = True
            if append:
                if len(limits_it) == 0:
                    limits_it.append(feature.location.start.position)
                    reverse_limits.append(feature.location.end.position)
                else:
                    limits_it.append(feature.location.end.position)
                    reverse_limits.append(feature.location.start.position)
                    break
            
        limits[id] = limits_it
        # Get complete FASTA file
        with open(
            find_annotation_paths(
                [id], 
                data_config['paths']['annotations'],
                format = '.fna'
            )[0],
            'r'
        ) as instream:
            sequences[id] = next(SeqIO.parse(instream, 'fasta')).seq[limits_it[0]:limits_it[1]+1]

    lengths = {
        k: len(v) for k, v in sequences.items()
    }
    # Save sequences as FASTA
    with open(
        os.path.join(
            args.out[0],
            description + '.fa'
        ),
        'w'
    ) as outfile:
        SeqIO.write(
            [SeqRecord(seq, id, name=id, description = description) for id, seq in sequences.items()],
            outfile,
            'fasta'
        )

    score_matrix = pd.DataFrame(
        np.zeros(
            (len(sequences), len(sequences))
        ),
        index = sequences.keys(),
        columns = sequences.keys()
    )
    for n, k in enumerate(sequences.keys()):
        for m, l in enumerate(list(sequences.keys())[n:]):
            score_matrix[l][k] = get_alignment_score(
                sequences[k], sequences[l]
            )
            score_matrix[k][l] = score_matrix[l][k]

    print(score_matrix.round(4)*100)