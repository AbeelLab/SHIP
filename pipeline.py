import argparse
import os
import subprocess
from Bio.SeqIO import parse, write
import warnings
import glob
import shutil

parser = argparse.ArgumentParser(
    description = '''
Apply SHIP to a set of contigs resulting from plasmid assembly with SCAPP.
'''
)
parser.add_argument(
    '--input', '-i', nargs = 1, required = True,
    help = 'File containing the output from SCAPP. Should be named final.contigs.confident_cycs.fasta.'
)
parser.add_argument(
    '--amrfinder', '-a', action='store_true',
    help='If set, runs AMRFinder+ on the plasmid sequences.'
)
parser.add_argument(
    '--prokka', '-p', action='store_true',
    help='If set, runs Prokka on the plasmid sequences.'
)
parser.add_argument(
    '--blast-db', '-b', nargs = 1, required = False,
    help = 'FASTA or GenBank file containing products to be used as references in the Prokka annotation.'
)
parser.add_argument(
    '--cdhit', '-c', action='store_true',
    help='If set, runs CD-HIT on the plasmid sequence features.'
)
parser.add_argument(
    '--ship', '-s', action='store_true',
    help='If set, runs SHIP on the plasmid sequences.'
)
parser.add_argument(
    '--out', '-o', nargs = 1, required = False,
    help = 'Output directory.'
)

def mkdir(path:str):
    if not os.path.exists(path): os.mkdir(path)

SHIP_PATH = os.path.join('ship.py')

if __name__=='__main__':

    args = parser.parse_args()
    args.out = args.out[0]
    try: args.blast_db = args.blast_db[0]
    except: pass
    args.input = args.input[0]

    # Make out dir
    tmp_path = os.path.join(args.out, 'tmp')
    mkdir(args.out)
    mkdir(tmp_path)

    amrfinder_dir = os.path.join(args.out, 'amrfinder_output')
    mkdir(amrfinder_dir)
    if args.amrfinder:
        print('Running AMRFinder+.')
        subprocess.run(['amrfinder', '-n', args.input, '-o', os.path.join(amrfinder_dir, 'amrfinder_output.txt'),
                        '--log', os.path.join(amrfinder_dir, 'amrfinder.log')])

    # Split FASTA file into one file per plasmid
    tmp_fasta_dir = os.path.join(tmp_path, 'plasmid_fasta')
    mkdir(tmp_fasta_dir)
    with open(args.input) as instream:
        content = parse(instream, 'fasta')
        for plasmid in content.records:
            with open(os.path.join(tmp_fasta_dir, f'{plasmid.id}.fa'), 'w') as outstream:
                write(plasmid, outstream, 'fasta')
    
    annotations_dir = os.path.join(args.out, 'prokka_output')
    mkdir(annotations_dir)

    if args.prokka:
        print('Running Prokka.')
        if not hasattr(args, 'blast_db'):
            warnings.warn('No BLAST database was passed. Resorting to the base Prokka database. Annotations are unlikely to be correct or informative.')
            for plasmid in glob.glob(f'{tmp_fasta_dir}/*.fa'):
                subprocess.run(['prokka', '--outdir', annotations_dir, plasmid])
        else:
            for plasmid in glob.glob(f'{tmp_fasta_dir}/*.fa'):
                plasmid_id = plasmid.split('/')[-1][:-3].replace('.', '_')
                subprocess.run(['prokka', '--outdir', os.path.join(annotations_dir, plasmid_id), '--proteins', args.blast_db, plasmid])

    cdhit_dir = os.path.join(args.out, 'cdhit_output')
    mkdir(cdhit_dir)
    if args.cdhit:
        print('Running CD-HIT.')
        with open(os.path.join(tmp_path, 'concat_cds.faa'), 'w') as file: file.write('')
        for plasmid in glob.glob(f'{annotations_dir}/*/*.faa'):
            content = parse(plasmid, 'fasta')
            with open(os.path.join(tmp_path, 'concat_cds.faa'), 'a') as outstream:
                write(content, outstream, 'fasta')
        subprocess.run(['cdhit', '-i', os.path.join(tmp_path, 'concat_cds.faa'), '-c', '0.9', '-n', '5', '-o', os.path.join(cdhit_dir, 'cdhit_output')])

    if args.ship:
        print('Running SHIP.')
        ship_dir = os.path.join(args.out, 'ship_output')
        mkdir(ship_dir)
        subprocess.run(['python', SHIP_PATH, '-a', annotations_dir, '-o', ship_dir, '-c', 
            os.path.join(cdhit_dir, 'cdhit_output.clstr'), '-r', amrfinder_dir, '-m', '0.1',
            '-l', '5', '-L', '9', '-n', '3'])

    shutil.rmtree(tmp_path)
    print('Done!')