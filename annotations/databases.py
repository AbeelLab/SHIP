#%%
import os
from typing import Iterable
import warnings
from Bio import SeqIO
import subprocess


def cat_genbank(
    path: str,
    out_path: str,
    format: str = 'gbff',
    bash = True
) -> None:
    '''
    Concatenates (using bash commands or not) all files of type 'format' inside
    'path'. Stores the result in 'out_path'.
    '''
    all_files = [
        x[2] 
        for x in os.walk(path)
    ][1:]

    subdirs = [
        x[1] 
        for x in os.walk(path)
    ][0]
    
    files = []
    if bash:
        for dir_files, dirs in zip(all_files, subdirs):
            files += [
                os.path.join(
                    path,
                    dirs,
                    x
                )
                for x in dir_files if x.endswith('.'+format)
            ]

        command = ['cat'] + files + ['>'] + [out_path] + ['|'] + ['less']
        process = subprocess.run(command)
        if process.returncode != 0:
            warnings.warn(f'Subprocess exited with exit code {process.returncode}.')
    else:
        if not os.path.exists(out_path):
            with(open(out_path, 'x')) as out_stream:
                out_stream.write('')
        for dir_files, dirs in zip(all_files, subdirs):
            for x in dir_files:
                if x.endswith('.'+format):
                    with open(
                        os.path.join(
                            path,
                            dirs,
                            x
                        ),
                        'r'
                    ) as gb_stream:
                        to_write = gb_stream.read()
                    with open(out_path, 'a') as out_stream:
                        out_stream.write(to_write+'\n')


def concat_genbank_files(
    path: str,
    output: str
) -> None:
    '''
    Concatenates all .gb files in 'path' into a single file.
    Writes that file into the path specified in 'output'.
    '''

    all_files = [
        x[2] 
        for x in os.walk(path)
    ][0]

    all_files = [os.path.x for x in all_files if x.endswith('.gb')]

    content = []
    for file in all_files:
        file = os.path.join(path, file)
        with open(file, 'r') as stream:
            content.append(stream.read())
    
    with open(output, 'x') as stream:
        stream.write(
            '\n'.join(content)
        )

def genbank_to_embl(
    path: str,
    output_dir: str
) -> None:

    bio_convert(
        path,
        output_dir,
        'genbank',
        'embl'
    )

def bio_convert(
    path: str,
    output_dir: str,
    in_fmt: str,
    out_fmt: str
) -> None:
    content = SeqIO.parse(path, in_fmt)
    SeqIO.write(content, output_dir, out_fmt)

def fasta_to_genbank(
    path: str,
    output_dir: str
):
    content = SeqIO.parse(path, 'fasta')
    for sequence in content:
        sequence.annotations['molecule_type'] = 'DNA'
    with open(output_dir, 'w') as stream:
        SeqIO.write(content, stream, 'genbank')

def concat_protein_annotations(
    annotations_path: str,
    paths: Iterable,
    out_dir: str
) -> None:

    for n, folder in enumerate(paths):
        files = [
            x[2] 
            for x in os.walk(
                os.path.join(
                    annotations_path,
                    folder
                )
            )
        ][0]
        file = [x for x in files if x.endswith('.faa')][0]
        with open(
            os.path.join(
                annotations_path, 
                folder, 
                file
            ), 
            'r'
        ) as stream:
            content = stream.read()
        if n == 0:
            open_type = 'x'
        else:
            open_type = 'a'
        with open(out_dir, open_type) as stream:
            stream.write(content+'\n')
# %%
