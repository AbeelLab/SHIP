# Utility functions to process FASTA files
import os
import numpy as np
import pandas as pd
from warnings import warn

def split_fastas(
    path_to_fastas: str,
    output_path: str
) -> None:

    original_files = [
        filename
        for filename in [
                x[2] 
                for x in os.walk(path_to_fastas)
            ][0]
        if not filename.endswith('Identifier')
    ]

    for filename in original_files:

        if not os.path.exists(
            os.path.join(
                output_path,
                '_'.join(filename.split('_')[:2])
            )
        ):
            os.mkdir(
                os.path.join(
                    output_path,
                    '_'.join(filename.split('_')[:2])
                )
            )

        with open(
            os.path.join(
                path_to_fastas,
                filename
            ),
            'r'
        ) as file:
            annotations = file.read().split('\n>')

        for n, annotation in enumerate(annotations):
            accession = annotation.split()[0].split('.')[0]
            if n == 0: accession = accession[1:]
            else:
                annotation = '>' + annotation
            with open(
                os.path.join(
                    output_path,
                    '_'.join(filename.split('_')[:2]),
                    accession + '.fa'
                ),
                'w'
            ) as write_file:
                write_file.write(
                    annotation
                )

def add_plasmid_tag(
    path_to_fastas: str
)-> None:

    all_files = [
        x[2] 
        for x in os.walk(path_to_fastas)
    ][1:]
    folders = [
        x[0]
        for x in os.walk(path_to_fastas)
    ][1:]

    for files, folder in zip(
        all_files,
        folders
    ):
        for file in files:
            with open(
                os.path.join(
                    folder,
                    file
                ),
                'r'
            ) as stream:
                sequence = stream.read()
                header = sequence.split('\n')[0]
                sequence = '\n'.join(sequence.split('\n')[1:])
                accession = header.split()[0]
                header = ' '.join([accession, '[location=plasmid]'])
            with open(
                os.path.join(
                    folder,
                    file
                ),
                'w'
            ) as stream:
                stream.write('\n'.join([header, sequence]))

def add_tags(
    path_to_fastas: str,
    path_to_genbank: str
) -> None:

    all_files = [
        x[2] 
        for x in os.walk(path_to_fastas)
    ][1:]
    folders = [
        x[0]
        for x in os.walk(path_to_fastas)
    ][1:]

    for files, folder in zip(
        all_files,
        folders
    ):
        with open(
            os.path.join(
                path_to_genbank,
                folder.split('/')[-1] + '_genbank.gb'
            ),
            'r'
        ) as stream:

            # Break info into samples
            all_gb = stream.read().split('LOCUS')
            headers = [
                text.split('\n')[0]
                for text in all_gb[1:]
            ]
            accessions = [
                text.split()[0]
                for text in headers
            ]
            locus_info = []
            for text, accession in zip(
                headers,
                accessions
            ):
                split_header = text.split()
                if 'circular' in split_header:
                    locus_info.append('circular')
                elif 'linear' in split_header:
                    locus_info.append('linear')
                else:
                    raise TypeError(f'No valid topology information found in {accession}.')

        locus_info = np.array(locus_info)
        accessions = np.array(accessions)

        for file in files:
            if file.endswith('.fa'):
                topology = locus_info[accessions == file[:-3]][0]
                with open(
                    os.path.join(
                        folder,
                        file
                    ),
                    'r'
                ) as stream:
                    sequence = stream.read()
                    header = sequence.split('\n')[0]
                    sequence = '\n'.join(sequence.split('\n')[1:])
                    accession = header.split()[0]
                    header = ' '.join([accession, f'[location=plasmid] [topology={topology}]'])
            
                with open(
                    os.path.join(
                        folder,
                        file
                    ),
                    'w'
                ) as stream:
                    stream.write('\n'.join([header, sequence]))

def read_plasmid_info(
    path_to_info: str
) -> pd.DataFrame:

    return pd.read_csv(
        path_to_info,
        sep = '\t',
        index_col = 0
    )

def add_info_tags(
    path_to_info: str,
    path_to_plasmids: str
) -> None:

    all_files = [
        x[2] 
        for x in os.walk(path_to_plasmids)
    ][1:]
    folders = [
        x[0]
        for x in os.walk(path_to_plasmids)
    ][1:]

    info = read_plasmid_info(path_to_info)
    info.index = [x.split('.')[0] for x in info.index]

    for files, folder in zip(
        all_files,
        folders
    ):
        for file in files:

            accession = file.split('.')[0]
            try:
                plasmid_info = info.loc[accession]
            except:
                warn(f'Could not find info for accession {accession}.')
                continue
            finally:

                with open(
                    os.path.join(
                        folder,
                        accession+'.fa'
                    ),
                    'r'
                ) as stream:
                    content = stream.read().split('\n')
                    header = content[0]
                    content = '\n'.join(content[1:])
                    header += f' [organism={plasmid_info["Organism"]}] [strain={plasmid_info["Strain"]}] [plasmid-name={plasmid_info["Plasmid"]}]'
                    
                with open(
                    os.path.join(
                        folder,
                        accession+'.fa'
                    ),
                    'w'
                ) as stream:

                    stream.write(
                        '\n'.join([
                            header,
                            content
                        ])
                    )

                build_submol_yaml(
                    os.path.join(
                        folder,
                        accession+'_submol.yaml'
                    ),
                    species = plasmid_info['Organism'],
                    strain = plasmid_info['Strain']
                )

                build_generic_yaml(
                    os.path.join(
                            folder,
                            accession+'_generic.yaml'
                    ),
                    os.path.join(
                        accession+'_submol.yaml'
                    )
                    ,
                    os.path.join(
                        accession+'.fa'
                    )
                )
                
def build_submol_yaml(
    output_dir: str,
    species: str,
    strain: str
):
    assert output_dir.endswith('.yaml'), 'Invalid file extension. Should be .yaml.'
    with open(
        output_dir,
        'x'
    ) as stream:
        stream.write(
f'''organism:
    genus_species: {species}
    strain: {strain}
'''
        )

def build_generic_yaml(
    output_dir: str,
    path_to_submol: str,
    path_to_fasta: str
):
    assert output_dir.endswith('.yaml'), 'Invalid file extension. Should be .yaml.'
    with open(
        output_dir,
        'x'
    ) as stream:
        stream.write(
f'''fasta:
    class: File
    location: {path_to_fasta}
submol:
    class: File
    location: {path_to_submol}
'''
        )       
# %%
