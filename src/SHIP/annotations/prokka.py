import os
from tqdm import tqdm
import subprocess
import warnings
from typing import Union

def run_annotation(
    path: str,
    protein_file: Union[str, None],
    output_dir: str
) -> None:

    if protein_file is not None:
        command = [
            "prokka",
            "--outdir",
            output_dir,
            "--proteins",
            protein_file,
            path
        ]
    else:
        command = [
            "prokka",
            "--outdir",
            output_dir,
            path
        ]  

    process = subprocess.run(command)
    if process.returncode != 0:
        warnings.warn(f'Prokka annotation for fine in {path} exited with code {process.returncode}.')

def batch_annotate(
    path: str,
    protein_file: Union[str, None],
    output_dir: str,
    skip_done: bool = True
) -> None:
    """
    Runs Prokka annotations for each file in path. Assumes a directory structure of:

    path\n
    |\n
    L Species one\n
    |   L_ sequence1.fa\n
    |   L_ sequence2.fa\n
    |\n
    L Species two\n
        L_ sequence3.fa\n
        L_ (...)\n

    Writes the result to output_dir. If included, uses protein_file as argument for 
    --protein when calling Prokka.
    """

    all_files = [
        x[2] 
        for x in os.walk(path)
    ][1:]
    folders = [
        x[0] 
        for x in os.walk(path)
    ][1:]

    if skip_done and os.path.exists(
        output_dir
    ):
        completed_files_list = [
            x[2] 
            for x in os.walk(output_dir)
        ][1:]
        completed_files = []
        for files in completed_files_list:
            completed_files += files
    else:
        completed_files = []

    for files, folder in tqdm(
        zip(
            all_files,
            folders
        ),
        total = len(folders),
        desc = 'Running batch annotations...',
        position = 0,
        colour = 'blue'
    ):
        if not os.path.exists(
            os.path.join(
                output_dir,
                folder.split('/')[-1]
            )
        ):
            os.mkdir(
                os.path.join(
                    output_dir,
                    folder.split('/')[-1]
                )
            )
        for file in tqdm(
            files,
            desc = f'Annotating for species {folder.split("/")[-1]}...',
            total = len(files),
            colour = 'red',
            position=1
        ):
            if file not in completed_files:

                out_dir = os.path.join(
                    output_dir,
                    folder.split('/')[-1],
                    file.split('.')[0]
                )

                if file.endswith('.fa'):
                    run_annotation(
                        os.path.join(
                            folder,
                            file
                        ),
                        protein_file,
                        out_dir
                    )
            
            else:
                print(f'Found annotations for {file}. Skipping...')