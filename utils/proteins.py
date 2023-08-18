'''
Class to process CD-HIT outputs and build homolog families.
'''
#%%
import os
import warnings
import pandas as pd
import yaml
from typing import Tuple, Union, Iterable
import numpy as np
import joblib
from copy import deepcopy

class ProteinClusters:

    def __init__(
        self,
        path_clustering: str
    ):

        self.path_clstr = path_clustering
        with open(self.path_clstr, 'r') as stream: self.__cluster_txt = stream.read()

        self.clusters = self.__parse_clstr_file(self.__cluster_txt)

    def __find_representative(
        self,
        lines: str
    ) -> Union[str, None]:
        for line in lines:
            if len(line.split('\t')[1].split()) == 3:
                return self.__get_name(line)
        return None

    def __get_name(
        self,
        line: str
    ) -> str:
        return line.split('\t')[1].split()[1][1:][:-3]

    def __get_protein_similarity(
        self,
        line: str
    ) -> float:
        if line.split('\t')[1].split()[2] == '*':
            return float(100)
        else:
            return float(line.split('\t')[1].split()[3][:-1])
        
    def __parse_cluster(
        self,
        txt: str
    ) -> pd.DataFrame:
        lines = txt.split('\n')
        cluster_number = lines[0].split()[-1]
        lines = lines[1:]
        representative = self.__find_representative(lines)

        for n, line in enumerate(lines):
            protein_name = self.__get_name(line)
            similarity = self.__get_protein_similarity(line)
            if n == 0:
                cluster_df = pd.DataFrame(
                    [
                        [
                            cluster_number,
                            protein_name,
                            similarity,
                            representative
                        ]
                    ],
                    columns = [
                        'Cluster',
                        'Protein',
                        'Similarity',
                        'Representative'
                    ]
                )
            else:
                cluster_df = pd.concat(
                    [
                        cluster_df,
                        pd.DataFrame(
                            [
                                [
                                    cluster_number,
                                    protein_name,
                                    similarity,
                                    representative
                                ]
                            ],
                            columns = [
                                'Cluster',
                                'Protein',
                                'Similarity',
                                'Representative'
                            ]
                        )
                    ],
                    ignore_index = True
                )
        return cluster_df

    def __parse_clstr_file(
        self,
        txt: str
    ) -> pd.DataFrame:

        clusters_txt = txt.split('\n>')[:-1]
        for n, cluster in enumerate(clusters_txt):
            if n == 0:
                clusters = self.__parse_cluster(cluster)
            else:
                clusters = pd.concat(
                    [
                        clusters,
                        self.__parse_cluster(cluster)
                    ],
                    ignore_index = True
                )
        clusters.index = clusters['Protein']
        return clusters

    def to_tsv(
        self,
        path: str
    ) -> None:
        '''
        Writes the clusters DataFrame into a tsv file.
        '''
        self.clusters.to_csv(
            path+'.tsv',
            sep = '\t'
        )

    def __protein_list_from_fasta(
        self,
        path: str
    ) -> list:
        with open(path, 'r') as stream:
            fasta = stream.read()
            if not len(fasta): return []
            proteins = fasta.split('\n>')
        
        protein_names = []
        for protein in proteins:
            if protein.startswith('>'):
                protein = protein[1:]
            protein_names.append(
                protein.split('\n')[0].split()[0]
            )
        return protein_names

    def get_raw_proteins(
        self,
        paths: Iterable[str]
    ) -> pd.DataFrame:
        '''
        Builds a set of protein IDs present in each plasmid. 

        Parameters
        ----------
        
        paths: Iterable
            List of paths to plasmid FASTA files.

        Returns
        -------
            
            Pandas DataFrame of proteins in each plasmid. Index are protein IDs
            and the column 'Plasmid' holds the plasmid ID in which they are
            present
        '''
        accessions = []
        all_proteins = []
        for path in paths:
            accession = path.split('/')[-2]
            proteins = self.__protein_list_from_fasta(path)
            accessions += [accession] * len(proteins)
            all_proteins += proteins

        return pd.DataFrame(
            [[accession] for accession in accessions],
            index = all_proteins,
            columns = ['Plasmid']
        )

    def __find_fasta(
        self,
        folder: str
    ) -> str:
        files = [
            x[2]
            for x in os.walk(folder)
        ][0]
        filename = [x for x in files if x.endswith('.faa')][0]
        return os.path.join(
            folder,
            filename
        )
        
    def replace_with_cluster(
        self,
        raw_proteins: pd.DataFrame
    ) -> pd.DataFrame:

        self.clusters, raw_proteins = self.__filter_for_common_idx(
            self.clusters,
            raw_proteins
        )
        #self.clusters, raw_proteins = self.clusters.drop_duplicates(keep='first'), raw_proteins.drop_duplicates(keep = 'first')
        raw_proteins = raw_proteins.join(
            self.clusters['Representative']
        )
        
        # Group DataFrame according to plasmid id.
        # Join all representative proteins in a list
        self.clustered_proteins = raw_proteins.groupby('Plasmid').apply(
            lambda x: list(x['Representative'].values)
        ).to_frame()
        self.clustered_proteins.columns = ['Proteins']

        return self.clustered_proteins

    def clustered_proteins_to_tsv(
        self,
        path: str
    ):

        self.clustered_proteins.to_csv(
            path+'.tsv',
            sep = '\t'
        )
        joblib.dump(
            deepcopy(self.clustered_proteins),
            path+'.pkl',
            compress = 3
        )


    def build_protein_db(
        self,
        path: str
    ):
        with open(
            self.path_clstr.split('.')[0],
            'r'
        ) as instream:
            reference_protein_file = instream.read()
        
        protein_ids = [
            protein.split()[0]
            for protein in reference_protein_file.split('\n>')
        ]
        protein_ids[0] = protein_ids[0][1:] # Remove first >

        protein_names = [
            ' '.join(protein.split('\n')[0].split()[1:])
            for protein in reference_protein_file.split('\n>')
        ]
        for n, (protein_id, protein_name) in enumerate(
            zip(
                protein_ids,
                protein_names
            )
        ):
            
            protein_names[n] = protein_names[n] + '\n' + protein_id
        
        protein_series = pd.Series(
            protein_names,
            index = protein_ids,
            name = 'Protein Names'
        )
        protein_series.to_csv(
            path,
            sep = '\n'
        )
        return protein_series

    def __get_common_idx(
        self,
        a: pd.DataFrame,
        b: pd.DataFrame
    ) -> list:
        '''
        Returns the index appearing both in DataFrame a and b.
        '''

        return list(
            set(
                a.index
            ).intersection(
                b.index
            )
        )

    def __filter_for_common_idx(
        self,
        a: pd.DataFrame,
        b: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        '''
        Removes all lines with non-common indices in both dataframe
        a and b.

        Parameters
        ----------

        a, b: Pandas DataFrames
            DataFrames to filter.

        Returns
        -------

        a_filtered, b_filtered: Pandas DataFrames
        '''
        idx = self.__get_common_idx(a, b)
        return a.loc[idx], b.loc[idx]


def get_number_cds(
    accessions: Iterable,
    path_to_annotations: str
) -> pd.Series:
    '''
    Given an iterable of plasmid accessions, returns the corresponding
    number of CDS found by Prokka. This information is present in 
    path-to-annotations/path-to-plasmid/PROKKA_xxxxxxxx.txt.
    '''
    paths = find_annotation_paths(accessions, path_to_annotations, format = '.txt')
    data = []
    idx = []
    for path in paths:
        txt_file = pd.read_csv(
                path,
                sep = ': ',
                index_col = 0,
                header = None,
                engine='python'
            )
        try:
            data.append(txt_file.loc['CDS'][1])
        except:
            data.append(0)
            print(f'CDS for accession {path.split("/")[-2]} not found. Set to 0.')
        idx.append(path.split('/')[-2])
    return pd.Series(
        data,
        index = idx
    ).astype(int)
    

def find_annotation_paths(
    accessions: Iterable,
    path_to_annotations: str,
    format: str = '.gff'
):
    walk = list(os.walk(path_to_annotations))
    folders = [
        os.path.join(
            path_to_annotations,
            folder
        )
        for folder in walk[0][1]
    ]
    subfolders = [
        list(
            os.walk(folder)
        )[0][1]
        for folder in folders
    ]
    all_paths = []
    for accession in accessions:
        found = False
        for folder, isubfolders in zip(
            folders,
            subfolders
        ):
            if accession in isubfolders:
                for file in list(os.walk(os.path.join(folder, accession)))[0][2]:
                    if file.endswith(format):
                        all_paths.append(
                            os.path.join(
                                folder,
                                accession,
                                file
                            )
                        )
                        found = True
                        break
        if not found:
            warnings.warn(f'Could not find {format} annotation file for accession {accession}. Skipping this plasmid...')

    return all_paths
# %%
