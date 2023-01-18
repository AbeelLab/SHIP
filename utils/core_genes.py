'''
Classes for panplasmidome analysis and phylogeny inference. Require SynerClust and/or Roary.
'''

from typing import Iterable, Union
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, to_tree, dendrogram
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import subprocess
import logging
import shutil


class PhylogeneticTree:

    def __init__(
        self
    ):
        pass

    def from_distance(
        self,
        distances: pd.DataFrame
    ):
        self.leaf_names = distances.index.to_list()
        condensed_distance = squareform(1-distances.values, checks=False)
        self.links = linkage(condensed_distance)
        self.scipy_tree = to_tree(self.links, False)

        self.newick_tree = self.__get_newick(
            self.scipy_tree, 
            self.scipy_tree.dist, 
            self.leaf_names
        )

        return self.newick_tree

    def __get_newick(
        self,
        node, 
        parent_dist, 
        leaf_names, 
        newick=''
    ) -> str:
        """
        Convert sciply.cluster.hierarchy.to_tree()-output to Newick format. Taken from
        https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format

        Parameters
        ----------
        
        node: 
            output of sciply.cluster.hierarchy.to_tree()
        parent_dist: 
            output of sciply.cluster.hierarchy.to_tree().dist
        leaf_names: 
            list of leaf names
        newick: 
            leave empty, this variable is used in recursion.
        
        Returns
        -------
            Tree in Newick format
        """
        if node.is_leaf():
            return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = "):%.2f%s" % (parent_dist - node.dist, newick)
            else:
                newick = ");"
            newick = self.__get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
            newick = self.__get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
            newick = "(%s" % (newick)
            return newick

    def plot_dendogram(
        self
    ):
        with plt.style.context('ggplot'):
            fig = plt.figure(
                figsize = (10,7)
            )
            dendrogram(
                self.links,
                labels = self.leaf_names,
                leaf_rotation = 90,
                distance_sort = 'descending'
            )
            plt.show()
        return fig

class SynerClustInput:

    def __init__(
        self,
        accessions: Iterable,
        path_to_annotations: str,
        path_to_sequences: str,
        path_to_output: str,
    ):
        '''
        Prepares the SynerClust input files by building the catalog text file with the paths to
        plasmid sequences and GFF annotations.

        Arguments
        ---------

        accessions: Iterable
            Iterable of accession numbers (without the trailing .1) of plasmids to include in the
            analysis.

        path_to_annotations: str
            Directory containing all annotation subfolders.

        path_to_sequences: str
            Directory containing all FASTA sequences.

        path_to_output: str
            Directory in which to write the input files for SynerClust.
        '''
        self.paths = {
            'annotations': path_to_annotations,
            'output': path_to_output,
            'sequences': path_to_sequences
        }
        self.accessions = accessions

        if not os.path.exists(
            self.paths['output']
        ):
            os.mkdir(self.paths['output'])

    def find_sequence_paths(
        self,
        accessions: Iterable
    ) -> Iterable:

        new_paths = []
        walk = [x for x in os.walk(self.paths['sequences'])]
        folders = walk[0][1]

        for plasmid in accessions:
            query = plasmid+'.fa'
            found_file = False
            for folder in folders:
                test_path = os.path.join(
                    self.paths['sequences'],
                    folder,
                    query
                )
                if os.path.exists(test_path):
                    new_paths.append(test_path)
                    found_file = True
                    break

            if not found_file:
                warnings.warn(f'Could not find FASTA file for accession {plasmid}. Skipping this plasmid...')

        self.__sequences = new_paths
        return self.__sequences

    def find_annotation_paths(
        self,
        accessions: Iterable,
        format: str = '.gff'
    ):
        walk = list(os.walk(self.paths['annotations']))
        folders = [
            os.path.join(
                self.paths['annotations'],
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

        self.__annotations = all_paths
        return self.__annotations

    def write_catalog(
        self,
        accessions: Union[Iterable, None]
    ) -> None:
        '''
        Writes a text file containing the genome name (set to its accession), path to the FASTA
        sequences, and path to the GFF annotations for each plasmid. Acts as input for SynerClust.
        The text file is written in /output_path/Catalog/plasmid_catalog.txt

        Arguments
        ---------

        accessions: Iterable
            Iterable of accession numbers (without the trailing .1) of plasmids to include in the
            analysis.
        '''

        if not os.path.exists(
            os.path.join(
                self.paths['output'],
                'Catalog'
            )
        ):
            os.mkdir(
                os.path.join(
                    self.paths['output'],
                    'Catalog'
                )
            )
        self.paths['catalog'] = os.path.join(
            self.paths['output'],
            'Catalog',
            'plasmid_catalog.txt'
        )
        if not os.path.exists(self.paths['catalog']):
            open_mode = 'x'
        else:
            open_mode = 'w'

        if accessions is None: accessions = self.accessions
        # Holds all paths to plasmid annotations and sequences
        self.__annotations = self.find_annotation_paths(accessions)
        self.__sequences = self.find_sequence_paths(accessions)
        with open(
            self.paths['catalog'],
            open_mode
        ) as out_stream:
            
            for id, sequence, annotation in zip(
                accessions,
                self.__sequences,
                self.__annotations
            ):
                out_stream.write(
f'''\\\\
Genome  {id}
Sequence    {sequence}
Annotation  {annotation}
'''
                )

    def write_tree(
        self,
        newick: str
    ):
        '''
        Writes a text file containing a Newick tree, to act as input for SynerClust.
        The text file is written in /output_path/Catalog/tree.txt.

        Arguments
        ---------

        newick: str
            Initial guiding tree in newick format. Can be obtained by calling 
            PhylogeneticTree().from_distance().
        '''
        if not os.path.exists(
            os.path.join(
                self.paths['output'],
                'Catalog'
            )
        ):
            os.mkdir(
                os.path.join(
                    self.paths['output'],
                    'Catalog'
                )
            )

        self.paths['tree'] = os.path.join(
            self.paths['output'],
            'Catalog',
            'tree.txt'
        )
        if not os.path.exists(self.paths['tree']):
            open_mode = 'x'
        else:
            open_mode = 'w'

        with open(
            self.paths['tree'],
            open_mode
        ) as out_stream:
            out_stream.write(newick)


class SynerClust:
    '''
    Runs SynerClust with the selected plasmids. Uses the PhylogeneticTree and SynerClustInput
    classes internally.
    '''

    def __init__(
        self,
        accessions: Iterable,
        initial_affinity: pd.DataFrame,
        path_to_annotations: str,
        path_to_sequences: str,
        path_to_output: str 
    ):
        '''
        Runs SynerClust with the selected plasmids. Uses the PhylogeneticTree and SynerClustInput
        classes internally.

        Arguments
        ---------

        accessions: Iterable
            Iterable of accession numbers (without the trailing .1) of plasmids to include in the
            analysis.

        initial_affinity: Pandas DataFrame
            Affinity matrix with which to build the initial guiding tree for SynerClust.

        path_to_annotations: str
            Directory containing all annotation subfolders.

        path_to_sequences: str
            Directory containing all FASTA sequences.

        path_to_output: str
            Directory in which to write the input files for SynerClust.
        '''
        self.__path_to_output = path_to_output
        if not os.path.exists(
            self.__path_to_output
        ):
            os.mkdir(self.__path_to_output)

        self.__set_logging_options()
        self.SYNERCLUST_PATH = os.path.expanduser('~')+'/SynerClust-master/bin/synerclust.py'
        logging.debug(f'Set SynerClust path as {self.SYNERCLUST_PATH}.')
        
        logging.debug(f'Preparing input files for accessions {accessions}.')
        self.accessions = accessions
        logging.debug(f'Building Newick tree.')
        try:
            self.newick = PhylogeneticTree().from_distance(initial_affinity)
        except Exception as e:
            logging.error(f'Tree building failed with exception {e}.')
            raise e
        logging.debug(f'Instanciating SynerClustInput object.')
        try:
            self.input = SynerClustInput(
                accessions,
                path_to_annotations,
                path_to_sequences,
                path_to_output
            )
        except Exception as e:
            logging.error(f'Instanciating failed with exception {e}.')
            raise e

        logging.debug(f'Writing catalog file for SynerClust in {self.__path_to_output}/Catalog.')
        try:
            self.input.write_catalog(self.accessions)
        except Exception as e:
            logging.error(f'Catalog writing failed with exception {e}.')
            raise e
        logging.debug(f'Writing Newick tree txt file for SynerClust.')
        try:
            self.input.write_tree(self.newick)
        except Exception as e:
            logging.error(f'Tree writing failed with exception {e}.')
            raise e
        
        self.paths = self.input.paths

    def run(
        self,
        cwd: Union[str, None] = None,
        env: str = 'synerclust',
        base_env: str = 'DBL'
    ):
        if cwd is None: cwd = os.getcwd()

        command = [
            'conda',
            'activate',
            'synerclust',
            '&&',
            self.SYNERCLUST_PATH,
            '-r',
            self.paths['catalog'],
            '-t',
            self.paths['tree'],
            '-w',
            cwd+'/',
            '-n',
            str(12),
            '--run single',
            '&&',
            'conda',
            'activate',
            base_env
        ]

        logging.debug(f'Calling SynerClust using command {" ".join(command)}.')
        process = subprocess.run(command, shell = True)
        logging.debug(f'SynerClust return code was {process.returncode}.')
        if process.returncode != 0:
            logging.error(f'SynerClust returned with code {process.returncode}.')
            raise Exception(f'SynerClust returned with code {process.returncode}.')

    def __set_logging_options(
        self,
        filename: str = 'log.txt'
    ):
        self.log_path = os.path.join(
            self.__path_to_output,
            filename
        )
        logging.basicConfig(
            filename = self.log_path, 
            filemode='w',
            format ='%(process)d :: %(levelname)s :: %(message)s',
            level = logging.DEBUG
            )

class Roary:

    def __init__(
        self,
        accessions: Iterable,
        path_to_annotations: str,
        path_to_output: str,
    ):
        '''
        Prepares the Roary input files by copying and pasting all GFF annotations
        for a set of accession IDs into the same folder.

        Arguments
        ---------

        accessions: Iterable
            Iterable of accession numbers (without the trailing .1) of plasmids to include in the
            analysis.

        path_to_annotations: str
            Directory containing all annotation subfolders.

        path_to_output: str
            Directory in which to write the input files for Roary and its output.
        '''
        self.__path_to_output = path_to_output
        logging.debug('Creating output paths if missing.')
        if not os.path.exists(
            self.__path_to_output
        ):
            os.mkdir(self.__path_to_output)
            os.mkdir(
                os.path.join(
                    self.__path_to_output,
                    'GFF_Annotations'
                )
            )
            logging.debug(f'Created missing output paths in {self.__path_to_output}.')
        else:
            logging.debug(f'No output missing directories.')
            shutil.rmtree(
                os.path.join(
                    self.__path_to_output,
                    'GFF_Annotations'
                )
            )
            logging.debug(f'Removed existing annotations in {os.path.join(self.__path_to_output, "GFF_Annotations")}.')
        self.__set_logging_options()
        logging.debug('Setting paths.')
        self.paths = {
            'annotations': path_to_annotations,
            'output': path_to_output,
            'roary': os.path.join(
                self.__path_to_output,
                'Roary_Output'
            ),
            'gff': os.path.join(
                self.__path_to_output,
                'GFF_Annotations'
            )
        }
        logging.debug(f'Preparing Roary inputs with paths {self.paths}')
        self.accessions = accessions
        logging.debug(f'Will process accessions {accessions}.')
            

    def __find_annotation_paths(
        self,
        accessions: Iterable,
        format: str = '.gff'
    ):
        walk = list(os.walk(self.paths['annotations']))
        folders = [
            os.path.join(
                self.paths['annotations'],
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

        self.__annotations = all_paths
        return self.__annotations 

    def __prepare_input(
        self
    ):
        logging.debug('Looking for annotation GFF files.')

        self.__annotations = self.__find_annotation_paths(self.accessions)
        logging.debug(f'Copying annotations to {self.paths["gff"]}')
        for path in self.__annotations:
            accession = path.split('/')[-2]
            shutil.copyfile(
                path,
                os.path.join(
                    self.paths['gff'],
                    accession+'.gff'
                )
            )
        logging.debug('Finished copying GFF files.')

    def run(
        self
    ):
        self.__prepare_input()
        os.chdir(self.paths["gff"])
        command = [
            'cd',
            self.paths['gff'],
            '&&'
            'roary',
            '-e',
            '-f',
            self.paths['roary'],
            '-i',
            '90',
            '--mafft',
            '*.gff'
        ]

        logging.debug(f'Running Roary with command {command}')
        code = subprocess.run(command, shell=True).returncode
        if code != 0:
            logging.error(f'Roary returned with code {code}.')
            raise Exception(f'Roary returned with code {code}.')
        else:
            logging.debug('Roary returned with code 0.')
        logging.debug('Done!')

    def __set_logging_options(
        self,
        filename: str = 'log.txt'
    ):
        self.log_path = os.path.join(
            self.__path_to_output,
            filename
        )
        logging.basicConfig(
            filename = self.log_path, 
            filemode='x',
            format ='%(process)d :: %(levelname)s :: %(message)s',
            level = logging.DEBUG
            )