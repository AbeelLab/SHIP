'''
Implements the GraphGenome and FastBreakDistance classes, used to calculate the
pairwise distance of plasmids.
'''
#%%
from copy import deepcopy
import pandas as pd
import numpy as np
from ship_plasmid.utils.amrfinder import filenames_from_contig_ids
from ship_plasmid.utils.files import find_annotation_paths
from BCBio.GFF import parse
from typing import Iterable, Tuple, Union
import warnings
import networkx as nx
from pyvis.network import Network

class GraphGenome:

    def __init__(
        self,
        id: str,
        path_to_annotations: str,
        path_to_representatives: str,
        path_to_amr: str,
    ):
        self.id = id
        self.path_to_annotations = path_to_annotations
        self.path_to_representatives = path_to_representatives
        self.path_to_amr = path_to_amr
        self.build()

    def __len__(self):
        return len(self.get_gene_set())

    def get_representatives(
        self
    ) -> pd.DataFrame:
        '''
        Builds a Pandas DataFrame containing the representative protein of all proteins in the
        dataset, according to CD-Hit clustering.
        '''

        self.representatives = pd.read_csv(
            self.path_to_representatives,
            sep = '\t',
            index_col = 0
        )
        self.representatives = self.representatives[
            ~self.representatives.index.duplicated(keep='first')
        ]

        return self.representatives

    def get_annotations(
        self
    ) -> pd.DataFrame:
        '''
        Returns a Pandas DataFrame with the gene annotations for the plasmid, as well as strand
        information about each gene.
        '''
        self.path_to_gff = find_annotation_paths(
            accessions = [self.id],
            annotations_dir = self.path_to_annotations,
            format = ['.gff', '.gff3']
        )[0]

        with open(self.path_to_gff, 'r') as file:
            contents = next(parse(file)).features
        
        self.strands = pd.DataFrame(
            [[x.strand, x.qualifiers['product'][0]] for x in contents],
            index = [x.id for x in contents],
            columns = ['Strands', 'Product']
        )

        return self.strands

    def get_amr_gene_names(self) -> list:
        '''
        Returns the gene annotations for all AMR genes in the plasmid.
        '''

        self.amr = pd.read_csv(self.path_to_amr, sep = '\t', index_col = 0)
        contig_ids = filenames_from_contig_ids(self.path_to_annotations)
        self.amr.index = [contig_ids[x] for x in self.amr.index]
        try:
            self.amr = self.amr.loc[self.id]
            if type(self.amr['Sequence name'])==pd.Series:
                self.__amr_names = self.amr['Sequence name'].values
            else:
                self.__amr_names = [self.amr['Sequence name']]
        except KeyError:
            self.__amr_names = []

        return self.__amr_names

    def join_gene_info(
        self,
        warn: bool = True
    ) -> pd.DataFrame:
        '''
        Joins the information about the representative protein in the cluster, strand, and
        AMR status for each gene in the plasmid into a single Pandas DataFrame.
        '''

        common_idx = list(set(self.representatives.index).intersection(
            self.strands.index
        ))
        not_found = list(set(self.strands.index).difference(
            self.representatives.index
        ))

        if warn:
            if len(not_found) > 0:
                warnings.warn(f'Could not find CD-Hit information for proteins {not_found}. These proteins will be removed from the genome.')
            if len(common_idx) == 0:
                raise Exception('Could not find any codified protein with CD-Hit info.')

        # [x for x in self.strands.index if x in common_idx] seems like bad coding
        # but it preserves the original sequence encoded in self.strands.index. 
        # The list common_index is not sorted according to genome synteny. Thus, 
        # using it to sort the DataFrame would yield incorrect results
        self.representatives = pd.concat(
            [
                self.representatives.loc[
                    common_idx
                ],
                self.strands
            ],
            join = 'inner',
            axis = 'columns'
        ).drop(['Cluster'], axis = 'columns').loc[
            [x for x in self.strands.index if x in common_idx]
        ]
        
        # AMRFinder annotations may differ from those in Prokka, for the same genes. 
        # Compare the positions of the AMR hits with the ORFs to detect resistant genes. 
        # Allows for a 200 bp mismatch in either side
        if type(self.amr) == pd.Series:
            self.amr = self.amr.to_frame().transpose()
        is_amr = self.representatives['Representative'].apply(
            self.amr_func
        )
    
        is_amr.name = 'Is AMR'

        self.representatives = pd.concat(
            [self.representatives, is_amr],
            axis = 'columns'
        )

        return self.representatives

    def amr_func(self, x):
        out = []
        for n in range(len(self.amr)):
            out.append(
                np.logical_and(
                    np.abs((self.positions.loc[x]['Start'] - self.amr.iloc[n]['Start'])) < 200,
                    np.abs((self.positions.loc[x]['End'] - self.amr.iloc[n]['Stop'])) < 200
                )
            )

        return np.any(out)

    def remove_duplicates(
        self
    ) -> dict:
        '''
        Deletes repeated genes from the representatives DataFrame (containing annotations,
        strand information and AMR status for all proteins in the plasmid). Calculates the
        number of total duplications, stored in the attribute duplicates, as well as the
        number of extra copies for each gene in duplicates_counts.

        Returns
        -------

        duplicates_counts: dict
        '''
        self.representatives_with_duplicates = deepcopy(self.representatives)
        self.duplicates = 0
        self.duplicates_counts = {}
        unique_prots, counts = np.unique(
            self.representatives['Representative'],
            return_counts = True
        )
        duplicate_prots = unique_prots[counts > 1]
        for prot in duplicate_prots:
            duplicate_df = self.representatives[self.representatives['Representative']==prot]
            duplicate_idx = duplicate_df.index[1:]
            self.representatives.drop(duplicate_idx, inplace=True)
            self.duplicates += len(duplicate_idx)
            self.duplicates_counts[prot] = len(duplicate_idx)
    
        return self.duplicates_counts

    def get_gene_sizes(
        self,
        include_duplicates = False
    ) -> pd.DataFrame:
        '''
        Returns the positions and gene sizes in base pairs for each gene in the plasmid.
        '''
        with open(self.path_to_gff, 'r') as file:
            contents = next(parse(file)).features
        
        if include_duplicates:
            index_df = self.representatives_with_duplicates
        else:
            index_df = self.representatives

        try:
            positions = pd.DataFrame(
                [
                    [
                        x.location.start.position, 
                        x.location.end.position,
                        x.location.end.position - x.location.start.position
                    ] for x in contents if x.id in index_df.index
                ],
                index = [
                    index_df.loc[x.id]['Representative'] 
                    for x in contents 
                    if x.id in index_df.index
                ],
                columns = ['Start', 'End', 'Size']
            )
        except:
            positions = pd.DataFrame(
                [
                    [
                        x.location.start.position, 
                        x.location.end.position,
                        x.location.end.position - x.location.start.position
                    ] for x in contents if x.id in index_df.index
                ],
                index = [
                    index_df.loc[x.id]['Representative'].iloc[0]
                    for x in contents 
                    if x.id in index_df.index
                ],
                columns = ['Start', 'End', 'Size']
            )

        if include_duplicates:
            self.positions_with_duplicates = positions
        else:
            self.positions = positions

        return positions

    def get_gene_set(
        self
    ):
        '''
        Returns the ids of all representative proteins in the plasmid as a list,
        with duplicates removed.
        '''
        return [
            self.representatives.iloc[n]['Representative']
            for n in range(len(self.representatives))
        ]

    def build_graph(
        self
    ):

        self.genome = nx.Graph(name = self.id)
        if len(self.representatives) > 1:
            for gene_n in range(0, len(self.representatives)-1):
                self.__add_gene(
                    self.representatives.iloc[gene_n]['Representative'],
                    self.representatives.index.to_list()[gene_n]+': '+self.representatives.iloc[gene_n]['Product'],
                    self.representatives.iloc[gene_n]['Is AMR']
                )
                self.__add_edge(
                    self.representatives.iloc[gene_n]['Representative'],
                    self.representatives.iloc[gene_n]['Strands'],
                    self.representatives.iloc[gene_n+1]['Representative'],
                    self.representatives.iloc[gene_n+1]['Strands']
                )
            self.__add_gene(
                self.representatives.iloc[-1]['Representative'],
                self.representatives.index.to_list()[-1]+': '+self.representatives.iloc[-1]['Product'],
                self.representatives.iloc[-1]['Is AMR']
            )
            self.__add_edge(
                self.representatives.iloc[-1]['Representative'],
                self.representatives.iloc[-1]['Strands'],
                self.representatives.iloc[0]['Representative'],
                self.representatives.iloc[0]['Strands']
            )

        elif len(self.representatives) == 1:
            self.__add_gene(
                self.representatives.iloc[0]['Representative'],
                self.representatives.index.to_list()[0]+': '+self.representatives.iloc[0]['Product'],
                self.representatives.iloc[0]['Is AMR']
            )
            self.__add_edge(
                self.representatives.iloc[-1]['Representative'],
                self.representatives.iloc[-1]['Strands'],
                self.representatives.iloc[0]['Representative'],
                self.representatives.iloc[0]['Strands']
            )
        else:
            warnings.warn(f'This plasmid ({self.id}) has no CDS.')

        self.synteny_edges = [
            (u, v) for u, v, att in self.genome.edges(data=True) 
            if att['color'] == 'synteny'
        ]
        self.gene_edges = [
            (u, v) for u, v, att in self.genome.edges(data=True) 
            if att['color'] == 'gene'
        ]
        
    def build(
        self
    ):
        '''
        Constructs the graph of the plasmid genome.
        '''
        self.get_representatives()
        self.get_annotations()
        self.get_gene_sizes()
        self.get_amr_gene_names()
        self.join_gene_info()
        self.remove_duplicates()

        self.build_graph()

    def __add_gene(
        self,
        name,
        product,
        is_amr
    ):
        '''
        Adds a gene to the genome graph, with the required edges and nodes.
        '''
        for side in ['Head', 'Tail']:
            self.genome.add_node(
                ' '.join([name, side]),
                title = product,
                is_amr = is_amr
            )
        self.genome.add_edge(
            name + ' Head',
            name + ' Tail',
            label = name,
            color = 'gene',
            title = product,
            is_amr = is_amr
        )

    def __add_edge(
        self,
        start_gene,
        start_strand,
        end_gene,
        end_strand
    ):
        '''
        Adds an edge to the genome graph.
        '''

        if start_strand > 0:
            start_gene = start_gene + ' Tail'
        else:
            start_gene = start_gene + ' Head'

        if end_strand > 0:
            end_gene = end_gene + ' Head'
        else:
            end_gene = end_gene + ' Tail'

        self.genome.add_edge(
            start_gene,
            end_gene,
            color = 'synteny',
            label = self.id
        )

class FastBreakDistance:

    def __init__(
        self,
        k: int = 3,
        indel_gene_penalty: Union[float, None] = 0.5
    ):
        self.k = k
        self.__called = False
        assert self.k in [2, 3], f'NOT IMPLEMENTED: k must be either 2 or 3. Got {k}.'
        self.MAX_DISTANCE = 1000
        self.distances = {
            'mutation': 0,
            'mobile element': 0,
            'recombination': 0,
            'reversal': 0,
            'transposition': 0,
            'duplicate': 0
        }

    def extract_duplicates(
        self,
        u: GraphGenome,
        v: GraphGenome
    ):
        '''
        Gets the duplicates_counts attributes from the GraphGenome objects.
        '''
        self.duplicates = {
            k.id: deepcopy(k.duplicates_counts)
            for k in [u, v]
        }

    def score_indels(
        self,
        u: GraphGenome,
        v: GraphGenome
    ):

        # Get genes present u bot not in v
        exclusive_genes = list(set(u.get_gene_set()).difference(v.get_gene_set()))
        
        # Build a subgraph of u containing only the different genes
        subgraph = nx.subgraph(
            u.genome,
            [
                ' '.join([x, 'Head']) for x in exclusive_genes
            ] + [
                ' '.join([x, 'Tail']) for x in exclusive_genes
            ]
        )

        # Iterate over the graph connected components. Each connected component corresponds
        # to a syntenic block
        components = nx.connected_components(subgraph)
        # End position of the last parsed component
        for n, component in enumerate(components):
            # Find start and end of the components by searching for single degree nodes
            degree = [x for x in nx.degree(subgraph, component)]
            limits = [x[0].split()[0] for x in degree if x[1] == 1]

            # If the component is a cycle, that means that there are no genes in common.
            no_overlap = False
            if len(limits) == 0:
                limits = [u.positions.iloc[0].name, u.positions.iloc[-1].name]
                no_overlap = True
            # Retrieve position in bp for each limit
            positions = self.__get_positions(u, limits, component, no_overlap)
            # Get component length in bp
            length = positions[0].iloc[-1]['End'] - positions[0].iloc[0]['Start']
            if len(positions) == 2:
                length += positions[1].iloc[-1]['End'] - positions[1].iloc[0]['Start']
            # Call scoring callback
            self.__score_callback(length, component, u, v)

            for position in positions:
                # Exclude duplicates inside the block from the duplicate count
                self.__duplicate_check_callback(u, position.iloc[0]['Start'], position.iloc[-1]['End'])

    def __get_positions(
        self,
        genome: GraphGenome,
        limits: Tuple,
        component: set,
        no_overlap: bool
    ) -> Union[Tuple[pd.DataFrame], Tuple[pd.DataFrame, pd.DataFrame]]:
        '''
        Computes the positions of the connected component. Identifies if it contains the end and
        start of the plasmid sequence. If it does, divides the positions DataFrame into two,
        containing the interval [sequence start, component start] and [component end, sequence end].
        '''
        # Sort the start and end genes by their position, in ascending order
        positions = genome.get_gene_sizes().loc[limits].sort_values('Start')
        if len(component) < 4:
            return [positions]
        elif len(component) == 4:
            if (positions.iloc[-1]['End'] - positions.iloc[0]['Start']) == (genome.positions.iloc[-1]['End'] - genome.positions.iloc[0]['Start']) and not no_overlap:
                # The component spans the start
                return (
                    [
                        pd.concat(
                            [
                                positions.iloc[-1:],
                                genome.positions.iloc[-1:]
                            ]
                        ),
                        pd.concat(
                            [
                                genome.positions.iloc[0:1],
                                positions.iloc[0:1]
                            ]
                        ),
                    ]
                )
            else:
                return [positions]
        # Take another gene from the component. If it falls between the start and
        # end genes, the component does not span the begining/end of the gene
        # sequence
        probe_gene = list(set(x.split()[0] for x in component).difference(limits).difference(list(genome.duplicates_counts.values())))[0]
        start_pos, end_pos = positions.iloc[0]['Start'], positions.iloc[-1]['End']
        probe_pos = genome.positions.loc[probe_gene]['Start']
        if not (probe_pos > start_pos and probe_pos < end_pos):
            # The component includes the start/end of the sequence. Split the
            # positions DataFrame into two: [start, component start] U [component
            # end, end]
            return (
                [
                    pd.concat(
                        [
                            positions.iloc[-1:],
                            genome.positions.iloc[-1:]
                        ]
                    ),
                    pd.concat(
                        [
                            genome.positions.iloc[0:1],
                            positions.iloc[0:1]
                        ]
                    ),
                ]
            )
        else:
            return [positions]

    def __score_callback(
        self,
        length: int,
        component: set,
        u: GraphGenome,
        v: GraphGenome
    ):
        '''
        Updates the recombination, mutation and mobile element distance according to the length
        of an indel block.

        The thresholds were computed using MaxLikelihood, gene sizes and ISFinder databases. The
        original code is available in /Scripts/Unused/Quality Check and 
        Testing/insertion-sequence-lengths.py
        '''

        # Length thresholds
        self.T_MUT = 915
        self.T_HR = 2752

        if length < self.T_MUT:
            if self.__is_mutation(
                component,
                u,
                v
            ):
                self.distances['mutation'] += 1
            else:
                self.distances['mobile element'] += 1
        elif length < self.T_HR:
            self.distances['mobile element'] += 1
        else:
            self.distances['recombination'] += 1

    def __duplicate_check_callback(
        self,
        genome: GraphGenome,
        start: int,
        end: int
    ):
        '''
        Checks if there are genes with multiple copies in the synteny block. If there are, deducts the 
        number of copies in the block to the number of duplicates in the plasmid. If a gene is less than
        1000 bp from the synteny block, it is considered inside it.

        This reduces the chance of including IS, especially IS6, proteins likely gained through HR or 
        mobile elements in the number of duplications.
        '''
        margin = 000
        positions = genome.get_gene_sizes(include_duplicates=True)
        for gene in genome.duplicates_counts:
            gene_positions = positions.loc[gene]
            included = gene_positions.loc[
                np.logical_and(
                    (start - gene_positions['Start']) < margin,
                    (gene_positions['End'] - end) < margin
                )
            ]
            self.duplicates[genome.id][gene] -= len(included)
            if self.duplicates[genome.id][gene] == -1:
                self.duplicates[genome.id][gene] = 0
            elif self.duplicates[genome.id][gene] < -1:
                self.duplicates[genome.id][gene] = 0
                warnings.warn(f'Found duplicate distance for gene {gene} in plasmid {genome.id} when comparing {self.u.id} and {self.v.id} less than 0: {self.duplicates[genome.id][gene]}. Positions: {[start, end]}. Setting the number of duplications to 0.')

    def __is_mutation(
        self,
        gene: set,
        u: GraphGenome,
        v: GraphGenome
    ) -> bool:
        '''
        If there is more than one gene in the block, return False.
        Checks if there is any gene in genome v which has the same neighbors as in the gene
        in genome u. If the query gene is X and it is present in u as B-X-C, the method
        returns true if and only if there is a sequence B-Y-C, where Y != X, in v.
        '''
        
        if len(gene) > 1: return False
        gene = list(gene)[0]

        # Number of genes shared by u and v, in each side of the query gene, needed 
        # to consider the gene a result of extensive mutation.
        flank_size = 1 # Do not change. It will not work otherwise
        # Find the position of the gene in genome u
        gene_idx = np.argmax(u.representatives['Representative'] == gene)
        # Get the flanking genes
        flanks = [
            u.representatives.iloc[gene_idx-1]['Representative'],
            u.representatives.iloc[gene_idx+1]['Representative']
        ]
        # Find the flanking gene in v. As the duplicates have been removed, there is only one
        # instance of the flanking gene in v
        v_pos = np.argmax(v.representatives['Representative'] == flanks[0])
        # Check the opposite flank
        if v.representatives.iloc[v_pos+1][
            'Representative'
        ] == flanks[1] or v.representatives.iloc[v_pos-1][
            'Representative'
        ] == flanks[1]:
            return True
        else:
            return False

    def drop_non_common(
        self,
        u: GraphGenome,
        v: GraphGenome
    ) -> Tuple[GraphGenome, GraphGenome]:
        '''
        Creates new v and u graphs, with only the genes in common between the two
        '''

        # Discard genes not in common
        common_genes = set(u.get_gene_set()).intersection(v.get_gene_set())
        new_graphs = []
        for n, graph in enumerate([u, v]):
            new_graph = deepcopy(graph)
            # Preserves the original order
            new_graph.representatives = graph.representatives.drop(
                graph.representatives.index.to_numpy()[
                    [
                        x not in common_genes 
                        for x in graph.representatives['Representative']
                    ]
                ]
            )
            new_graph.build_graph()
            new_graphs.append(new_graph)
        return new_graphs

    def make_synteny_loops(
        self,
        u: GraphGenome,
        v: GraphGenome
    ) -> nx.Graph:
        '''
        Joins two genome graphs, deleting the edges connecting genes. This allows for the
        identification of syntenic edge loops, needed for calculating break distances.

        Parameters
        ----------

        u, v: NetworkX Graphs
            Genome graphs with insertions and deletions removed. U and V should have the same
            nodes, but can have different edges. They can be obtained by calling 
            self.remove_indels().

        Returns
        -------

        joined_graph: NetworkX Graph
            Graph with all synteny edges in u and v.
        '''

        joined_graph = nx.Graph(name = ' '.join([u.id, v.id]))
        joined_graph.add_nodes_from(u.genome.nodes(data=True))
        for graph in [u.genome, v.genome]:
            synteny_edges = [
                (a, b) 
                for a, b, data in graph.edges(data=True) 
                if data['color']=='synteny'
            ]

            joined_graph.add_edges_from(
                synteny_edges,
                color = graph.name+' Synteny'
            )
        
        return joined_graph

    def __get_common_components(
        self,
        graph: Union[nx.MultiGraph, nx.Graph],
        return_sizes: bool = True
    ) -> Union[Tuple[list, list], list]: 

        components = list(nx.connected_components(graph))
        sizes = [
            len(
                graph.subgraph(x).edges()
            ) 
            for x in components
        ] 
        if return_sizes: return components, sizes
        else: return components

    def __is_odd(
        self,
        x: Iterable
    ):
        return [xx % 2 != 0 for xx in x]

    def three_break_distance(
        self,
        common_graph: nx.MultiGraph
    ) -> int:
        common_components, sizes = self.__get_common_components(
            common_graph,
            return_sizes = True
        )
        is_odd = self.__is_odd(sizes)
        n_odd_components = sum(is_odd)
        n_nodes = len(common_graph.nodes())

        return np.ceil((n_nodes/2 - n_odd_components)/2)

    def two_break_distance(
        self,
        common_graph: nx.MultiGraph
    ) -> int:
        common_components = self.__get_common_components(
            common_graph,
            return_sizes = False
        )
        n_nodes = len(common_graph.nodes())

        return n_nodes/2 - len(common_components)

    def duplicate_distance(
        self,
        u: GraphGenome,
        v: GraphGenome
    ):

        distance = 0

        common_duplicates = set(
            self.duplicates[u.id].keys()
        ).intersection(
            self.duplicates[v.id].keys()
        )

        for gene in common_duplicates:
            distance += np.abs(self.duplicates[u.id][gene] - self.duplicates[v.id][gene])

        return distance

    def __call__(
        self,
        u: GraphGenome,
        v: GraphGenome
    ) -> float:
        self.u = u
        self.v = v
        self.__called = True
        self.distance_ = self.calculate(
            u,
            v
        )
        return self.distance_

    def calculate(
        self,
        u: GraphGenome,
        v: GraphGenome
    ) -> float:

        # Reset distances
        self.distances = {
            'mutation': 0,
            'mobile element': 0,
            'recombination': 0,
            'reversal': 0,
            'transposition': 0,
            'duplicate': 0
        }
        self.extract_duplicates(u, v)
        self.score_indels(u, v)
        self.score_indels(v, u)
        common_u, common_v = self.drop_non_common(u, v)
        self.distances['duplicate'] = self.duplicate_distance(u, v)
        self.common_graph = self.make_synteny_loops(common_u, common_v)
        if len(common_u.genome.nodes) <= 2: 
            self.distances['reversal'] = 0
            self.distances['transposition'] = 0
        else:
            self.distances['reversal'] = self.two_break_distance(self.common_graph)
            if self.k == 3:
                self.distances['transposition'] = self.three_break_distance(self.common_graph)
            elif self.k == 2:
                self.distances['transposition'] = 0
            else:
                raise NotImplementedError(f'Only implemented k are 2 and 3. Got {self.k}.')

        self.distances['mutation'] /= 2
        return sum([v for v in self.distances.values()])
