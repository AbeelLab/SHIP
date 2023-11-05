'''
Classes to compute pairwise distances of plasmids and display detailed plasmid similarity
networks.
'''

from copy import deepcopy
from typing import Callable, Iterable, Tuple, Union
import joblib
import pandas as pd
import os
import numpy as np
from ship_plasmid.utils.jaccard import jaccard_distance
from ship_plasmid.utils.genome import GraphGenome
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pyvis.network import Network
import webbrowser
from time import sleep
import pymc as pm
import arviz as az
import networkx as nx

class PlasmidDistance:

    def __init__(
        self,
        accessions: Iterable,
        path_to_annotations: str,
        path_to_representatives: str,
        path_to_amr: str,
        clustering_method: Callable,
        learn_weights: bool = True,
        learning_method: str = 'ratio',
        weights: Union[list, None] = None
    ):
        '''
        Parameters
        ----------

        accessions: Iterable
            Accession numbers of the plasmids to include in the dendogram.

        path_to_annotations: str
            Path to the directory containing annotations in GFF format. Usually
            Data/Annotations.

        path_to_representatives: str
            Path to the TSV file containing the protein cluster membership comming 
            from CD-HIT, and their representative protein. Usually, this file is 
            in InitialClustering/Protein Cluster Membership/protein_cluster_membership.tsv.

        clustering_method: Callable
            Clustering object implementing a fit_predict function taking an affinity matrix.
            Can be scikit learn's implementation of Agglomerative Clustering or other.

        learn_weights: bool. Default is True
            If True, each distance module d_k will be weighted according to the weights learned 
            by the specified learning method.

        learning_method: str. Default is "ratio"
            Method according to which to update the weights given to each distance component.
            Must be either "ratio", "predefined", or "max-likelihood". Ignored if learn_weights
            is False.

                "ratio": Sets w_i = mean(N_i/N_total), i.e., the mean proportion of changes of
                type i.

                "predefined": User defined weights, as a dict, passed to the parameter "weights".

                "maximum-likelihood": Parameters estimated through MLE, sampling from the posterior,
                assuming Poisson distributions for all distance components.

        weights: dict or None. Default is None
            Ignored if learning_method is not "predefined".

        '''

        self.accessions = accessions
        self.path_to_annotations = path_to_annotations
        self.path_to_representatives = path_to_representatives
        self.clustering_method = deepcopy(clustering_method)
        if weights is None:
            self.weights = {
                'mutation': 0,
                'mobile element': 0,
                'recombination': 0,
                'reversal': 0,
                'transposition': 0,
                'duplicate': 0
            }
        else:
            self.weights = weights
        self.__learn = learn_weights
        self.path_to_amr = path_to_amr
        self.learning_method = learning_method
        self.HR_DIST = 1.1

    def __build_affinity_matrix(
        self,
        dist_f: Callable
    ):
        '''
        Parameters
        ----------

        dist_f: Callable.
            Distance function taking two graph genomes (GraphGenome objects) as input
            and returning a number (distance between the two genomes).

        Returns
        -------

        affinity_matrix: Numpy array
        '''

        # Initialize affinity matrix
        affinity = np.zeros(
            (len(self.accessions), len(self.accessions))
        )
        self.jaccard = np.zeros_like(affinity)
        self.distances = {
            name_ : np.zeros_like(affinity)
            for name_ in self.weights
        }
        self.sizes = np.zeros_like(affinity)

        graph_genomes = {}
        for accession in self.accessions:
            graph_genomes[accession] = GraphGenome(
                accession,
                path_to_representatives = self.path_to_representatives,
                path_to_annotations = self.path_to_annotations,
                path_to_amr = self.path_to_amr
            )
        
        total_n_genomes = len(self.accessions)
        for n, u_accession in enumerate(self.accessions):
            for m, v_accession in enumerate(self.accessions):
                if m > n:
                    print(f'Computing pairwise distances. {int(((n*total_n_genomes+m)/(total_n_genomes**2))*100)}% Done')
                    affinity[n, m] = dist_f(
                        graph_genomes[u_accession],
                        graph_genomes[v_accession]
                    )
                    affinity[m, n] = affinity[n, m]
                    self.jaccard[n, m] = jaccard_distance(
                        graph_genomes[u_accession].get_gene_set(),
                        graph_genomes[v_accession].get_gene_set()
                    )
                    self.jaccard[m, n] = self.jaccard[n, m]

                    self.sizes[n,m] = (len(
                        graph_genomes[u_accession].representatives
                    ) + len(
                        graph_genomes[v_accession].representatives
                    ))/2
                    self.sizes[m,n] = self.sizes[n,m]

                    for k in self.distances:
                        self.distances[k][n, m] = dist_f.distances[k]
                        self.distances[k][m, n] = self.distances[k][n, m]

        if self.__learn:
            self.weights = self.__learn_weights()
            affinity = np.zeros_like(self.distances['mobile element'])
            for k in self.distances:
                if k != 'recombination':
                    affinity += self.distances[k]*self.weights[k]/ (
                        len(graph_genomes[u_accession]) + len(graph_genomes[v_accession])
                    )
            affinity += (self.distances['recombination']>0)*self.weights['recombination']*np.max(affinity)*self.jaccard
            
                
        return affinity

    def __learn_weights(
        self
    ):
        if self.learning_method == "ratio":

            totals = np.zeros_like(self.distances['mobile element'])
            for k in self.distances:
                totals += self.distances[k]*self.weights[k]

            new_weights = deepcopy(self.weights)
            
            for k in self.weights:
                new_weights[k] = np.sum(totals)/np.sum(self.distances[k])

            return new_weights
        
        elif self.learning_method == "precomputed":
            return self.weights

        elif self.learning_method == "max-likelihood":
            model = pm.Model()

            with model:

                n_samples = len(self.distances['mobile element'].flatten())
                n_classes = len(self.distances) - 1
                # Priors
                ts = pm.Uniform('ts', lower=0, upper=100, shape=n_samples)
                w = pm.HalfNormal('w', sigma=10, shape=n_classes)

                # Poisson's lambda
                lam_reversal = ts*w[0]*self.sizes.flatten()
                lam_duplication = ts*w[1]*self.sizes.flatten()
                lam_mobile = ts*w[2]*self.sizes.flatten()
                lam_transposition = ts*w[3]*self.sizes.flatten()
                lam_mutation = ts*w[4]*self.sizes.flatten()

                # Likelihood
                dist_reversal = pm.Poisson('dist_reversal', mu = lam_reversal, observed = self.distances['reversal'].flatten())
                dist_duplication = pm.Poisson('dist_duplication', mu = lam_duplication, observed = self.distances['duplicate'].flatten())
                dist_mobile = pm.Poisson('dist_mobile', mu = lam_mobile, observed = self.distances['mobile element'].flatten())
                dist_transposition = pm.Poisson('dist_transposition', mu = lam_transposition, observed = self.distances['transposition'].flatten())
                dist_mutation = pm.Poisson('dist_mutation', mu = lam_mutation, observed = self.distances['mutation'].flatten())

                _ = pm.sample(draws = 1)
                map_estimation = pm.find_MAP()

            return {
                'reversal': 1/(map_estimation['w'][0]+1e-10),
                'duplicate': 1/(map_estimation['w'][1]+1e-10),
                'mobile element': 1/(map_estimation['w'][2]+1e-10),
                'transposition': 1/(map_estimation['w'][3]+1e-10),
                'mutation': 1/(map_estimation['w'][4]+1e-10),
                'recombination': self.HR_DIST
            }

    def fit(
        self,
        dist_f: Callable
    ):
        '''
        Parameters
        ----------

        dist_f: Callable.
            Distance function taking two graph genomes (GraphGenome objects) as input
            and returning a number (distance between the two genomes).
        '''

        self.affinity = self.__build_affinity_matrix(dist_f)
        self.clusters = self.clustering_method.fit_predict(self.affinity)

    def fit_predict(
        self,
        dist_f: Callable
    ):
        '''
        Parameters
        ----------

        dist_f: Callable.
            Distance function taking two graph genomes (GraphGenome objects) as input
            and returning a number (distance between the two genomes).

        Returns
        -------

        clusters: list
            Cluster labels for each point.
        '''

        self.fit(dist_f)
        return self.clusters

    def get_clusters_as_series(
        self
    ) -> pd.Series:

        return pd.Series(
            self.clusters,
            self.accessions,
            name = 'Cluster'
        )

    def get_affinity_as_frame(
        self
    ) -> pd.DataFrame:

        return pd.DataFrame(
            self.affinity,
            index = self.accessions,
            columns = self.accessions
        )

    def plot_dendogram(
        self,
        figsize = (10,7)
    ):
        '''
        Adapted from https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html.
        '''
        # create the counts of samples under each node
        counts = np.zeros(self.clustering_method.children_.shape[0])
        n_samples = len(self.clustering_method.labels_)
        for i, merge in enumerate(self.clustering_method.children_):
            current_count = 0
            for child_idx in merge:
                if child_idx < n_samples:
                    current_count += 1  # leaf node
                else:
                    current_count += counts[child_idx - n_samples]
            counts[i] = current_count

        self.linkage_matrix = np.column_stack(
            [
                self.clustering_method.children_, 
                self.clustering_method.distances_, 
                counts
            ]
        ).astype(float)

        # Plot the corresponding dendrogram
        with plt.style.context('ggplot'):

            fig, ax = plt.subplots(
                1,
                1,
                figsize = figsize,
                constrained_layout = True
            )
            dendrogram(
                self.linkage_matrix,
                labels = self.accessions,
                distance_sort = 'ascending',
                leaf_label_func = lambda x: ' '.join(
                    [
                        self.accessions[x],
                        f'({self.clusters[x]})'
                    ]
                ),
                leaf_rotation = 90,
                ax = ax
            )

            ax.set_xlabel('Plasmid ID (Cluster)')
            ax.set_ylabel('Syntenic Distance')
            ax.grid(False)
            return ax

class Pangenome:
    
    def __init__(
        self,
        genomes: Iterable,
        path_to_protein_names: str,
        figsize: Tuple = ('1200px', '900px'),
        path: Union[str, None] = None
    ):
        self.genomes = genomes
        self.__w = figsize[0]
        self.__h = figsize[1]
        self.path = path
        self.colors = cm.rainbow(np.linspace(0, 1, len(self.genomes)))
        self.colors = [self.__rbga_to_hex(x) for x in self.colors]
        self.protein_names = self.__read_protein_names(path_to_protein_names)
        self.AMR_COLOR = '#0d9b18'
        self.SIZES = {'node': 2, 'edges': 10}

    def __read_protein_names(
        self,
        path: str
    ):
        with open(path, 'r') as instream:
            content = instream.read().split('\"')
        content[0] = content[0].split('\n')[-2]
        content = content[:-1]

        protein_ids = [content[x].replace('\n', '') for x in range(0, len(content), 2)]
        protein_names = [content[x].replace('\n', ' ') for x in range(1, len(content), 2)]

        unique_ids = np.unique(protein_ids)
        values = []
        indices = []
        for id_, name_ in zip(
            protein_ids,
            protein_names
        ):
            if id_ not in indices:
                values.append(name_)
                indices.append(id_)

        return pd.Series(
            values,
            indices,
            name = 'Protein Names'
        )

    def __light_build(
        self
    ):
        '''
        Iterates over all genomes, finding all genes and edges to add to the network.
        For each gene, it finds the original (Prokka-derived, before CD-HIT clustering)
        gene annotations, so it can display them in the gene edges.

        Unlike __build(), only one edge is added per gene.
        '''
        graph_dict = {
            'nodes': {},
            'edges': {'genes': {}, 'synteny': []}
        }
        added_gene_nodes = []
        added_gene_edges = []
        amr_genes = []

        # Find all the nodes. Add them to the graph dict
        for genome in self.genomes:
            for node, data_ in genome.genome.nodes(data=True):
                if data_['is_amr']:
                    color = self.AMR_COLOR
                    amr_genes.append(node.split()[0])
                else:
                    color = 'black'

                if node not in added_gene_nodes:
                    graph_dict['nodes'][node] = {
                        'id': node,
                        'title': node,
                        'color': color,
                        'size': 1
                    }
                    added_gene_nodes.append(node)

        for genome, color in zip(self.genomes, self.colors):
            for e1, e2, data_ in genome.genome.edges(data=True):
                if data_['color'] == 'gene':
                    if e1.split()[0] in amr_genes:
                        color_ = self.AMR_COLOR
                    else:
                        color_ = 'black'
                    
                    if (e1, e2) not in added_gene_edges:
                        graph_dict['edges']['genes'][e1.split()[0]] = {
                            'from': e1,
                            'to': e2,
                            'color': color_,
                            'title': e1.split()[0] + '\n' + ' | '.join(
                                [genome.id, data_['title']]
                            ),
                            'size': 1
                        }
                        added_gene_edges.append((e1, e2))
                        added_gene_edges.append((e2, e1))

                    else:
                        # If the edge was already added, add the original annotations
                        # to the title and increment the size by one
                        graph_dict['edges']['genes'][e1.split()[0]]['size'] += 1
                        graph_dict['edges']['genes'][e1.split()[0]]['title'] = '\n'.join(
                            [
                                graph_dict['edges']['genes'][e1.split()[0]]['title'],
                                ' | '.join([genome.id, data_['title']])
                            ]
                        )
                else:
                    # Synteny edges
                    graph_dict['edges']['synteny'].append(
                        {
                            'from': e1,
                            'to': e2,
                            'title': genome.id,
                            'color': color,
                            'size': 1
                        }
                    )

        # Build pyVis Network
        self.G = Network(
            height = self.__h,
            width = self.__w,
            bgcolor = '#FFFFFF',
            directed = True
        )

        for node in graph_dict['nodes'].values():
            self.G.add_node(
                node['id'],
                color = node['color'],
                size = self.SIZES['node'] * node['size'],
                title = node['title']
            )

        for edge in graph_dict['edges']['genes'].values():
            self.G.add_edge(
                edge['from'],
                edge['to'],
                color = edge['color'],
                title = edge['title'],
                arrows = 'none',
                width = edge['size'] * self.SIZES['edges']
            )
        for edge in graph_dict['edges']['synteny']:
            self.G.add_edge(
                edge['from'],
                edge['to'],
                color = edge['color'],
                title = edge['title'],
                arrows = 'none',
                width = edge['size'] * self.SIZES['edges']
            )

        self.G.toggle_physics(False)
        self.G.set_edge_smooth('dynamic')
        self.G.show_buttons(['physics', 'nodes'])
        self.G.force_atlas_2based(
            gravity = -50,
            central_gravity = 0.005,
            spring_length = 20,
            spring_strength = 1.,
            damping = .5
        )

        return self.G

    def __build(self):
        
        self.G = Network(
            height = self.__h,
            width = self.__w,
            bgcolor = '#FFFFFF',
            directed = True
        )

        added_gene_edges = []
        added_gene_title = []
        amr_genes = []

        for genome in self.genomes:
            for n, d in genome.genome.nodes(data=True):
                if d['is_amr']:
                    amr_genes.append(n.split()[0])

        for genome, color in zip(self.genomes, self.colors):

            for node, data_ in genome.genome.nodes(data=True):

                if node.split()[0] in amr_genes:
                    node_color = self.AMR_COLOR
                else:
                    node_color = 'black'

                self.G.add_node(
                    node,
                    color = node_color,
                    size = 10,
                    title = node
                )
            
            for e1, e2, d in genome.genome.edges(data=True):
                if d['color']=='gene':
                    # Is the gene an AMR gene? If yes, color green
                    if e1.split()[0] in amr_genes:
                        edge_color = self.AMR_COLOR
                    else:
                        edge_color = 'black'

                    # Each line in the gene edge hover box takes the form [Accession ID] | [Gene annotation]
                    title = ' | '.join(
                        [genome.id, e1.split()[0], d['title']]
                    )
                    # Add the edge to the already added gene edges
                    if (e1, e2) not in added_gene_edges:
                        self.G.add_edge(
                            e1,
                            e2,
                            color = edge_color,
                            width = 22,
                            title = title,
                            arrows = 'none'
                        )
                    else:
                        titles = np.array(added_gene_title)[
                            np.all(np.array(added_gene_edges) == (e1, e2), axis=-1)
                        ]
                        self.G.add_edge(
                            e1,
                            e2,
                            color = edge_color,
                            width = 22*len(titles),
                            title = '\n'.join(np.hstack([titles, title])),
                            arrows = 'none'
                        )
                    added_gene_edges.append((e1, e2))
                    added_gene_edges.append((e2, e1))
                    for _ in range(2):
                        added_gene_title.append(title)

                elif d['color'] != 'gene':
                    
                    self.G.add_edge(
                        e1,
                        e2,
                        color = color,
                        width = 10,
                        title = genome.id,
                        arrows = 'none'
                    )
                    
            
        self.G.toggle_physics(False)
        self.G.set_edge_smooth('dynamic')
        self.G.show_buttons(['physics', 'nodes'])
        self.G.force_atlas_2based(
            gravity = -50,
            central_gravity = 0.005,
            spring_length = 20,
            spring_strength = 1.,
            damping = .5
        )

    def show(
        self, 
        light_build = False, 
        with_duplicates = False, 
        hidden = False
    ):
        
        if with_duplicates:
            self.show_original_sequences(self.genomes, self.path)
        elif light_build:
            self.__light_build()
        else:
            self.__build()
        if self.path is None:
            path_ = os.path.join(
                os.getcwd(),
                'Pangenome-Graph-tmp.html'
            )
        else:
            path_ = self.path

        if not hidden:
            print(path_)

            self.G.show(path_)

            webbrowser.open(
                'file://///wsl.localhost/Ubuntu'+path_,
                new = 1
            )

            sleep(5)

            if self.path is None:
                os.remove(path_)
        
        return self.G

    def __rbga_to_hex(
        self,
        rgba_:  Iterable
    ) -> str:

        rgba = rgba_
        # Convert to [0, 255]
        rgba *= 255
        rgba = rgba[:-1].astype(int)
        rgba = np.minimum(255, rgba)

        return "#{0:02x}{1:02x}{2:02x}".format(*rgba)

    def to_networkx(
        self
    ):

        self.nx_graph = nx.MultiGraph()

        self.nx_graph.add_nodes_from(
            [
                (x['label'], {k: x[k] for k in ['color', 'size', 'id']})
                for x in self.G.nodes
            ]
        )
        self.nx_graph.add_edges_from(
            [
                (
                    x['from'],
                    x['to'],
                    {k: x[k] for k in ['color', 'width', 'title', 'arrows']}
                )
                for x in self.G.edges
                if x['color'] != 'black' and x['color'] != '#0d9b18'
            ]
        )

        # Add gene edges. These are not inherited from the pyVis
        # Network so as to only have a single edge per gene
        genes = [x['label'].split()[0] for x in self.G.nodes if x['color']=='black']
        amr_genes = [x['label'].split()[0] for x in self.G.nodes if x['color']!='black']
        genes = np.unique(genes)
        amr_genes = np.unique(amr_genes)

        self.nx_graph.add_edges_from(
            [
                (' '.join([x, 'Head']), ' '.join([x, 'Tail']))
                for x in genes
            ],
            color = 'black',
            width = 22,
            arrows = 'none'
        )
        self.nx_graph.add_edges_from(
            [
                (' '.join([x, 'Head']), ' '.join([x, 'Tail']))
                for x in amr_genes
            ],
            color = '#0d9b18',
            width = 22,
            arrows = 'none'
        )

        return self.nx_graph

    def show_original_sequences(
        self,
        genomes: Iterable,
        path_: str
    ) -> Network:
        '''
        Shows the original plasmid sequences, with duplicates, by building the pangenome from
        the GFF annotation files instead of the GraphGenome object graph.
        '''

        graph = Network(
            directed=True,
            width = '1200px',
            height = '900px',
            bgcolor = '#FFFFFF'
        )

        added_genes = []
        genes = {}
        for plasmid in genomes:
            # Get sequence information from GFF3 files, including duplicates
            sequence = plasmid.get_gene_sizes(include_duplicates=True).index.to_numpy()

            # Set the index of the genome.representatives DataFrame as the representative
            # proteins to aid querying
            protein_info = deepcopy(plasmid.representatives)
            protein_info = protein_info.set_index('Representative', drop=True)
            protein_info = pd.concat(
                [
                    protein_info,
                    pd.Series(
                        plasmid.representatives.index.to_numpy(),
                        index = protein_info.index,
                        name = 'Protein ID'
                    )
                ],
                axis = 'columns'
            )

            for gene in sequence:

                if protein_info.loc[gene]['Is AMR']: color_ = 'green'
                else: color_ = 'black'

                if gene not in added_genes:
                    added_genes.append(gene)
                    genes[gene] = {
                        'title': f'{gene}\nPlasmid {plasmid.id} | Gene {protein_info.loc[gene]["Protein ID"]} | {protein_info.loc[gene]["Product"]}',
                        'label': gene,
                        'color': color_,
                        'size': 1
                    }
                else:
                    # Append the product to the title
                    genes[gene]['title'] = '\n'.join(
                        [
                            genes[gene]['title'], 
                            f'Plasmid {plasmid.id} | Gene {protein_info.loc[gene]["Protein ID"]} | {protein_info.loc[gene]["Product"]}'
                        ]
                    )
                    genes[gene]['size'] += 1
            
        # Add gene nodes to the graph, as well as gene edges
        for gene, gene_info in genes.items():
            for suffix in [' Head', ' Tail']:
                graph.add_node(
                    gene + suffix,
                    label = gene_info['label'] + suffix,
                    size = 5 + gene_info['size']*0.5,
                    title = gene_info['title'],
                    color = gene_info['color']
                )
            graph.add_edge(
                gene+' Head',
                gene+' Tail',
                arrows = 'none',
                color = gene_info['color'],
                width = 10 + gene_info['size']*5,
                title = gene_info['title']
            )
        
        colors = cm.rainbow(np.linspace(0,1,len(genomes)))
        for plasmid, color in zip(genomes, colors):
            # Get sequence information from GFF3 files, including duplicates
            sequence = plasmid.get_gene_sizes(include_duplicates=True).index.to_list()
            # Set the index of the genome.representatives DataFrame as the representative
            # proteins to aid querying
            protein_info = deepcopy(plasmid.representatives).set_index('Representative', drop=True)
            protein_info['Protein ID'] = plasmid.representatives.index.to_numpy()

            for e1, e2 in zip(sequence, sequence[1:]+[sequence[0]]):

                # Set edge direction according to strand information
                if protein_info.loc[e1]['Strands'] > 0:
                    e1_suffix = ' Tail'
                else:
                    e1_suffix = ' Head'
                if protein_info.loc[e2]['Strands'] > 0:
                    e2_suffix = ' Head'
                else:
                    e2_suffix = ' Tail'
                # Add synteny edges
                graph.add_edge(
                    e1 + e1_suffix, e2 + e2_suffix,
                    color = self.__rbga_to_hex(deepcopy(color)),
                    width = 10,
                    arrows = 'none',
                    title = plasmid.id
                )

        graph.show_buttons(['physics', 'nodes'])
        graph.toggle_physics(False)
        graph.set_edge_smooth('dynamic')
        graph.force_atlas_2based(
            central_gravity = 0.005,
            spring_length = 50,
            spring_strength = 1
        )
        self.G = graph

        return graph


def plot_pangenome(
    ids: Iterable,
    tmp_path: str,
    annotations_path: str,
    with_duplicates = True,
    hidden = False
) -> Pangenome:
    pangenome = Pangenome(
        [
            GraphGenome(x, annotations_path, os.path.join(tmp_path, 'protein_cluster_membership.tsv'),
                os.path.join(tmp_path, 'amrfinder_output.tsv'))
            for x in ids
        ],
        path_to_protein_names = os.path.join(tmp_path, 'Protein Clusters.tsv'),
        path = None
    )
    pangenome.show(with_duplicates=with_duplicates, hidden = hidden)
    print(f'Showing pangenomes for plasmids {ids}.')

    return pangenome
    