'''
Classes to find conserved regions with AMR genes in plasmids.
'''
from time import sleep
import warnings
import webbrowser
import numpy as np
import pandas as pd
import os
import networkx as nx
from pyvis.network import Network
from copy import deepcopy
from typing import Iterable, Union
from utils.phylogeny import Pangenome
import joblib
from tqdm import tqdm
from utils.plasmid_typing import LabeledNetwork
from utils.phylogeny import Pangenome, plot_pangenome, SimpleDendogram
import matplotlib.cm as cm
from utils.files import find_annotation_paths
from BCBio.GFF import parse
from scipy.stats import kruskal
from matplotlib.colors import Normalize
from utils.interpretation import get_protein_function

class GeneMotifs:

    def __init__(
        self,
        phylo: SimpleDendogram,
        data_config_: dict,
        phylo_config_: dict,
        *,
        min_n_edges: int = 5,
        graph_gene_size: float = 60.0,
        graph_edge_size: float = 10.0
    ):

        self.min_n_edges = min_n_edges
        self.graph_gene_size = graph_gene_size
        self.graph_edge_size = graph_edge_size
        self.data_config = data_config_
        self.phylo_config = phylo_config_
        self.phylo = phylo

    def from_ids(
        self,
        pangenome: Pangenome,
        multigraph: bool = False
    ):
        '''
        Builds a graph containing sets of neighboring genes occuring multiple
        times across a pangenome. This allows for the identification of common
        motifs across distant plasmids. This graph will include duplicate genes,
        unlike the one obtained with from_pangenome(). Because of this, no prior
        call to Pangenome.show() is needed.
        
        In this graph, nodes represent genes and edges link genes present 
        consecutively in at least min_n_edges plasmids. The nodes are weighted 
        by the number of plasmids in which they appear consecutively. Strand 
        information is not taken into account.

        Parameters
        ----------

        pangenome: Pangenome
            Pangenome object for the set of plasmids in which to search for motifs.

        Returns
        -------

        motif_graph: pyVis Network
            Graph containing the motif graph in a pyVis Network object.

        See Also
        --------

        GeneMotifs.show_motifs(): Displays the motif graph in a browser tab.
        GeneMotifs.from_pangenome(): Builds a motif graph without gene duplicates,
            from a Pangenome object.
        '''
        self.__constructed_w_pangenome = False
        self.pangenome = pangenome
        self.protein_names = pangenome.protein_names
        self.protein_names = self.protein_names.apply(lambda x: ' '.join(x.split()[:-1]))
        # Plasmid accession numbers
        ids = [x.id for x in pangenome.genomes]
        # Get GFF3 annotation paths
        annotation_paths = find_annotation_paths(
            ids, 
            pangenome.genomes[0].path_to_annotations,
            '.gff'
        )

        # Get AMR info about the genes
        self.is_amr = pd.Series(
            [x[-1]['color']!='black' for x in pangenome.to_networkx().nodes(data=True)],
            index = [x.split()[0] for x in pangenome.to_networkx().nodes()],
            name = 'Is AMR'
        )
        self.is_amr = self.is_amr[self.is_amr.index.duplicated(keep = 'first')]

        protein_clusters = pd.read_csv(
            self.phylo_config['paths']['representative-proteins'], 
            sep='\t', 
            index_col = 0
        )['Representative']

        graph_dict = {'nodes': {}, 'edges': {}}
        added_nodes, added_edges = [], []
        for path_ in annotation_paths:
            with open(path_) as file:
                features = next(parse(file)).features
            for feature in features:
                try:
                    gene = protein_clusters[feature.id]
                except KeyError:
                    gene = feature.id
                if gene not in added_nodes:
                    try:
                        if self.is_amr[gene]:
                            color_ = 'green'
                        else:
                            color_ = 'gray'
                    except KeyError:
                        color_ = 'gray'

                    graph_dict['nodes'][gene] = {
                        'id': gene,
                        'title': gene,
                        'label': feature.qualifiers['product'][0],
                        'color': color_
                    }
                    added_nodes.append(gene)
                else:
                    graph_dict['nodes'][gene]['label'] = '\n'.join(
                        [graph_dict['nodes'][gene]['label'], feature.qualifiers['product'][0]]
                    )
            
            for e1, e2 in zip(features, [features[-1]]+features[:-1]):
                try:
                    e1_ = protein_clusters[e1.id]
                except KeyError:
                    e1_ = e1.id
                try:
                    e2_ = protein_clusters[e2.id]
                except KeyError:
                    e2_ = e2.id
                edge = np.sort([e1_, e2_])
                if '-'.join(edge) not in added_edges:
                    added_edges.append('-'.join(edge))
                    graph_dict['edges'][tuple(edge)] = {
                        'from': edge[0],
                        'to': edge[1],
                        'color': 'gray',
                        'value': 1
                    }
                else:
                    graph_dict['edges'][tuple(edge)]['value'] += 1
        
        # Convert dict to pyvis
        if multigraph: directed = True
        else: directed = False
        graph = Network(
            height = '900px',
            width = '1200px',
            directed = directed
        )            
        for edge in graph_dict['edges'].values():
            if edge['value'] >= self.min_n_edges:
                for node_id in [edge['from'], edge['to']]:
                    node = graph_dict['nodes'][node_id]
                    graph.add_node(
                        node['id'],
                        color = node['color'],
                        title = node['label'],
                        label = node['title'],
                        size = 9
                    )
                graph.add_edge(
                    edge['from'],
                    edge['to'],
                    value = edge['value'],
                    color = edge['color'],
                    title = ' ',
                    label = ' ',
                    arrows = 'none'
                )

        self.motif_graph = graph

    def from_pangenome(
        self,
        pangenome: Pangenome
    ) -> Network:

        '''
        Builds a graph containing sets of neighboring genes occuring multiple
        times across a pangenome. This allows for the identification of common
        motifs across distant plasmids.
        
        In this graph, nodes represent genes and edges link genes present 
        consecutively in at least min_n_edges plasmids. The nodes are weighted 
        by the number of plasmids in which they appear consecutively. Strand 
        information is not taken into account.

        Parameters
        ----------

        pangenome: Pangenome
            Pangenome object for the set of plasmids in which to search for motifs.
            The Pangenome method show() should be called before passing it to this
            class method.

        Returns
        -------

        motif_graph: pyVis Network
            Graph containing the motif graph in a pyVis Network object.

        See Also
        --------

        GeneMotifs.show_motifs(): Displays the motif graph in a browser tab.

        '''
        self.__constructed_w_pangenome = True
        self.pangenome = deepcopy(pangenome)
        pangenome_nx = pangenome.to_networkx()
        self.protein_names = pangenome.protein_names
        self.protein_names = self.protein_names.apply(lambda x: ' '.join(x.split()[:-1]))

        # Graph with a single node per gene (no strand information)
        self.merged_graph = nx.MultiGraph()
        for from_, to, data in pangenome_nx.edges(data=True):
            if data['color'] not in ['black', '#0d9b18'] and from_.split()[0] != to.split()[0]:
                self.merged_graph.add_edge(
                    from_.split()[0],
                    to.split()[0]
                )

        deg = [x for x in self.merged_graph.degree()]
        gene_names, gene_degrees = np.array(
            [x[0] for x in deg]
        ), np.array(
            [x[1] for x in deg]
        ) 
        self.is_amr = pd.Series(
            [x[-1]['color']!='black' for x in pangenome_nx.nodes(data=True)],
            index = [x.split()[0] for x in pangenome_nx.nodes()],
            name = 'Is AMR'
        )
        self.is_amr = self.is_amr[self.is_amr.index.duplicated(keep = 'first')]

        idx = np.argsort(gene_degrees)[::-1]
        gene_degrees = gene_degrees[idx]
        gene_names = gene_names[idx]

        self.degrees = pd.Series(
            gene_degrees,
            index = gene_names,
            name = 'Degree'
        ).apply(lambda x: x/2)

        self.motif_graph = Network(
            height = '900px',
            width = '1200px',
            directed = True
        )

        unique_edges, weights = np.unique(
            [
                np.sort([x[0], x[1]]) 
                for x in self.merged_graph.edges(data=True)
            ],
            return_counts = True,
            axis = 0
        )

        for edge, weight in zip(unique_edges, weights):
            normalized_weight = weight/(self.degrees[edge].sum())
            #if normalized_weight > THRESHOLD and weight > 2:
            # I don't think we should be aiming for specificity, so it should
            # be better to not normalize the number of co-occurences by the
            # frequency of both genes in the dataset
            if weight > self.min_n_edges:
                for gene in edge:
                    if self.is_amr[gene]: color = '#0d9b18'
                    else: color = 'gray'
                    self.motif_graph.add_node(
                        gene,
                        size = self.graph_gene_size/2,
                        color = color,
                        title = f'{gene} | {self.protein_names.loc[gene]}'
                    )
                self.motif_graph.add_edge(
                    edge[0], edge[1],
                    title = f'No. co-occurences: {weight}\nNormalized weight: {normalized_weight}',
                    value = float(weight),
                    arrows = 'none',
                    color = 'gray'
                )

        return self.motif_graph

    def show_motifs(
        self,
        path_: Union[str, None] = None
    ) -> str:
        '''
        Shows the graph of gene motifs in a browser window. Optionally, stores the 
        HTML file in the directory specified in path_. In this graph, nodes represent
        genes and edges link genes present consecutively in at least min_n_edges 
        plasmids. The nodes are weighted by the number of plasmids in which they
        appear consecutively. Strand information is not taken into account.

        Parameters
        ----------

        path_: str or None. Default is None
            Path in which to store the resulting HTML pyVis file. If None, the file is
            temporarily stored in the current working directory, and deleted after it
            is shown in a browser tab.

        Returns
        -------

            path_: str
                The path to the resulting HTML file.

        '''

        if path_ is None:
            path_ = os.path.join(
                            os.getcwd(),
                            'Pangenome-Graph-tmp.html'
                        )
          
        self.motif_graph.toggle_physics(False)
        self.motif_graph.show_buttons('physics')
        self.motif_graph.show(path_)
        webbrowser.open(
            'file://///wsl.localhost/Ubuntu'+path_,
            new = 1
        )
        sleep(5)

        if path_ is None:
            os.remove(path_)

        return path_

    def __get_genomes(
        self,
        path_to_representatives: str
    ) -> pd.Series:
        '''
        Returns the gene sequences of the pangenome plasmids as a Pandas Series.
        '''
        self.representative_sequences = joblib.load(path_to_representatives).loc[
            [x.id for x in self.pangenome.genomes]
        ]
        return self.representative_sequences

    def has_motif(
        self,
        query_sequence: list,
        motif: list
    ) -> bool:
        '''
        Checks if a gene motif is present in a given gene sequence. Returns true if the
        motif is found, and false otherwise.

        Parameters
        ----------

        query_sequence: List of str
            List containing the ordered gene sequence of the plasmid in which to check for
            the presenece of the motif.

        motif: List of str
            List contiaing the ordered genes that constitute the motif.

        Returns
        -------

        found: bool
        '''

        motif = self.filter_query_for_duplicates(motif) 
        # Start by searching if all genes in the motif are in the query sequence
        #if not all(x in query_sequence for x in motif): return False

        # Append motif-sized chunks from the end of the query to its start,
        # to get a circular list-like behaviour
        query = query_sequence[-len(motif):] + query_sequence
        reverse_motif = motif[::-1]
        for m in [motif, reverse_motif]:
            for i in range(len(query)-len(m)):
                if query[i:i+len(m)] == m: return True

        return False

    def get_motif_graph_as_nx(
        self
    ) -> nx.Graph:

        self.motif_nx = nx.MultiGraph()
        self.motif_nx.add_nodes_from(
            [
                (x['id'], x)
                for x in self.motif_graph.nodes
            ]
        )
        self.motif_nx.add_edges_from(
            [
                (x['from'], x['to'], {k: x[k] for k in ['title', 'value', 'arrows', 'color']})
                for x in self.motif_graph.edges
            ]
        )

        return self.motif_nx

    def get_representative_w_duplicates_removed(self):
        return {
            x.id: [
                xx 
                for xx in x.representatives_with_duplicates['Representative'].values
                if np.sum(x.representatives_with_duplicates['Representative'].values == xx) == 1 and not self.is_amr[xx]
            ]   
            for x in self.pangenome.genomes
        }

    def filter_query_for_duplicates(
        self,
        query: Iterable[str]
    ):
        all_genes = []
        for v in self.representative_sequences.values(): all_genes += v
        
        return [x for x in query if x in all_genes]
    
    def find_regions(
        self,
        query_gene: str,
        max_distance: int,
        min_distance: int = 5
    ):
        self.representative_sequences = {x.id: list(x.representatives_with_duplicates['Representative'].values) for x in self.pangenome.genomes}

        # Break each sequence into chunks of size k
        chunks_per_plasmid = {id_: [] for id_ in self.representative_sequences.keys()}
        for k in range(min_distance, max_distance+1):
            for id_, gene_seq in self.representative_sequences.items():
                i = 0
                if query_gene in gene_seq:
                    extend_gene_seq = gene_seq + gene_seq[:k]
                    while i < len(extend_gene_seq)-k:
                        chunk = extend_gene_seq[i:i+k]
                        if query_gene in chunk: chunks_per_plasmid[id_].append(chunk)
                        i += 1

        unique_chunks = []
        for v in chunks_per_plasmid.values(): unique_chunks += v
        unique_chunks = np.unique(unique_chunks)

        # Chunks present in > 1 plasmid and the ids of plasmids containing each chunk
        regions_out = []
        for chunk in unique_chunks:
            plasmids_w_chunk = []
            for k, v in chunks_per_plasmid.items():
                if chunk in v or chunk[::-1] in v:
                    plasmids_w_chunk.append(k)
            if len(plasmids_w_chunk) > 1:
                regions_out.append((chunk, plasmids_w_chunk))

        return regions_out
    
    def __call__(
        self,
        query_gene: str,
        distance_matrix_: pd.DataFrame,
        max_distance: int,
        min_plasmid_distance: float = 0.0,
        *,
        show: bool = True
    ):
        # Normalize the distance matrix
        distance_matrix = deepcopy(distance_matrix_) / np.max(distance_matrix_.values)
        self.query_graph = nx.MultiGraph()
        regions = self.find_regions(query_gene, max_distance, 5)
        self.n_paths = 0
        self.report = None
        for path, plasmids_with_motif in regions:
            avg_distance = np.mean(distance_matrix.loc[plasmids_with_motif][plasmids_with_motif].values)
            if avg_distance >= min_plasmid_distance:
                # Add path to the resulting graph, with weight equal to the average distance
                nx.add_path(self.query_graph, path, value = avg_distance, weight = avg_distance, path_n = self.n_paths)
                data = [[path, len(path), len(distance_matrix.loc[plasmids_with_motif][plasmids_with_motif]),
                    avg_distance, plasmids_with_motif,
                    [self.protein_names.loc[x] for x in path if 'integrase' in self.protein_names.loc[x]],
                    [self.protein_names.loc[x] for x in path if 'transposase' in self.protein_names.loc[x]]
                ]]
                if self.n_paths == 0:
                    self.report = pd.DataFrame(
                        data,
                        columns = ['Genes', 'Length', 'Number of Plasmids', 'Average Distance', 'Plasmids', 'Integrases', 'Transposases'],
                        index = [self.n_paths]
                    )
                else:
                    self.report = pd.concat([self.report, pd.DataFrame(data, columns = self.report.columns, index = [self.n_paths])])
                self.n_paths += 1
    '''
    def __call__(
        self,
        query_gene: str,
        distance_matrix_: pd.DataFrame,
        max_distance: int,
        min_plasmid_distance: float = 0.0,
        *,
        show: bool = True
    ):

        # Normalize the distance matrix
        distance_matrix = deepcopy(distance_matrix_) / np.max(distance_matrix_.values)
        self.query_graph = nx.MultiGraph()
        self.representative_sequences = self.get_representative_w_duplicates_removed()
        self.motif_nx = self.get_motif_graph_as_nx()
        # Iterate over connected components
        components = nx.connected_components(self.motif_nx)
        self.n_paths = 0
        for component in components:
            # Check if the query gene is in the component
            if query_gene in component:
                # Keep only the nodes with distance to the query gene <= d_MAX
                candidate_nodes = [
                    x for x in component
                    if nx.shortest_path_length(
                        self.motif_nx,
                        query_gene,
                        x
                    ) <= max_distance
                ]
                # Find all paths between candidate nodes, including the query gene

                for n, start in enumerate(candidate_nodes):
                    for end in candidate_nodes[n+1:]:
                        possible_paths = nx.all_simple_paths(self.motif_nx, start, end, cutoff=max_distance)
                        for path in possible_paths:
                            if query_gene in path:
                                # For each possible path (motif), compute the average distance of
                                # plasmids containing that motif
                                plasmids_with_motif = [
                                    plasmid
                                    for plasmid, sequence in self.representative_sequences.items()
                                    if self.has_motif(sequence, path)
                                ]
                                if len(plasmids_with_motif):
                                    avg_distance = np.mean(
                                        distance_matrix.loc[plasmids_with_motif][plasmids_with_motif].values
                                    )
                                    if avg_distance >= min_plasmid_distance:
                                        # Add path to the resulting graph, with weight equal to the average distance
                                        nx.add_path(
                                            self.query_graph, 
                                            path, 
                                            value = avg_distance, 
                                            weight = avg_distance, 
                                            path_n = self.n_paths
                                        )
                                        data = [[
                                            path,
                                            len(path),
                                            len(distance_matrix.loc[plasmids_with_motif][plasmids_with_motif]),
                                            avg_distance,
                                            plasmids_with_motif,
                                            [
                                                self.protein_names.loc[x] 
                                                for x in path 
                                                if 'integrase' in self.protein_names.loc[x]
                                            ],
                                            [
                                                self.protein_names.loc[x]
                                                for x in path
                                                if 'transposase' in self.protein_names.loc[x]
                                            ]
                                        ]]
                                        if self.n_paths == 0:
                                            self.report = pd.DataFrame(
                                                data,
                                                columns = ['Genes', 'Length', 'Number of Plasmids', 'Average Distance', 'Plasmids', 'Integrases', 'Transposases'],
                                                index = [self.n_paths]
                                            )
                                        else:
                                            self.report = pd.concat(
                                                [
                                                    self.report,
                                                    pd.DataFrame(
                                                        data,
                                                        columns = self.report.columns,
                                                        index = [self.n_paths]
                                                    )
                                                ]
                                            )
                                        self.n_paths += 1

        try:
            self.report = self.report.sort_values(['Average Distance', 'Length', 'Number of Plasmids'], ascending=False)
        except:
            self.report = None
        if show:
            self.report = self.__show_table(self.report)
        return self.query_graph
    '''

    def show_query_graph(
        self,
        path_: Union[str, None] = None
    ) -> str:
        '''
        Shows the resulting graph after the call to GeneMotifs() as a pyVis Network.
        '''

        colors = cm.rainbow(np.linspace(0, 1, self.n_paths))
        colors = [self.__rbga_to_hex(x) for x in colors]

        net = Network(
            directed=True,
            height = '900px',
            width='1200px'
        )
        for gene in self.query_graph.nodes:
            if self.is_amr[gene]: color = '#0d9b18'
            else: color = 'gray'
            net.add_node(
                gene,
                size = self.graph_gene_size/2,
                color = color,
                title = f'{gene} | {self.protein_names.loc[gene]}'
            )
        for edge in self.query_graph.edges(data=True):
            net.add_edge(
                edge[0], edge[1],
                title = f'Path Number {edge[2]["path_n"]}\nAverage plasmid distance: {edge[2]["weight"]:.5}',
                value = edge[2]['weight'],
                arrows = 'none',
                color = colors[edge[2]['path_n']]
            )
        
        if path_ is None:
            path_ = os.path.join(
                os.getcwd(),
                'Pangenome-Graph-tmp.html'
            )
          
        net.toggle_physics(False)
        net.show_buttons('physics')
        net.show(path_)
        webbrowser.open(
            'file://///wsl.localhost/Ubuntu'+path_,
            new = 1
        )
        sleep(5)

        if path_ is None:
            os.remove(path_)

        return path_

    def __show_table(
        self,
        df: pd.DataFrame
    ):

        filtered_df = df.drop(['Genes', 'Plasmids'], axis = 'columns')
        tmp_df = deepcopy(filtered_df)
        tmp_df.index = [f'{n_plasmids}-{avg_d}' for n_plasmids, avg_d in zip(tmp_df['Number of Plasmids'], tmp_df['Average Distance'])]
        idx = ~tmp_df.index.duplicated(keep = 'first')
        filtered_df = df.loc[idx].reset_index(drop=True)[['Average Distance', 'Length', 'Number of Plasmids', 'Integrases', 'Transposases', 'Genes', 'Plasmids']]

        html_table = filtered_df.to_html()
        document = f'''
<!DOCTYPE html>
<html>
    <body>'''+html_table+'</body>\n</html>'

        path_ = os.path.join(
            os.getcwd(),
            'Paths-Table-tmp.html'
        )
          
        with open(path_, 'x') as outstream:
            outstream.write(document)
        webbrowser.open(
            'file://///wsl.localhost/Ubuntu'+path_,
            new = 1
        )
        sleep(1)

        os.remove(path_)

        return filtered_df


    def __rbga_to_hex(
        self,
        rgba_: Iterable
    ) -> str:

        rgba = deepcopy(rgba_)
        # Convert to [0, 255]
        rgba *= 255
        rgba = rgba[:-1].astype(int)
        rgba = np.minimum(255, rgba)

        return "#{0:02x}{1:02x}{2:02x}".format(*rgba)                         

    def show_plasmids_from_path(
        self,
        path_number: int
    ) -> list:

        
        self.result_graph = plot_pangenome(
            self.report.loc[path_number]['Plasmids'],
            self.data_config,
            self.phylo_config
        )

        return self.report.loc[path_number]['Plasmids']

    def highlight_plasmids(
        self,
        plasmids: list,
        path: Union[None, str] = None,
        distance_threshold: float = .5,
        colors = None
    ):
        assert colors is None or len(colors) == 2, 'If specified, colors be length two.'
        network = LabeledNetwork(
            distance_threshold=distance_threshold, 
            colors=colors
        )
        network.from_distance(
            self.phylo.get_affinity_as_frame(),
            pd.Series(
                [(x in plasmids)*'With Motif'+(x not in plasmids)*'Without Motif' for x in self.phylo.get_affinity_as_frame().index],
                index = self.phylo.get_affinity_as_frame().index,
                name = 'Motif Presence'
            )
        )
        network.show(path)
        return path

    def show_path(
        self,
        path_number: int,
        path_: str = None
    ):
        path = self.report.loc[path_number]['Genes']
        self.motif_pyvis = Network(
            directed=False,
            height = '900px',
            width = '1200px'
        )
        for node in path:
            pangenome_node = [x for x in self.result_graph.G.nodes if x['id'].split()[0]==node][0]
            # Find gene edges and extract information
            gene_edge = [x for x in self.result_graph.G.edges if x['from'].split()[0] == node and x['to'].split()[0] == node][0]
            self.motif_pyvis.add_node(
                node,
                title = gene_edge['title'],
                label = gene_edge['title'].split('\n')[1].split(':')[-1],
                color = pangenome_node['color'],
                size = 7
            )
            
        for e1, e2 in zip(path[:-1], path[1:]):
            self.motif_pyvis.add_edge(
                e1,
                e2,
                color = 'gray',
                width = self.graph_edge_size
            )

        self.motif_pyvis.toggle_physics(False)
        self.motif_pyvis.set_edge_smooth('dynamic')
        self.motif_pyvis.show_buttons(['physics', 'nodes'])
        if path_ is None:
            path_ = os.path.join(
                            os.getcwd(),
                            'Motif-Graph-tmp.html'
                        )
        self.motif_pyvis.show(path_)
        webbrowser.open(
            'file://///wsl.localhost/Ubuntu'+path_,
            new = 1
        )
        sleep(5)

        if path_ is None:
            os.remove(path_)

        return self.motif_pyvis

    def show_original_sequences(
        self,
        selected_path_number: int,
        path_: str = None
    ) -> Network:
        '''
        Shows the selected motif in the original gene sequences, with duplicate genes.
        '''

        ids = self.report['Plasmids'].loc[selected_path_number]
        motif = self.report['Genes'].loc[selected_path_number]

        colors = cm.rainbow(np.linspace(0, 1, len(ids)))

        graph = Network(
            directed=True,
            width = '1200px',
            height = '900px',
            bgcolor = '#FFFFFF'
        )

        protein_clusters = pd.read_csv(
            self.phylo_config['paths']['representative-proteins'], 
            sep='\t', 
            index_col = 0
        )

        added_nodes = []
        for plasmid, color in zip(ids, colors):
            # Get sequence information from GFF3 files
            gff_path = find_annotation_paths([plasmid], self.data_config['paths']['annotations'], format='.gff')[0]
            with open(gff_path, 'r') as file:
                features = next(parse(file)).features
            original_sequence = [x.id for x in features]
            sequence = []
            # Find the representative protein (from CD-Hit)
            for protein in original_sequence:
                try:
                    sequence.append(protein_clusters.loc[protein]['Representative'])
                except KeyError:
                    sequence.append(protein)
            
            for e1, e2 in zip(sequence, sequence[1:]+[sequence[0]]):
                for node in [e1, e2]:
                    # Add node if missing from graph
                    if node not in added_nodes:
                        if node in motif:
                            node_color = 'red'
                        else:
                            node_color = 'gray'

                        try:
                            title = self.protein_names.loc[node]
                        except KeyError:
                            title = f'Could not find name for protein {node}.'

                        graph.add_node(
                            node,
                            color = node_color,
                            size = 9,
                            label = node,
                            title = title
                        )
                        added_nodes.append(node)

                # Add synteny edges
                graph.add_edge(
                    e1, e2,
                    color = self.__rbga_to_hex(deepcopy(color)),
                    width = self.graph_edge_size/2,
                    arrows = 'none',
                    title = plasmid
                )

        graph.show_buttons(['physics'])
        graph.toggle_physics(False)
        graph.set_edge_smooth('dynamic')
        if path_ is None:
            path_ = os.path.join(
                            os.getcwd(),
                            'Motif-Graph-w-Duplicates-tmp.html'
                        )
        graph.show(path_)
        webbrowser.open(
            'file://///wsl.localhost/Ubuntu'+path_,
            new = 1
        )
        sleep(5)

        if path_ is None:
            os.remove(path_)

        return graph

class MotifFinder(GeneMotifs):
    '''
    Wraps a GeneMotifs object to aid the exploration of AMR gene cassetes and
    transposable elements.
    '''
    def __init__(
        self, 
        phylo: SimpleDendogram, 
        data_config_: dict, 
        phylo_config_: dict,
        min_n_edges: int = 5, 
        salamzade: bool = False,
        *,
        graph_gene_size: float = 60, 
        graph_edge_size: float = 10,
        graph_max_distance: float = .4
    ):
        super().__init__(
            phylo, 
            data_config_, 
            phylo_config_,
            min_n_edges = min_n_edges,
            graph_gene_size = graph_gene_size,
            graph_edge_size = graph_edge_size
        )
        self.graph_max_distance = graph_max_distance
        self.salamzade = salamzade

    def __prepare_pangenome(
        self,
        ids: Iterable,
        salamzade: bool = False
    ):
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.pangenome = plot_pangenome(
                ids,
                self.data_config,
                self.phylo_config,
                with_duplicates = False,
                salamzade = salamzade
            )
        self.from_pangenome(self.pangenome)
    
    def __select_motifs(
        self
    ):
        print('\nShowing pangenome connected components.')
        self.show_motifs()

        options = [x['id'] for x in self.motif_graph.nodes]
        selected_id = input(
            'Select a gene id to find motifs centered around it: '
        )
        while selected_id not in options:
            selected_id = input(
                f'Invalid ID {selected_id}.Select a gene id to find motifs centered around it: '
            )
        return selected_id

    def __select_search_options(
        self
    ):
        valid_length = False
        while not valid_length:
            max_length = input(
                'Select a maximum motif length: '
            )
            try:
                max_length = int(max_length)
                if max_length <= 0:
                    print(f'Invalid length (non-positive): {max_length}.')
                else:
                    valid_length = True
            except:
                print(f'Invalid length. Got {max_length}.')

        valid_threshold = False
        while not valid_threshold:
            min_score = input(
                'Select a minimum average plasmid distance: '
            )
            try:
                min_score = float(min_score)
                if min_score < 0 or min_score > 1:
                    print(f'Invalid threshold. Must be in [0, 1]. Got {min_score}.')
                else:
                    valid_threshold = True
            except:
                print(f'Invalid threshold. Must be in [0, 1]. Got {min_score}.')

        return max_length, min_score

    def __select_path(
        self
    ):
        print(f'\nShowing motifs found for gene {self.gene}: ')
        
        options = self.report.index.to_list()
        selected_path = input(
            'Select a path number for further analysis: '
        )
        selected_path = int(selected_path)
        while selected_path not in options:
            selected_path = input(
                f'Invalid path {selected_path}.Select a path number for further analysis: '
            )
        return selected_path

    def __binary_input(
        self,
        msg: str
    ):
        while True:
            user_input = input(msg)
            user_input = user_input.lower()
            if user_input == 'y':
                return True
            elif user_input == 'n':
                return False
            else:
                print(f'Invalid input, must be Y or N. Got {user_input}.\n')

    def start(
        self,
        ids: Iterable,
        call_pangenome: bool = True
    ):
        # Allows searching for multiple paths
        continue_search = True

        self.ids = ids
        if call_pangenome: self.__prepare_pangenome(ids, self.salamzade)
        self.gene = self.__select_motifs()
        self.max_length, self.min_plasmid_distance = self.__select_search_options()
        self.result = self(
            self.gene, 
            self.phylo.get_affinity_as_frame(), 
            self.max_length, 
            self.min_plasmid_distance
        )
        self.show_query_graph()
        print(
            self.report[['Length', 'Number of Plasmids', 'Average Distance']]
        )
        while continue_search:
            self.selected_path = self.__select_path()
            self.plasmids_in_path = self.show_plasmids_from_path(self.selected_path)
            self.highlight_plasmids(
                self.plasmids_in_path,
                distance_threshold  = self.graph_max_distance,
                colors = np.vstack([
                    np.array([33/255, 111/255, 221/255, 1]),
                    np.array([242/255, 242/255, 242/255, 1]),
                ])
            )
            self.show_path(self.selected_path)
            self.show_original_sequences(
                self.selected_path
            )
            continue_search = self.__binary_input('Would you like to examine another path? (Y | N): ')
        
        if self.__binary_input('Would you like to examine another gene in the same cluster? (Y | N): '):
            self.start(
                ids,
                call_pangenome = False
            )

class BulkMotifFinder:

    def __init__(
        self,
        ids: Iterable,
        phylo: SimpleDendogram,
        data_config: dict,
        phylo_config: dict,
        *,
        min_n_edges: int = 2,
        salamzade: bool = False
    ):

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.pangenome = plot_pangenome(
                ids,
                data_config,
                phylo_config,
                with_duplicates = True, #
                hidden = True,
                salamzade=salamzade
            )
        self.ids = ids
        self.phylo = phylo
        self.data_config = data_config
        self.phylo_config = phylo_config
        self.min_n_edges = min_n_edges

        self.motif_finder = MotifFinder(
            self.phylo,
            self.data_config,
            self.phylo_config,
            self.min_n_edges,
            salamzade
        )

        #self.motif_finder.from_pangenome(self.pangenome)
        self.motif_finder.from_ids(self.pangenome)

    def search(
        self,
        min_distance: float,
        min_length: int,
        max_length: int,
        min_n_plasmids: int = 3
    ) -> pd.DataFrame:

        # Find AMR genes
        self.amr_genes = self.motif_finder.is_amr[self.motif_finder.is_amr].index.to_numpy()
        built = False
        self.report = None

        for n, gene in tqdm(enumerate(self.amr_genes), total = len(self.amr_genes)):

            self.motif_finder(
                gene,
                self.phylo.get_affinity_as_frame(),
                max_length,
                min_distance,
                show = False
            )

            report = self.motif_finder.report
            if report is not None:
                report = report[report['Length'] >= min_length]
                report = report[report['Number of Plasmids'] >= min_n_plasmids]

                # Add query gene to the report DataFrame
                report['Query Gene ID'] = [gene]*len(report)
                gene_annotation = self.motif_finder.protein_names[gene]
                if len(gene_annotation) == 0:
                    gene_annotation = 'Unknown'
                report['Query Gene Annotation'] = [gene_annotation]*len(report)
                report['Gene Annotations'] = [
                    list(self.motif_finder.protein_names[genes_in_path].values.flatten())
                    for genes_in_path in report['Genes']
                ]

                if not built:
                    self.report = report
                    built = True
                else:
                    self.report = pd.concat(
                        [
                            self.report,
                            report
                        ]
                    )
        return self.report

    def save_report(
        self,
        path: str
    ):

        self.report.to_csv(
            path,
            sep = '\t'
        )

def is_region_in_plasmid(
    plasmid_gene_seq: Iterable[str],
    query: Iterable[str]
) -> bool:
    '''
    Searches for an ordered list of genes/homolog groups, referenced by their ID, in 
    a plasmid gene sequence. Returns True if the query is found. The search is repeated
    with the query list inverted

    Parameters
    ----------
    plasmid_gene_seq: Iterable of str
        Gene or representative gene IDs in the plasmid, ordered.
    query: Iterable of str
        Genes or representative genes IDs to search for.

    Returns
    -------
    found: bool
    '''

    rev_query = query[::-1]
    query_len = len(query)

    for q in [query, rev_query]:

        ref_cursor_pos = 0
        query_cursor_pos = 0

        while ref_cursor_pos < len(plasmid_gene_seq):

            if plasmid_gene_seq[ref_cursor_pos] == q[query_cursor_pos]:
                query_cursor_pos += 1
            else:
                # Start from the begining of the query again
                query_cursor_pos = 0
            
            if query_cursor_pos >= query_len: return True

            ref_cursor_pos += 1

    return False