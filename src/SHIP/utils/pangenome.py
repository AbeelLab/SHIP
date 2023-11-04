'''
Class to plot panplasmidomes.
'''

from time import sleep
import webbrowser
import networkx as nx
import pandas as pd
import numpy as np
from dataclasses import dataclass
import os 
from scipy.stats import kruskal
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from pyvis.network import Network
from typing import Iterable
from copy import deepcopy
from utils.interpretation import get_protein_function

@dataclass
class PangenomeAnalysis:
    '''
    Computes several stats from a pangenome graph (without strand information) obtained
    using a GeneMotifs or MotifFinder object. This graph corresponds to the pangenome
    graph of plasmids in a Jaccard megacluster, without edges present in only one plasmid.
    '''
    graph: nx.Graph
    motif_graph: Network
    phylo_config: dict
    
    def degree_report(
        self,
        show_centrality: bool = True,
        show_communities: bool = True
    ):
        '''
        Compares the degree of nodes in different regions of the filtered pangenome graph.

        Defines communities in the filtered pangenome, using the Louvain algorithm. The degree of 
        each node in a community is computed and these communities are compared using Kruskal-Wallis
        H-test to confirm the difference in median degree.

        This allows for the identification of plasmid regions showing more frequent genomic 
        alterations, quantified by the node degrees.
        '''
        partitions = nx.algorithms.community.louvain_communities(
            self.graph,
            resolution = 1
        )
        self.community_annotations = self.__community_annotation_content(partitions)
        degrees = {
            n: [nx.degree(self.graph, node) for node in partition]
            for n, partition in enumerate(partitions)
            if len(partition) > 4
        }
        self.mean_degrees = {k: np.mean(v) for k, v in degrees.items()}
        _, pval = kruskal(*list(degrees.values()))
        degrees = {
            n: [nx.degree(self.graph, node) for node in partition]
            for n, partition in enumerate(partitions)
        }
        self.mean_degrees = {k: np.mean(v) for k, v in degrees.items()}

        print(f'Found {len(partitions)} communities. The differences in degree are {"not "*(pval>.05)}statistically significant (p-value of {pval:.5}, Kruskal-Wallis).')

        if show_communities:
            colors, communities = {}, {}
            colormap = cm.inferno
            norm = Normalize(np.min(list(self.mean_degrees.values())), np.max(list(self.mean_degrees.values())))
            net = Network(
                height = '900px',
                width = '1200px'
            )
            for n, community in enumerate(partitions):
                for node in community:
                    colors[node] = np.array(colormap(norm(self.mean_degrees[n])))
                    communities[node] = n

            for node in self.motif_graph.nodes:
                net.add_node(
                    node['id'],
                    size = node['size'],
                    title = str(communities[node['id']]) + '\n' + node['title'],
                    label = node['label'],
                    color = self.__rgba_to_hex(colors[node['id']])
                )
            for edge in self.motif_graph.edges:
                net.add_edge(
                    edge['from'],
                    edge['to'],
                    title = edge['title'],
                    color = 'gray',
                    value = edge['value']
                )

            path_ = os.path.join(
                os.getcwd(),
                'Paths-Table-tmp.html'
            )
            net.toggle_physics(False)
            net.show_buttons('physics')
            net.force_atlas_2based(
                central_gravity = 0.005,
                spring_length = 20,
                spring_strength = 0.6,
                overlap = 0.9
            )
            net.show(path_)
            webbrowser.open(
                'file://///wsl.localhost/Ubuntu'+path_,
                new = 1
            )
            sleep(5)

            os.remove(path_)

        if show_centrality:
            centrality = nx.degree_centrality(self.graph)
            net = Network(
                height = '900px',
                width = '1200px'
            )

            # Maps the centrality values to [0, 1]
            norm = Normalize(min(centrality.values()), max(centrality.values()))
            for node in self.motif_graph.nodes:
                net.add_node(
                    node['id'],
                    size = node['size'],
                    title = node['title'],
                    label = node['label'],
                    color = self.__rgba_to_hex(
                        np.array(
                            cm.inferno(norm(centrality[node['id']]))
                        )
                    )
                )
            for edge in self.motif_graph.edges:
                net.add_edge(
                    edge['from'],
                    edge['to'],
                    title = edge['title'],
                    color = 'gray',
                    value = edge['value']
                )

            path_ = os.path.join(
                os.getcwd(),
                'Paths-Table-tmp.html'
            )
            net.toggle_physics(False)
            net.show_buttons('physics')
            net.force_atlas_2based(
                central_gravity = 0.005,
                spring_length = 20,
                spring_strength = 0.6,
                overlap = 0.9
            )
            net.show(path_)
            webbrowser.open(
                'file://///wsl.localhost/Ubuntu'+path_,
                new = 1
            )
            sleep(5)

            os.remove(path_)

        return self.mean_degrees, self.community_annotations

    def get_component_stats(
        self
    ):
        components = {}
        components['nodes'] = [x for x in nx.connected_components(self.graph)]
        components['number'] = len(components['nodes'])
        components['sizes'] = [len(x) for x in components['nodes']]

        return components


    def __community_annotation_content(
        self,
        communities: Iterable
    ) -> pd.DataFrame:
        '''
        Finds all gene annotations present in a community given by node_ids.
        '''
        for n, node_ids in enumerate(communities):
            annotations = get_protein_function(
                node_ids,
                self.phylo_config['paths']['protein-names']
            )
            # Remove inverted commas from the annotations. Drop the final accession
            annotations = annotations.apply(lambda x: ' '.join(x[1:-1].split()[:-1]))
            # Count the number of instances of each annotation
            annotations, counts = np.unique(annotations, return_counts=True)
            community_df = pd.DataFrame(
                [
                    [
                        n,
                        annotation,
                        count
                    ]
                    for annotation, count in zip(annotations, counts)
                ],
                columns = ['Community', 'Annotation', 'Count']
            )
            if n == 0:
                # Construct a DataFrame with all communities
                annotations_df = community_df
            else:
                annotations_df = pd.concat(
                    [annotations_df, community_df]
                )

        return annotations_df

    def __rgba_to_hex(
        self,
        rgba_: Iterable
    ) -> str:

        rgba = deepcopy(rgba_)
        # Convert to [0, 255]
        rgba *= 255
        rgba = rgba[:-1].astype(int)
        rgba = np.minimum(255, rgba)

        return "#{0:02x}{1:02x}{2:02x}".format(*rgba)
