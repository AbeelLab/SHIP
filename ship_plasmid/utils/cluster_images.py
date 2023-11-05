'''
Classes for clustering analysis. Draws networks and plots with plasmid
network statistics.
'''
from copy import deepcopy
from itertools import accumulate
import os
from time import sleep
import webbrowser
import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Iterable, Any, Union, Tuple
from pyvis.network import Network
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import joblib
from scipy.cluster.hierarchy import dendrogram
from matplotlib.patches import Rectangle
from ship_plasmid.utils.phylogeny import SimpleDendogram
from ship_plasmid.utils.jaccard import jaccard
from scipy.stats import pearsonr

class GraphPlot:

    def __init__(
        self,
        X: Union[pd.DataFrame, None],
        clusters: pd.DataFrame,
        affinity: pd.DataFrame,
        info_path: str,
        amr_path: str
    ):
        '''
        
        Arguments
        ----------

        X: Pandas DataFrame
            Data in which clustering was performed.
        
        clusters: Pandas DataFrame
            DataFrame with clusters as indices and a column 'Plasmids' with all plasmids
            belonging to each cluster.

        affinity: Pandas DataFrame
            Affinity matrix.

        info_path: str
            Path to plasmid info.

        amr_path:
            Path to info about AMR gene hits from AMRFinder.
        '''

        assert type(clusters) == pd.DataFrame, 'Clusters must be in Pandas DataFrame format.'
        assert type(affinity) == pd.DataFrame, 'Affinity matrix must be in Pandas DataFrame format.'

        self.clusters = clusters
        self.affinity = affinity
        self.info = pd.read_csv(
            info_path,
            sep = '\t',
            index_col = 0
        )
        self.info.index = [x.split('.')[0] for x in self.info.index]
        self.amr_info = pd.read_csv(
            amr_path,
            sep = '\t',
            index_col = 0
        )
        self.amr_info.index = [x.split('.')[0] for x in self.amr_info.index]
        self.amr_info.index.name = 'Contig id'
        self.X = X
        self.cluster_by_plasmid = self.__build_cluster_by_plasmid(self.clusters)
        self.fontsizes = {}

    def draw_graph(
        self,
        path: str,
        size = ('700px', '700px'),
        background_color = '#ffffff',
        color_cycle = None,
        edge_threshold: float = 0.0,
        show_sizes = False,
        min_size = 0
    ):
        self.min_size = min_size
        self.network = self.__create_network(
            height = size[1],
            width = size[0],
            bgcolor = background_color
        )
        with plt.style.context('fivethirtyeight'):
            self.__clusters_to_nodes(
                self.cluster_by_plasmid, 
                len(self.clusters), 
                color_cycle, 
                show_sizes=show_sizes
            )
            self.__affinity_to_edges(self.affinity, edge_threshold)

        self.network.toggle_physics(False)
        self.network.show_buttons(['physics', 'nodes'])
        self.network.show(path)
        
        return self.network

    
    def __build_cluster_by_plasmid(
        self,
        clusters: pd.DataFrame
    ) -> pd.DataFrame:

        for n, cluster in enumerate(clusters.index):

            it_df = pd.DataFrame(
                    [n]*len(clusters.loc[cluster]['Plasmids']),
                    index = clusters.loc[cluster]['Plasmids'],
                    columns = ['Cluster']
                )

            if n == 0:
                cluster_by_plasmid = it_df
            else:
                cluster_by_plasmid = pd.concat(
                    [
                        cluster_by_plasmid,
                        it_df
                    ]
                )

        return cluster_by_plasmid


    def __create_network(
        self,
        height,
        width,
        bgcolor = '#ffffff'
    ) -> Network:
        return Network(
            height = height,
            width = width,
            bgcolor = bgcolor
        )

    def __clusters_to_nodes(
        self,
        clusters: pd.DataFrame,
        n_clusters: int,
        color_cycle = None,
        show_sizes = False
    ):
        if color_cycle is None:
            self.color_cycle = np.vstack(
                [
                    plt.cm.rainbow(np.linspace(0, 1, 20)), 
                    plt.cm.rainbow(np.linspace(0, 1, n_clusters-20))
                ]
            )
        else:
            self.color_cycle = color_cycle

        nodes = []
        sizes = []  # Number of genes
        colors = [] # By cluster
        labels = [] # Species, plasmid name, AMR genes
        self.__amr_genes = self.amr_info.groupby(
            ['Contig id', 'Gene symbol']
        ).size()
        cluster_sizes = clusters.groupby('Cluster').size()

        for plasmid in self.affinity.index:
            plasmid_cluster = clusters.loc[plasmid]['Cluster']
            if cluster_sizes[plasmid_cluster] >= self.min_size:
                try:
                    plasmid_info = self.info.loc[plasmid]
                    nodes.append(plasmid)
                    colors.append(
                        'rgb(' + ','.join(
                            list(
                                (
                                    self.color_cycle[
                                        clusters.loc[plasmid]['Cluster']
                                    ] * [255, 255, 255, 1]
                                ).astype(int).astype(str)
                            )
                        ) + ')'
                    )
                    if self.X is not None:
                        sizes.append(
                            len(
                                self.X.loc[plasmid][self.X.columns.to_list()[0]]
                            )
                        )
                    try:
                        amr_hits = self.__amr_genes.loc[plasmid].index.to_list()
                    except KeyError:
                        amr_hits = []
                    if self.X is not None:
                        labels.append(
                            ' '.join(
                                [
                                    'Organism:',
                                    plasmid_info['Organism'],
                                    '\n',
                                    'Plasmid:',
                                    plasmid_info['Plasmid'],
                                    '\n',
                                    'AMR Genes:',
                                ] + amr_hits + [
                                    '\n',
                                    'Total Number of Genes',
                                    str(len(
                                        self.X.loc[plasmid][self.X.columns.to_list()[0]]
                                    )),
                                    '\n',
                                    'Cluster: ',
                                    str(clusters.loc[plasmid]['Cluster'])
                                ]
                            )
                        )
                    else:
                        labels.append(
                            ' '.join(
                                [
                                    'Organism:',
                                    plasmid_info['Organism'],
                                    '\n',
                                    'Plasmid:',
                                    plasmid_info['Plasmid'],
                                    '\n',
                                    'AMR Genes:',
                                ] + amr_hits + [
                                    '\n',
                                    'Cluster: ',
                                    str(clusters.loc[plasmid]['Cluster'])
                                ]
                            )
                        )
                except KeyError:
                    pass
        
        if show_sizes:
            self.network.add_nodes(
                nodes = nodes,
                size = sizes,
                color = colors,
                title = labels,
                label = nodes
            )
        else:
            self.network.add_nodes(
                nodes = nodes,
                color = colors,
                title = labels,
                label = nodes
            )

    def __affinity_to_edges(
        self,
        affinity: pd.DataFrame,
        threshold: float = 0.0
    ):
        for node_a in self.network.get_nodes():
            for node_b in self.network.get_nodes():
                if affinity.loc[node_a][node_b] > threshold and node_a != node_b:
                    self.network.add_edge(
                        node_a,
                        node_b, 
                        value = affinity.loc[node_a][node_b]
                    )

    def __plot_species_by_cluster(
        self,
        clusters,
        info,
        ax,
        barw = .70
    ):
        n_clusters = len(clusters)
        plot_data = {
            k: np.zeros(n_clusters)
            for k in np.unique(info['Organism'])
        }
        for n, cluster in enumerate(clusters.index):
            species, counts = np.unique(
                info.loc[
                    clusters.loc[cluster]['Plasmids']
                ]['Organism'],
                return_counts = True
            )
            for species_k, n_plasmids in zip(
                species,
                counts
            ):
                plot_data[species_k][n] = n_plasmids
        
        accumulate = np.zeros(n_clusters)
        for species, counts in plot_data.items():
            if np.any(counts):
                ax.bar(
                    np.arange(n_clusters),
                    counts,
                    bottom = accumulate,
                    width = barw,
                    label = '\n'.join(species.split())
                )
                accumulate += counts

        ax.legend(
            loc = 'upper right',
            fontsize = self.fontsizes['stats'],
            frameon = True,
            facecolor = ax.get_facecolor(),
            framealpha = 1,
            edgecolor = ax.get_facecolor()
        )
        ax.set_xticks(
            np.arange(n_clusters)
        )
        ax.tick_params(
            axis = 'y',
            labelsize = self.fontsizes['stats']
        )
        ax.set_xticklabels(
            [
                f'Cluster {n+1}'
                for n in range(n_clusters)
            ],
            fontsize = self.fontsizes['stats'],
            ha = 'center'
        )
        ax.set_ylabel(
            'Number of Plasmids',
            fontsize = self.fontsizes['stats']
        )
        ax.grid(
            False,
            axis = 'x'
        )

    def __plot_amr_per_cluster(
        self,
        amr,
        clusters,
        ax,
        bw = .5
    ):
        amr_classes = np.unique(amr['Class'])
        n_classes = len(amr_classes)
        n_clusters = len(clusters)
        plot_data = pd.DataFrame(
            np.zeros((n_classes+1, n_clusters)),
            index = list(amr_classes) + ['No Resistance'],
            columns = clusters.index
        )

        for n in clusters.index:
            for m, amr_class in enumerate(amr_classes):
                cluster_plasmid_ids = clusters.loc[n]['Plasmids']
                common_ids = set(
                    amr.index
                ).intersection(
                    cluster_plasmid_ids
                )
                plot_data.loc['No Resistance'][n] = len(
                    set(
                        cluster_plasmid_ids
                    ).difference(
                        amr.index
                    )
                )
             
                try:
                    plot_data.loc[amr_class][n] = len(
                        amr.loc[
                            common_ids
                        ].groupby(
                            ['Class', 'Contig id']
                        ).size().loc[amr_class].index
                    )
                except KeyError as e:
                    pass

        plot_data.drop(
            plot_data.index[
                ~np.any(plot_data, axis=1)
            ],
            axis = 'index',
            inplace = True
        )
        plot_data = plot_data.loc[
            plot_data.sum(axis = 1).sort_values(ascending=False).index
        ]
        for n, m in enumerate(plot_data.columns):
            
            if n_clusters > 6:
                color_cycle = plt.cm.rainbow(np.linspace(0, 1, n_clusters-6))
            if n < 6:
                ax.bar(
                    np.arange(plot_data.shape[0]),
                    plot_data[m],
                    width = bw,
                    label = f'Cluster {n+1}',
                    bottom = np.sum(plot_data.iloc[:,:n], axis=1)
                )
            else:
                ax.bar(
                    np.arange(plot_data.shape[0]),
                    plot_data[m],
                    width = bw,
                    label = f'Cluster {n+1}',
                    bottom = np.sum(plot_data.iloc[:,:n], axis=1),
                    color = color_cycle[n-6]
                )
                    
        ax.set_xticks(
            np.arange(plot_data.shape[0])
        )
        ax.set_xticklabels(
            [x.lower().capitalize() for x in plot_data.index],
            ha = 'right',
            rotation = 45,
            fontsize = 10
        )
        ax.set_ylabel(
            'Number of Plasmids with Resistance',
            fontsize = self.fontsizes['stats']
        )
        ax.tick_params(
            axis = 'y',
            labelsize = self.fontsizes['stats']
        )
        ax.set_xlabel(
            'Antimicrobial Class',
            fontsize = self.fontsizes['stats']
        )
        ax.legend(
            loc = 'upper right',
            fontsize = self.fontsizes['stats'],
            ncol = 2,
            frameon = True,
            facecolor = ax.get_facecolor(),
            framealpha = 1,
            edgecolor = ax.get_facecolor()
        )
        ax.grid(
            False,
            axis = 'x'
        )

    def plot_stats(
        self,
        figsize = (12,7),
        constrained_layout = True,
        species_by_cluster_bw = .70,
        amr_per_cluster_bw = .70,
        fontsize = 10,
        heatmap = False
    ):
        self.fontsizes['stats'] = fontsize

        self.stats_figure, axs = plt.subplots(
            1,2,
            figsize = figsize,
            constrained_layout = constrained_layout
        )
        if heatmap:
            self.__plot_species_by_cluster_heatmap(
                self.clusters,
                self.info,
                axs[0]
            )
            self.__plot_amr_per_cluster_heatmap(
                self.amr_info,
                self.clusters,
                axs[1]
            )
        else:
            self.__plot_species_by_cluster(
                self.clusters,
                self.info,
                axs[0],
                species_by_cluster_bw
            )

            self.__plot_amr_per_cluster(
                self.amr_info,
                self.clusters,
                axs[1],
                amr_per_cluster_bw
            )
        plt.plot()

    def plot_proteins(
        self,
        proteins_by_plasmid,
        protein_names: pd.Series,
        fontsize = 10,
        figsize = (12, 7),
        clusters_to_include = [0, 1],
        threshold = 0
    ):
        clusters = self.clusters
        n_clusters = len(clusters_to_include)
        self.fontsizes['proteins'] = fontsize

        fig, axs = plt.subplots(
            int(np.ceil(n_clusters/4)),
            np.minimum(4, n_clusters),
            figsize = figsize,
            constrained_layout = True
        )
        for cluster, ax in zip(
            clusters_to_include,
            axs.flatten()
        ):
            n_plasmids = len(
                clusters.loc[cluster]['Plasmids']
            )
            all_proteins = proteins_by_plasmid.loc[
                clusters.loc[cluster]['Plasmids']
            ]['Proteins']
            for n, proteins in enumerate(all_proteins):
                all_proteins.iloc[n] = np.unique(
                    proteins
                )
            all_proteins = np.hstack([x for x in all_proteins])
            unique_proteins, protein_counts = np.unique(
                all_proteins,
                return_counts = True
            )
            proteins_in_other_clusters = proteins_by_plasmid.loc[
                clusters.loc[
                    [
                        x for x in clusters.index
                        if x != cluster
                    ]
                ]['Plasmids'].sum()
            ].sum(axis = 0)['Proteins']
            unique_proteins_outside, protein_counts_outside = np.unique(
                proteins_in_other_clusters,
                return_counts = True
            )
            protein_df = pd.concat(
                [
                    pd.Series(
                        protein_counts,
                        index = unique_proteins
                    ),
                    pd.Series(
                        protein_counts_outside,
                        index = unique_proteins_outside
                    )
                ],
                axis = 1
            )
            protein_df[np.isnan(protein_df)] = 0
            protein_df.columns = ['In', 'Out']
            protein_df = protein_df.astype(int)

            proteins_in = protein_df['In']/n_plasmids
            proteins_in = proteins_in[proteins_in > threshold].sort_values()
            ax.scatter(
                proteins_in,
                np.arange(len(proteins_in)),
                label = 'Inside the Cluster'
            )
            ax.scatter(
                protein_df['Out'].loc[proteins_in.index]/(len(proteins_by_plasmid) - n_plasmids),
                np.arange(len(proteins_in)),
                label = 'Outside the Cluster'
            )
            ax.set_xlabel(
                'Number of Occurences (% Plasmids)',
                fontsize = fontsize
            )

            ax.set_yticks(
                np.arange(len(proteins_in))
            )
            ax.set_yticklabels(
                protein_names.loc[proteins_in.index],
                fontsize = fontsize,
                ha = 'right'
            )

            ax.set_title(
                f'Cluster {cluster+1}',
                fontsize = fontsize*1.3
            )

            ax.set_xticks(
                np.arange(
                    0,
                    1.1,
                    .1
                )
            )
            ax.set_xticklabels(
                (np.arange(
                    0,
                    1.1,
                    .1
                )*100).astype(int),
                fontsize = fontsize,
                ha = 'center'
            )

            ax.legend(
                facecolor = ax.get_facecolor(),
                fontsize = fontsize,
                edgecolor = ax.get_facecolor()
            )

        plt.plot()

    def __plot_species_by_cluster_heatmap(
        self,
        clusters,
        info,
        ax,
    ):
        selected_info = info.loc[self.affinity.index]['Organism']
        counts = pd.concat(
            [
                self.cluster_by_plasmid,
                selected_info
            ],
            join='outer',
            axis = 1
        ).groupby(['Organism', 'Cluster']).size()
        n_clusters = len(clusters)
        species = np.unique(
            counts.index.get_level_values('Organism')
        )
        X_img = np.zeros((len(species), n_clusters))
        for n, specie in enumerate(species):
            for cluster in counts.loc[specie].index.get_level_values('Cluster'):
                X_img[n, cluster] = counts.loc[specie][cluster]

        X_img = X_img.transpose()
        alphas = np.zeros_like(X_img)
        alphas[X_img > 0] = 1
        ax.imshow(
            X_img,
            cmap = 'inferno',
            alpha = alphas,
            aspect = 'auto'
        )
        for i in range(X_img.shape[0]):
            for j in range(X_img.shape[1]):
                if X_img[i, j] > 0.5*np.max(X_img):
                    color = 'black'
                else:
                    color = 'white'
                if X_img[i, j] > 0:
                    ax.text(
                        j,
                        i,
                        int(X_img[i, j]),
                        fontsize = 9,
                        ha = 'center',
                        va = 'center',
                        color = color,
                        weight = 'bold'
                    )

        ax.set_yticks(
            np.arange(0,n_clusters,2)
        )
        ax.tick_params(
            axis = 'x',
            labelsize = self.fontsizes['stats']
        )
        ax.set_yticklabels(
            [
                f'{n+1}'
                for n in range(0,n_clusters,2)
            ],
            fontsize = self.fontsizes['stats']
        )
        ax.set_xticks(
            np.arange(len(species))
        )
        ax.set_xticklabels(
            ['\n'.join(x.split()) for x in species],
            fontsize = self.fontsizes['stats'],
            ha = 'right',
            rotation = 45
        )
        ax.set_ylabel(
            'Cluster',
            fontsize = self.fontsizes['stats']
        )
        ax.set_yticks(
            np.arange(0.5, n_clusters+.5, 1),
            minor = True,
        )
        ax.set_xticks(
            np.arange(0.5, len(species)+.5, 1),
            minor = True,
        )
        ax.grid(
            False
        )
        ax.grid(
            True,
            axis = 'y',
            which = 'minor',
            c = ax.get_facecolor()
        )
        ax.grid(
            True,
            axis = 'x',
            which = 'minor',
            c = ax.get_facecolor()
        )

        plt.colorbar(
            mappable = ScalarMappable(
                norm = Normalize(
                    vmin=0, vmax=np.max(
                        counts.values.flatten()
                    )
                ), cmap = 'inferno'
            ),
            ax = ax,
            label = 'Number of Plasmids'
        )

    def __plot_amr_per_cluster_heatmap(
        self,
        amr,
        clusters,
        ax,
    ):
        amr_classes = np.unique(amr['Class'])
        n_classes = len(amr_classes)
        n_clusters = len(clusters)
        plot_data = pd.DataFrame(
            np.zeros((n_classes+1, n_clusters)),
            index = list(amr_classes) + ['No Resistance'],
            columns = clusters.index
        )

        for n in clusters.index:
            for m, amr_class in enumerate(amr_classes):
                cluster_plasmid_ids = clusters.loc[n]['Plasmids']
                common_ids = set(
                    amr.index
                ).intersection(
                    cluster_plasmid_ids
                )
                plot_data.loc['No Resistance'][n] = len(
                    set(
                        cluster_plasmid_ids
                    ).difference(
                        amr.index
                    )
                )
             
                try:
                    plot_data.loc[amr_class][n] = len(
                        amr.loc[
                            common_ids
                        ].groupby(
                            ['Class', 'Contig id']
                        ).size().loc[amr_class].index
                    )
                except KeyError as e:
                    pass

        plot_data.drop(
            plot_data.index[
                ~np.any(plot_data, axis=1)
            ],
            axis = 'index',
            inplace = True
        )
        plot_data = plot_data.loc[
            plot_data.sum(axis = 1).sort_values(ascending=False).index
        ]

        alphas = np.zeros_like(plot_data.values.transpose())
        alphas[plot_data.values.transpose() > 0] = 1
        ax.imshow(
            plot_data.values.transpose(),
            cmap = 'inferno',
            alpha = alphas,
            aspect = 'auto'
        )            
                    
        ax.set_xticks(
            np.arange(plot_data.shape[0])
        )
        ax.set_xticks(
            np.arange(.5, plot_data.shape[0]+.5),
            minor=True
        )
        ax.set_xticklabels(
            [x.lower().capitalize() for x in plot_data.index],
            ha = 'right',
            rotation = 45,
            fontsize = 10
        )
        ax.set_ylabel(
            'Cluster',
            fontsize = self.fontsizes['stats']
        )
        ax.tick_params(
            axis = 'y',
            labelsize = self.fontsizes['stats']
        )
        ax.set_yticks(
            np.arange(
                0,
                len(
                    plot_data.columns
                ),
                2
            )
        )
        ax.set_yticks(
            np.arange(.5, len(plot_data.columns)+.5),
            minor = True
        )
        ax.set_yticklabels(
            [x + 1 for x in list(plot_data.columns)[::2]]
        )
        ax.grid(
            False
        )
        ax.grid(
            True,
            axis = 'y',
            which = 'minor',
            c = ax.get_facecolor()
        )
        ax.grid(
            True,
            axis = 'x',
            which = 'minor',
            c = ax.get_facecolor()
        )
        plt.colorbar(
            mappable = ScalarMappable(
                norm = Normalize(
                    vmin = 0,
                    vmax = np.max(plot_data.values.flatten())
                ),
                cmap = 'inferno'
            ),
            ax = ax,
            label = 'Number of Plasmids'
        )
        self.amr_data = plot_data
# %%

class SubClusterImages:

    def __init__(
        self,
        phylo: SimpleDendogram,
        genomes: Iterable,
        path_to_info: str,
        path_to_amr: str,
        path_to_linkage: str,
        figsize: Tuple = (10, 7),
        fontsize: float = 10.,
        distance_threshold: float = 70.
    ):
        self.figsize = figsize
        self.fontsize = fontsize
        self.__ccycle = self.__get_color_cycle()
        self.path_to_info = path_to_info
        self.path_to_amr = path_to_amr
        self.linkage, self.labels = joblib.load(path_to_linkage)
        self.phylo = phylo
        self.distance_threshold = distance_threshold
        self.CMAP = 'Set1'
        self.WITH_NUMBERS = True
        self.LW = 1
        self.genomes = genomes

    def __get_color_cycle(
        self,
        style = 'ggplot'
    ):
        with plt.style.context(style):
            return plt.rcParams['axes.prop_cycle'].by_key()['color']

    def __set_axis_labels(
        self,
        ax,
        x: str,
        y: str
    ):

        ax.set_xlabel(
            x,
            fontsize = self.fontsize
        )
        ax.set_ylabel(
            y,
            fontsize = self.fontsize
        )

    def __boxplot(
        self,
        ax,
        x
    ):
        return ax.boxplot(
            x,
            patch_artist = True,
            boxprops={
                'fill': True,
                'facecolor': self.__ccycle[0],
                'edgecolor': self.__ccycle[0]
                },
            medianprops={
                'color': self.__ccycle[4],
                'linewidth': 3
                },
            capprops={'color': self.__ccycle[0]},
            flierprops={'markeredgecolor': self.__ccycle[0]},
            whiskerprops={'color': self.__ccycle[0]}
        )

    def plot_sizes(
        self,
        clusters
    ):

        sizes = clusters.groupby('Cluster').size().sort_index()

        with plt.style.context('ggplot'):

            fig, ax = plt.subplots(
                1, 2,
                figsize = (10,10),
                constrained_layout = True,
                gridspec_kw={'width_ratios': [5, 1]}
            )
            ax[0].bar(
                sizes.index,
                sizes.values
            )
            self.__set_axis_labels(
                ax[0],
                'Cluster',
                'Number of Plasmids'
            )
            ax[0].set_yticks(
                np.arange(0, max(sizes.values)+1, dtype = int)
            )
            if len(sizes) < 10:
                step = 1
            else:
                step = 2
            ax[0].set_xticks(
                np.arange(0, len(sizes), step)
            )
            ax[0].set_xticklabels(
                np.arange(1, len(sizes)+1, step)
            )
            ax[0].set_xlim(-0.5, len(sizes))

            self.__boxplot(ax[1], sizes.values)
            self.__set_axis_labels(
                ax[1],
                None,
                'Number of Plasmids'
            )
            ax[1].set_yticks(
                np.arange(0, max(sizes.values)+1, dtype = int)
            )

        plt.show()

    def __join_info(
        self,
        clusters
    ):

        info = pd.read_csv(
            self.path_to_info,
            sep = '\t',
            index_col = 0
        )
        info.index = [x.split('.')[0] for x in info.index]
        info = info.loc[clusters.index]

        self.clusters = pd.concat(
            [
                clusters,
                info
            ],
            join = 'inner',
            axis = 'columns'
        )

    def plot_dendrogram_info(
        self,
        clusters,
        linkage,
        labels
    ):
        with plt.style.context('ggplot'):
            fig, ax = plt.subplots(
                    1,
                    3,
                    figsize = self.figsize,
                    constrained_layout = True,
                    gridspec_kw={'width_ratios': [4, .5, 1]}
                )
            fig.set_facecolor('white')
            for x in range(1, len(ax)):
                ax[x].set_facecolor('white')
            self.dendrogram = self.plot_dendrogram(
                ax[0],
                clusters,
                linkage,
                labels
            )
            self.label_order = [x.split()[0] for x in self.dendrogram['ivl']]
            self.plot_species(
                ax[1],
                labels
            )
            self.plot_amr_classes(ax[2])
        plt.show()

    def plot_species(
        self,
        ax,
        labels
    ):
        species = self.clusters['Organism'].apply(
            lambda x: ' '.join(
                [
                    x.split()[0][0]+'.',
                    x.split()[1]
                ]
            )
        )
        species = species[self.label_order]
        unique_species = np.unique(species)
        n_species = len(unique_species)
        n_plasmids = len(species)

        species_matrix = np.array([
            [
                plasmid_species == species_name
                for species_name in unique_species   
            ]
            for plasmid_species in species
        ])[::-1] 
        
        self.colored_species_matrix = species_matrix * np.arange(1, n_species+1)
        ax.imshow(
            self.colored_species_matrix,
            cmap = self.CMAP,
            alpha = species_matrix.astype(float),
            aspect = 'auto'
        )
        ax.set_xticks(
            np.arange(n_species)
        )
        ax.set_xticklabels(
            unique_species,
            rotation = 90,
            ha = 'right'
        )
        ax.set_yticks(
            np.arange(n_plasmids) + .5
        )
        ax.set_xlim(-.5, n_species)
        
        ax.grid(
            True,
            axis = 'y',
            lw = self.LW,
            color = ax.get_facecolor()
        )
        ax.grid(
            False,
            axis = 'x'
        )
        
        ax.set_yticklabels([])
        ax.set_xlim(-.5,n_species-.5)
        ax.tick_params(axis="y", left=False)

    def __load_amr(
        self
    ):
        all_amr = pd.read_csv(
            self.path_to_amr,
            sep = '\t',
            index_col = 0
        )
        all_amr.index = [x.split('.')[0] for x in all_amr.index]

        unique_classes = np.unique(all_amr['Class'])
        n_classes = len(unique_classes)
        n_plasmids = len(self.label_order)

        self.amr = pd.DataFrame(
            np.zeros(
                (n_plasmids, n_classes),
                dtype = int
            ),
            index = self.label_order,
            columns = unique_classes
        )
        self.amr_counts = pd.DataFrame(
            np.zeros(
                (n_plasmids, n_classes),
                dtype = int
            ),
            index = self.label_order,
            columns = unique_classes
        )

        for plasmid in self.label_order:
            if plasmid in all_amr.index.to_list():
                plasmid_amr = all_amr.loc[plasmid]
                plasmid_classes, plasmid_class_counts = np.unique(
                    plasmid_amr['Class'],
                    return_counts = True
                )
                for class_, count in zip(
                    plasmid_classes,
                    plasmid_class_counts
                ):
                    self.amr.loc[plasmid][class_] = 1
                    self.amr_counts.loc[plasmid][class_] = count

        all_zero_columns = self.amr.columns[~np.any(self.amr, axis=0)]
        self.amr.drop(
            all_zero_columns,
            inplace=True,
            axis = 'columns'
        )
        self.amr_counts.drop(
            all_zero_columns,
            inplace = True,
            axis = 'columns'
        )

    def plot_amr_classes(
        self,
        ax
    ):
        self.__load_amr()

        n_plasmids, n_classes = self.amr.shape
        amr_vals = self.amr.values[::-1]
        counts_vals = self.amr_counts.values[::-1]
        ax.imshow(
            (amr_vals*(np.arange(n_classes)+1)).astype(int),
            cmap = 'gist_rainbow',
            alpha = amr_vals.astype(float),
            aspect = 'auto'
        )
        ax.set_xticks(
            np.arange(n_classes)
        )
        ax.set_xticklabels(
            [x.lower().capitalize() for x in self.amr.columns],
            rotation = 90,
            ha = 'right'
        )
        ax.set_yticks(
            np.arange(n_plasmids)+.4
        )
        ax.set_xticks(
            np.arange(n_classes)+.45,
            minor = True
        )
        ax.grid(False, axis = 'x', which = 'major')
        ax.grid(
            True,
            axis = 'y',
            lw = self.LW,
            color = ax.get_facecolor()
        )
        ax.grid(
            True,
            which = 'minor',
            axis = 'x',
            lw = self.LW,
            color = ax.get_facecolor()
        )
        
        ax.set_yticklabels([])
        ax.set_xlim(-.5,n_classes-.5)
        ax.tick_params(axis="y", left=False)

        if self.WITH_NUMBERS:
            for n in range(n_classes):
                for m in range(n_plasmids):

                    if counts_vals[m, n] > 0:

                        ax.text(
                            n, m,
                            counts_vals[m, n],
                            ha = 'center',
                            va = 'center',
                            fontweight = 'bold',
                            fontsize = self.figsize[1]*.4 - 1*(n_plasmids>40),
                            color = ax.get_facecolor()
                        )

        ax.tick_params(axis="x", bottom=False, which = 'minor')
        for x, color_ in zip(
            range(n_classes),
            matplotlib.cm.gist_rainbow(np.linspace(0,1,n_classes+1))[1:]
        ):
            HEIGHT = .2
            MARGIN = .015
            ax.add_patch(
                Rectangle(
                    (x-.5+MARGIN, n_plasmids-HEIGHT-0.5),
                    width = 1-2*MARGIN,
                    height = HEIGHT,
                    facecolor = color_,
                    zorder = -10
                )
            )

    def plot_dendrogram(
        self,
        ax,
        clusters,
        linkage,
        labels,
        info_field: str = 'Organism'
    ):

        with plt.style.context('ggplot'):
            
            r = dendrogram(
                linkage,
                labels = labels,
                distance_sort = 'ascending',
                leaf_label_func = lambda x: ' '.join(
                    [
                        labels[x],
                        f'({clusters.iloc[x]["Cluster"]})'
                    ]
                ),
                leaf_rotation = 0,
                ax = ax,
                color_threshold = self.distance_threshold,
                orientation = 'left',
                above_threshold_color = 'black'
            )

            ax.set_ylabel(f'{info_field} (Cluster)')
            ax.set_xlabel('Syntenic Distance')
            ax.grid(False)

            return r

    def plot_synteny_jaccard_ratio(
        self
    ):
        self.affinity = self.phylo.get_affinity_as_frame()

        self.jaccard = deepcopy(self.affinity)
        for i in self.affinity:
            for j in self.affinity:
                self.jaccard.loc[i][j] = jaccard(
                    self.genomes[i].get_gene_set(), 
                    self.genomes[j].get_gene_set()
                )

        with plt.style.context('ggplot'):
            fig, axs = plt.subplots(
                1,1,
                constrained_layout = True,
                figsize = (10,10)
            )
            ax_hist = axs
            h, _, _, _ = ax_hist.hist2d(
                self.affinity.values.flatten(),
                self.jaccard.values.flatten(),
                density = False,
                cmin = 1,
                cmap = 'inferno',
                bins = [25, 25]
            )
            ax_hist.set_xlabel('Syntenic Distance')
            ax_hist.set_ylabel('Jaccard Score')

            plt.colorbar(
                mappable = ScalarMappable(
                    norm = Normalize(
                        vmin = 0,
                        vmax = np.max(h)
                    ),
                    cmap = 'inferno'
                ),
                ax = ax_hist,
                label = 'Number of Plasmid Pairs'
            )

        plt.show()

    def plot_size_distance_correlation(
        self
    ):
        plasmid_sizes = {
            k: len(v.get_gene_set())
            for k, v in self.genomes.items()
        }

        affinity = self.phylo.get_affinity_as_frame()
        distances_, sizes = [], []
        for n, i in enumerate(affinity):
            for j in affinity.columns.to_list()[n+1:]:
                distances_.append(affinity[i][j])
                sizes.append(
                    plasmid_sizes[i]+plasmid_sizes[j]
                )

        distances_ = np.array(distances_)
        sizes= np.array(sizes)

        with plt.style.context('ggplot'):

            f, ax = plt.subplots(
                1, 1,
                constrained_layout = True,
                figsize = (10,10)
            )

            ax.scatter(
                sizes,
                distances_,
                alpha = .7
            )
            ax.set_ylabel('Syntenic Distance', fontsize=15)
            ax.set_xlabel('Total Number of Genes in Both Plasmids', fontsize=15)
            correlation, pval = pearsonr(
                sizes,
                distances_,
                alternative='greater'
            )
            ax.text(
                ax.get_xlim()[1]*.05,
                max(ax.get_ylim())*.9,
                f'One-sided Pearson\'s $r$: {correlation:.5}\n$p$-value: {pval:.5}',
                fontsize = 15,
                fontweight = 'bold'
            )
            plt.show()

    def plot_jaccard_distance_correlation(
        self
    ):
        with plt.style.context('ggplot'):

            f, ax = plt.subplots(
                1, 1,
                constrained_layout = True,
                figsize = (10,10)
            )

            ax.scatter(
                self.jaccard.values.flatten(),
                self.affinity.values.flatten(),
                alpha = .7
            )
            ax.set_ylabel('Syntenic Distance', fontsize=15)
            ax.set_xlabel('Jaccard Index', fontsize=15)
            correlation, pval = pearsonr(
                self.jaccard.values.flatten(),
                self.affinity.values.flatten()
            )
            ax.text(
                ax.get_xlim()[0]+(ax.get_xlim()[1]-ax.get_xlim()[0])*.05,
                max(ax.get_ylim())*.9,
                f'Two-sided Pearson\'s $r$: {correlation:.5}\n$p$-value: {pval:.5}',
                fontsize = 15,
                fontweight = 'bold'
            )
            plt.show()

    def show(
        self,
        clusters,
        *,
        sizes: bool = True,
        species: bool = True,
        ratio: bool = True,
        jaccard_distance_correlation: bool = True,
        size_correlation: bool = True,

    ):
        self.__join_info(clusters)

        if sizes:
            self.plot_sizes(clusters)
        if species:
            self.plot_dendrogram_info(self.clusters, self.linkage, self.labels)
        if ratio:
            self.plot_synteny_jaccard_ratio()
        if jaccard_distance_correlation:
            self.plot_jaccard_distance_correlation()
        if size_correlation:
            self.plot_size_distance_correlation()
        

# %%

class PlasmidNetwork:

    def __init__(
        self,
        path_to_info: str,
        figsize: Tuple[str, str] = ('700px', '700px'),
        edge_threshold: float = 0.0,
        hide_threshold: float = 0.9
    ):
        '''
        Parameters
        ----------

        path_to_info: str
            Path to plasmid_info.csv.

        figsize: Tuple of str.
            Canvas size.

        edge_threshold: float in [0, 1]. Default is 0
            Edges with weights below this value (fraction of the maximum edge weight as given
            by the affinity matrix passed to show()) will not be included in the graph and will
            not take part in the node layout construction.

        hide_threshold: float in [0, 1]. Default is 0.9
            Edges with weights less than this value will not be shown on the graph, but still 
            influence the layout construction. Useful to unclutter the plot if a large number of
            edges are present.
        '''

        self.info = pd.read_csv(path_to_info, sep='\t', index_col=0)
        self.info.index = [x.split('.')[0] for x in self.info.index]
        self.figsize = figsize
        self.NODE_SIZE = 20
        self.__threshold = edge_threshold
        self.__hide_t = hide_threshold

    def __scale_distance(
        self,
        distance: pd.DataFrame
    ):

        min_nonzero = np.min(distance.values[distance.values > 0])

        scaled_dist = deepcopy(distance)
        #scaled_dist.values[scaled_dist.values>0] = scaled_dist.values[scaled_dist.values>0] - min_nonzero + 1e-2
        max_val = np.max(scaled_dist.values)

        return scaled_dist/max_val

    def __distance_to_affinity(
        self,
        distance
    ):

        return 1.0-distance

    def __build(
        self,
        clusters: pd.Series,
        affinity: pd.DataFrame
    ):

        cluster_names = np.unique(clusters)
        n_clusters = len(cluster_names)
        colors = plt.cm.rainbow(np.linspace(0, 1, n_clusters))

        self.network = Network(
            self.figsize[1],
            self.figsize[0],
            bgcolor='#FFFFFF'
        )

        for cluster, color in zip(cluster_names, colors):
            plasmids_in_cluster = clusters[clusters==cluster].index.to_list()
            for plasmid in plasmids_in_cluster:
                try:
                    self.network.add_node(
                        plasmid,
                        label = plasmid,
                        title = '\n'.join(
                            [
                                plasmid,
                                'Species: ' + self.info.loc[plasmid]['Organism'],
                                'Strain: ' + str(self.info.loc[plasmid]['Strain']),
                                'Plasmid: ' + str(self.info.loc[plasmid]['Plasmid']),
                                'Cluster: ' + str(cluster)
                            ]
                        ),
                        color = self.__rbga_to_hex(color),
                        size = self.NODE_SIZE
                    )
                except KeyError:
                        self.network.add_node(
                        plasmid,
                        label = plasmid,
                        title = '\n'.join(
                            [
                                plasmid,
                                'Could not find this plasmid in the information database.',
                                'Cluster: ' + str(cluster)
                            ]
                        ),
                        color = self.__rbga_to_hex(color),
                        size = self.NODE_SIZE
                    )

        for n in range(len(affinity)):
            for m in range(n+1, len(affinity)):

                from_ = affinity.index.to_list()[n]
                to_ = affinity.index.to_list()[m]
                weight = affinity[from_][to_]

                if weight > self.__threshold:
                    if weight > self.__hide_t:
                        hide_ = False
                    else:
                        hide_ = True

                    self.network.add_edge(
                        from_,
                        to_,
                        value = weight,
                        hidden = hide_
                    )

    def __show(
        self,
        path
    ):
        if path is None:
            path_ = os.path.join(
                os.getcwd(),
                'temp_network_plot.html'
            )
        else:
            path_ = path

        self.network.show(
            path_
        )
        webbrowser.open(
            'file://///wsl.localhost/Ubuntu'+path_,
            new = 1
        )

        sleep(5)

        if path is None:
            os.remove(path_)

    def __set_options(
        self
    ):
        self.network.toggle_physics(False)
        self.network.force_atlas_2based()
        self.network.show_buttons(filter_ = ['nodes', 'physics'])

    def __rbga_to_hex(
        self,
        rgba_:  Iterable
    ) -> str:

        rgba = deepcopy(rgba_)
        # Convert to [0, 255]
        rgba *= 255
        rgba = rgba[:-1].astype(int)
        rgba = np.minimum(255, rgba)
        rgba = np.maximum(0, rgba)


        return "#{0:02x}{1:02x}{2:02x}".format(*rgba)

    def show(
        self,
        clusters: pd.Series,
        distance: pd.DataFrame,
        path: Union[str, None] = None
    ):
        '''
        Draws a pyVis Network graph of plasmids, given a distance matrix. The nodes represent
        plasmids, with edges weighted by the affinity matrix. These nodes will be coloured as
        specified in the clusters Pandas Series.

        The distances are rescaled into [0, 1] by dividing them by the maximum value and
        subtracting the minimum non-zero value to all non-zero values. Then, these are converted
        to an affinity matrix by subtracting the distance matrix to 1.

        Parameters
        ----------

        clusters: Pandas Series
            Series containing the cluster membership of each plasmid. Plasmid IDs should act as
            primary keys, with the cluster number or name being the Series value.

        distance: Pandas DataFrame
            Pandas DataFrame containing the distance matrix. Index and columns should contain 
            the same plasmid IDs, without the trailing .1 or .2. This DataFrame can be obtained
            from a SimpleDendogram object by calling get_affinity_as_frame().

        path: None or str. Default is None
            Path in which to save the graph HTML file. If None, the file is temporarily stored in
            the current working directory, and deleted after showing it in the web browser.
        '''

        distance = self.__scale_distance(distance)
        self.affinity = self.__distance_to_affinity(distance)
        self.__build(
            clusters,
            self.affinity
        )
        self.__set_options()
        self.__show(path)
# %%
