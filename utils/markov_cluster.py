'''
Implements Markov Clustering for Jaccard-based plasmid similarity networks.
'''

from argparse import ArgumentError, ArgumentTypeError
from turtle import distance

from sklearn.datasets import make_multilabel_classification
from utils.base_cluster import Clustering
from numpy import array
from markov_clustering import run_mcl, get_clusters, draw_graph, modularity
from utils.jaccard import jaccard
import pandas as pd
import numpy as np
from copy import deepcopy
import networkx as nx
from typing import Union, Iterable
import warnings
import os
import joblib
from sklearn.cluster import AgglomerativeClustering

class MarkovClustering(Clustering):

    def __init__(
        self,
        inflation: float,
        similarity_threshold: Union[float, str] = .5,
        return_type: str = 'numpy',
        reg_factor: float = .5,
        resolution: float = .8,
        score: str = 'core',
        core_threshold: float = 0.7
    ):
        '''
        Clusters a graph from an affinity matrix, using Markov Clustering.

        Arguments
        ----------

        inflation: float
            Inflation parameter for Markov Clustering. Higher values mean higher granularity.

        similarity_threshold: float | str
            Threshold for building the adjacency matrix. Edges with weights bellow or equal to
            this value will not be included in the adjacency matrix. If 'none', the affinity 
            matrix will be used as a starting point to the Markov Clustering, instead of the
            adjacency matrix.

        return_type: str
            If 'numpy', returns clusters as an array. If 'pandas', returns clustered communities
            as a Pandas DataFrame.

        reg_factor: float
            Regularization factor penalizing the number of communities. When score() is called,
            reg_factor * Number of communities / Number of edges is subtracted to the modularity.

        resolution: float
            Must be a positive number. If "score" is "modularity", then it represents the resolution 
            parameter for modularity. Numbers less than 1 favour a smaller number of clusters. 
            Numbers above 1 favour smaller clusters. If score is "core", then it represents the "eta"
            parameter, penalizing the number of clusters. Large numbers lead to a small number of clusters.
            Should be less than 1.

        score: str, "core" or "modularity". Default is "core"
            Default fitness function for hyperparameter selection. If "core", this function is
            set to $$-\\eta N_{clusters} + N_{clusters with zero core genes}$$. If "modularity",
            uses graph modularity.

        core_threshold: float. Default is 0.7.
            Fraction of plasmids in which a gene must be present to be considered a core gene.
        
        '''
        
        self.inflation = inflation
        self.sim_f = jaccard
        self.similarity_threshold = similarity_threshold
        self.reg_factor = reg_factor
        self.resolution = resolution
        self.score_ = score
        self.core_threshold = core_threshold
        super().__init__(return_type)
    
    def _Clustering__cluster(
        self, 
        affinity: array
    ) -> array:

        if self.similarity_threshold == 'none':
            self.adjacency = self.affinity
        else:
            self.adjacency = self.build_adjacency_matrix(
                as_sparse = False,
                threshold=self.similarity_threshold
            )
        self.__mcl_result = run_mcl(
            self.adjacency,
            inflation = self.inflation
        )
        self.clusters = get_clusters(self.__mcl_result)
        self.__np_clusters = deepcopy(self.clusters)
        return self.clusters

    def _Clustering__predictions_to_pandas(
        self,
        y
    ) -> pd.DataFrame:
        return pd.DataFrame(
            [  
                [[
                    self.index[plasmid_number]
                    for plasmid_number in cluster
                ]]
                for cluster in y
            ],
            index = np.arange(len(self.clusters)),
            columns = ['Plasmids']
        )

    def set_params(self, params):
        for k, v in params.items():
            if k == 'similarity_threshold':
                self.similarity_threshold = v
            elif k == 'inflation':
                self.inflation = v
            elif k == 'return_type':
                self.return_type = v
            else:
                raise ArgumentError(f'MarkovClustering has no attribute {k}.')

    def score(
        self,
        X
    ) -> float:
        if self.score_ == 'modularity':
            return nx.algorithms.community.modularity(
                nx.from_numpy_matrix(self.adjacency),
                self.__np_clusters,
                resolution = self.resolution
            ) - len(self.__np_clusters)/len(self.adjacency) * self.reg_factor
        elif self.score_ == 'core':
            n_core = []
            clusters = self._Clustering__predictions_to_pandas(self.clusters)
            for cluster in clusters.index:
                plasmids = clusters.loc[cluster]['Plasmids']
                proteins = self.X.loc[plasmids]['Proteins']
                unique_proteins = proteins.apply(lambda x: np.unique(x)).values
                unique_proteins = np.hstack(unique_proteins)
                unique_proteins, counts = np.unique(unique_proteins, return_counts=True)
                n_plasmids = len(plasmids)
                threshold = n_plasmids * self.core_threshold
                n_core.append(len(counts[counts>=threshold]))

            self.n_core = np.array(n_core)
            self.n_zero_core = sum(self.n_core==0)
            return -(len(self.clusters)*self.resolution + self.n_zero_core)
    
class JaccardClusters:

    def __init__(
        self,
        path_to_proteins: str,
        similarity: int = 9,
        k: int = 5,
        inflation: Union[float, Iterable] = 1.4,
        reg_factor: float = 1.0,
        aggregate: bool = True,
        min_cluster_size: int = 10,
        min_similarity: float = 0.02,
        drop_isolates = True,
        score: str = 'core',
        core_threshold: float = 0.7,
        resolution: float = 0.8,
        n_clusters: Union[int, None] = None,
        distance_threshold: float = 0.8,
        max_size: int = 80,
        linkage: str = 'complete'
    ):
        self.path = path_to_proteins
        self.similarity = similarity
        self.k = k
        self.inflation = inflation
        self.reg_factor = reg_factor
        self.__grid_search = (type(inflation) == list)
        self.aggregate = aggregate
        self.min_cluster_size = min_cluster_size
        self.min_similarity = min_similarity
        self.drop_isolates = drop_isolates
        self.score_ = score
        self.core_threshold = core_threshold
        self.resolution = resolution
        self.max_size = max_size
        self.distance_threshold = distance_threshold
        self.n_clusters = n_clusters
        self.linkage = linkage

        if self.__grid_search:
            self.clustering = MarkovClustering(
                inflation = 1.1,
                return_type = 'pandas',
                reg_factor = self.reg_factor,
                similarity_threshold = 'none',
                score = self.score_,
                core_threshold = self.core_threshold,
                resolution = self.resolution
            )
        else:
            self.clustering = MarkovClustering(
                inflation = self.inflation,
                return_type = 'pandas',
                reg_factor = self.reg_factor,
                similarity_threshold = 'none',
                score = self.score_,
                core_threshold = self.core_threshold,
                resolution = self.resolution
            )
        self.__get_proteins()

    def __get_proteins(self):
        self.proteins = joblib.load(
            os.path.join(
                self.path,
                f'Plasmids with Clustered Proteins_s{self.similarity}_k{self.k}.pkl'
            )
        )

    def cluster(
        self
    ):
        if self.__grid_search:
            self.clusters = self.clustering.grid_search(
                self.proteins,
                params = {
                    'inflation': self.inflation,
                }
            )
        else:
            self.clusters = self.clustering.fit_predict(
                self.proteins
            )

        self.affinity = pd.DataFrame(
            self.clustering.affinity,
            index = self.proteins.index,
            columns = self.proteins.index
        )

        if self.aggregate:
            aglomeration = Aglomeration(
                clusters = deepcopy(self.clusters),
                affinity = self.affinity,
                min_similarity = self.min_similarity,
                min_size = self.min_cluster_size
            )

            self.clusters = aglomeration.cluster().to_frame()
            self.clusters.columns = ['Plasmids']

            splitter = Splitter(
                self.clusters,
                self.affinity,
                distance_threshold = self.distance_threshold,
                linkage = self.linkage,
                n_clusters = self.n_clusters,
                max_size = self.max_size
            )
            self.clusters = splitter.split()

            if self.drop_isolates:
                self.clusters, membership = aglomeration.drop_isolates()
                self.affinity = aglomeration.affinity.loc[membership.index][membership.index]
        
        return self.clusters

    def to_tsv(
        self,
        path: str
    ):
        self.clusters.to_csv(
            os.path.join(
                path,
                f'Jaccard-Clusters-s{self.similarity}-k{self.k}.tsv'
            ),
            sep = '\t'
        )
        self.affinity.to_csv(
            os.path.join(
                path,
                f'Jaccard-Affinity-s{self.similarity}-k{self.k}.tsv'
            ),
            sep = '\t'
        )
        joblib.dump(
            self.clusters,
            os.path.join(
                path,
                f'Jaccard-Clusters-s{self.similarity}-k{self.k}.pkl'
            ),
            compress = 3
        )

class Aglomeration:

    def __init__(
        self,
        clusters: pd.DataFrame,
        affinity: pd.DataFrame,
        min_similarity: float = 0,
        min_size: int = 2
    ):
        self.membership = self.__build_membership(clusters)
        self.affinity = affinity
        self.min_similarity = min_similarity
        self.min_size = min_size

        self.cluster_counts = self.membership.groupby(['Cluster']).size()
        # Clusters to aggregate
        self.__to_fix = self.cluster_counts[self.cluster_counts < self.min_size].index.to_list()
        # Clusters without NNs with similarity > threshold
        self.__prepetual_isolates = []

    def __find_nn_cluster(
        self,
        affinities_: Union[pd.DataFrame, pd.Series]
    ):
        affinities = deepcopy(affinities_)
        if type(affinities) == pd.Series: 
            plasmid = [affinities.name]
        else: plasmid = affinities.index.to_list()

        # Make sure that it is not assigned its own cluster
        affinities[plasmid] = -1

        affinities = pd.concat(
            [affinities.loc[x] for x in affinities.index]
        )
        max_score = affinities.max()
        neighbor = affinities[affinities == max_score].index.to_list()[0]
        
        return self.membership.loc[neighbor]['Cluster'], max_score

    def __cluster_iteration(
        self
    ):
        for cluster in self.__to_fix:
            # Find all plasmids in that cluster
            plasmids = self.membership[self.membership['Cluster'] == cluster].index
            new_cluster, similarity = self.__find_nn_cluster(
                self.affinity.loc[plasmids]
            )
            if similarity > self.min_similarity:

                # Change membership of all plasmids in the old cluster
                # This way, if plasmid A is assigned to plasmid B cluster,
                # plasmid B is an isolate, and then plasmid B is assigned
                # to another cluster, A will follow B to the new assignment
                # and won't be turned into an isolate again
                self.membership.loc[plasmids] = new_cluster
            else:
                self.__prepetual_isolates.append(cluster)

    def cluster(
        self
    ):
        while len(self.__to_fix) > 0:

            self.__cluster_iteration()
            self.cluster_counts = self.membership.groupby(['Cluster']).size()
            # Clusters to aggregate
            self.__to_fix = self.cluster_counts[self.cluster_counts < self.min_size].index.to_list()
            self.__to_fix = [x for x in self.__to_fix if x not in self.__prepetual_isolates]

            self.clusters = self.__build_clusters()
            self.clusters.index.name = None
            self.clusters.columns = ['Plasmids']
        return self.clusters

    def __build_clusters(
        self
    ):
        return self.membership.groupby('Cluster').apply(lambda x: x.index.to_list())

    def __build_membership(
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

    def drop_isolates(
        self
    ):
        counts = self.cluster_counts.to_frame()
        counts.index.name = None
        counts.columns = ['Counts']
        clusters = self.clusters[counts['Counts'] > 1].to_frame()
        clusters.columns = ['Plasmids']

        membership = self.__build_membership(clusters)

        return clusters.reset_index(drop=True), membership

class Splitter:

    '''
    Class for big cluster subdivision. Takes the clusters DataFrame (cluster number vs cluster accessions
    in each cluster), as well as the full affinity matrix DataFrame.
    '''

    def __init__(
        self,
        clusters: pd.DataFrame,
        affinity: pd.DataFrame,
        distance_threshold: float = 0.8,
        linkage: str = 'average',
        n_clusters: Union[int, None] = None,
        max_size: int = 80
    ):
        self.clusters = clusters
        self.affinity = affinity
        self.linkage = linkage
        self.n_clusters = n_clusters
        self.distance_threshold = distance_threshold
        self.max_size = max_size

    def __split_cluster(
        self,
        cluster_number
    ):
        plasmids = self.clusters['Plasmids'].loc[cluster_number]
        cluster_affinity = self.affinity.loc[plasmids][plasmids]

        agg_clust = AgglomerativeClustering(
            n_clusters = self.n_clusters,
            affinity = 'precomputed',
            compute_full_tree = True,
            linkage = self.linkage,
            distance_threshold = self.distance_threshold
        )
        new_clusters = agg_clust.fit_predict((1-cluster_affinity).values)
        return pd.DataFrame(
            [
                [list(np.array(plasmids)[new_clusters == x])]
                for x in np.unique(new_clusters)
            ],
            index = [
                str(cluster_number)+f'.{n}'
                for n in np.unique(new_clusters)
            ],
            columns = ['Plasmids']
        )

    def split(
        self
    ):
        for cluster in self.clusters.index:
            if len(self.clusters['Plasmids'].loc[cluster]) > self.max_size:
                self.clusters = pd.concat(
                    [self.clusters, self.__split_cluster(cluster)]
                )
                self.clusters.drop(cluster, inplace = True)
            else:
                pass
        
        return self.clusters

# %%
