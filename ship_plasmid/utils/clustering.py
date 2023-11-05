'''
Implements classes useful for plasmid clustering and to build Jaccard-based
similarity networks.

NOTE: get_affinity_matrix() actually returns a distance matrix.
'''

from typing import Tuple, Union, Iterable
import numpy as np
import pandas as pd
from markov_clustering import run_mcl, get_clusters, modularity
from ship_plasmid.utils.jaccard import jaccard
import numpy.typing as npt
from copy import deepcopy
import os
from sklearn.cluster import AgglomerativeClustering
import joblib
from scipy import sparse
from abc import abstractmethod, ABC
from itertools import product
import matplotlib.pyplot as plt
import networkx as nx
import networkx.algorithms.community as nx_comm

class PlasmidClusters(ABC):

    def __init__(
        self,
        inflation: Union[float, Iterable] = 1.4,
        reg_factor: float = 1.0,
        aggregate: bool = True,
        min_cluster_size: int = 10,
        min_n_core: int = 5,
        drop_isolates = True,
        score: str = 'core',
        core_threshold: float = 0.7,
        resolution: float = 0.8,
        n_clusters: Union[int, None] = None,
        distance_threshold: float = 0.8,
        max_size: int = 80,
        linkage: str = 'complete'
    ):
        self.inflation = inflation
        self.reg_factor = reg_factor
        self.__grid_search = (type(inflation) == list)
        self.aggregate = aggregate
        self.min_cluster_size = min_cluster_size
        self.min_n_core = min_n_core
        self.drop_isolates = drop_isolates
        self.score_ = score
        self.core_threshold = core_threshold
        self.resolution = resolution
        self.max_size = max_size
        self.distance_threshold = distance_threshold
        self.n_clusters = n_clusters
        self.linkage = linkage

    def cluster(
        self,
        affinity: pd.DataFrame
    ) -> pd.DataFrame:

        self.affinty = affinity
        self.mcl = MarkovClustering(
            inflation = self.inflation
        )
        self.clusters = self.mcl.fit_predict(affinity)

        self.aggregation = CoreAggregate(
            min_n_core = self.min_n_core,
            min_cluster_size = self.min_cluster_size
        )
        self.clusters = self.aggregation.cluster(
            self.clusters,
            self.proteins
        )

        self.splitter = Splitter(
            clusters = self.clusters,
            affinity = self.affinty,
            distance_threshold = self.distance_threshold,
            linkage = self.linkage,
            n_clusters = self.n_clusters,
            max_size = self.max_size
        )
        self.clusters = self.splitter.split()

        if self.drop_isolates:
            self.clusters, idx = drop_isolates(self.clusters, return_idx=True)
            self.membership = get_membership(self.clusters)
            self.affinty = self.affinty.loc[self.membership.index][self.membership.index]
        return self.clusters

    def score(
        self
    ) -> float:
        if self.score_ == 'core':
            n_core = []
            for cluster in self.clusters.index:
                plasmids = self.clusters.loc[cluster]['Plasmids']
                proteins = self.proteins.loc[plasmids]['Proteins']
                unique_proteins = proteins.apply(lambda x: np.unique(x)).values
                unique_proteins = np.hstack(unique_proteins)
                unique_proteins, counts = np.unique(unique_proteins, return_counts=True)
                n_plasmids = len(plasmids)
                threshold = n_plasmids * self.core_threshold
                n_core.append(len(counts[counts>=threshold]))

            self.n_core = np.array(n_core)
            self.n_zero_core = sum(self.n_core==0)
            print(f'Score | Number of clusters without core genes: {self.n_zero_core}.')
            return -(len(self.clusters)*self.resolution + self.n_zero_core)
        elif self.score_ == 'modularity':

            return nx_comm.modularity(
                nx.Graph(self.mcl._MarkovClustering__mcl_result),
                self.mcl._MarkovClustering__np_clusters,
                resolution = self.resolution
            ) - self.reg_factor * len(self.mcl._MarkovClustering__np_clusters)

        else:
            raise TypeError('Invalid score type.')

    def grid_search(
        self,
        inflations: Iterable,
        affinity: pd.DataFrame,
        refit: bool = True,
        verbose: bool = True
    ):
        inflations_copy = deepcopy(inflations)
        self.search_scores = {}
        best_score, best_params = -np.inf, None
        for inflation in inflations:
            affinity_copy = deepcopy(affinity)
            self.inflation = inflation
            self.cluster(affinity_copy)
            self.search_scores[inflation] = self.score()
            if self.search_scores[inflation] > best_score:
                best_score = self.search_scores[inflation]
                best_params = inflation
            if verbose:
                print(f'MCL Grid Search | Inflation: {inflation} | Score: {self.search_scores[inflation]} | Best Inflation: {best_params} | Best Score: {best_score}')

        if refit:
            self.inflation = best_params
            self.cluster(affinity)

        self.inflation = inflations_copy
        if refit: return self.clusters

    def save(
        self,
        path: str
    ):
        self.clusters.to_csv(
            os.path.join(
                path,
                f'Clusters.tsv'
            ),
            sep = '\t'
        )
        self.affinity.to_csv(
            os.path.join(
                path,
                f'Affinity.tsv'
            ),
            sep = '\t'
        )
        joblib.dump(
            self.clusters,
            os.path.join(
                path,
                f'Clusters.pkl'
            ),
            compress = 3
        )
    

    @abstractmethod
    def build_affinity_matrix(
        self
    ) -> pd.DataFrame:
        pass

class MarkovClustering():

    def __init__(
        self,
        inflation: float,
    ):
        '''
        Clusters a graph from an affinity matrix, using Markov Clustering.

        Arguments
        ----------

        inflation: float
            Inflation parameter for Markov Clustering. Higher values mean higher granularity.        
        '''

        self.inflation = inflation

    def fit(
        self,
        X: Union[pd.DataFrame, npt.ArrayLike]
    ) -> pd.DataFrame:
        '''
        Parameters
        ----------

        X: Numpy Array or Pandas DataFrame.
            Affinity matrix in which to perform clustering.

        Returns
        -------

        clusters: Pandas DataFrame.
        '''

        self.adjacency = X
        self.__mcl_result = run_mcl(
            self.adjacency.values,
            inflation = self.inflation
        )
        # Convert to Numpy and store index in self.index
        if type(X) == pd.DataFrame: X = self.__pandas_to_numpy(X)
        self.clusters = get_clusters(self.__mcl_result)
        self.__np_clusters = deepcopy(self.clusters)
        self.clusters = self.__clusters_to_pandas(self.clusters)

    def fit_predict(
        self,
        X: Union[pd.DataFrame, npt.ArrayLike]
    ) -> pd.DataFrame:
        '''
        Parameters
        ----------

        X: Numpy Array or Pandas DataFrame.
            Affinity matrix in which to perform clustering.

        Returns
        -------

        clusters: Pandas DataFrame.
        '''
        self.fit(X)
        return self.clusters

    def __pandas_to_numpy(
        self,
        X: pd.DataFrame
    ):
        self.index = X.index
        return X.values

    def __clusters_to_pandas(
        self,
        clusters
    ) -> pd.DataFrame:
        return pd.DataFrame(
            [  
                [[
                    self.index[plasmid_number]
                    for plasmid_number in cluster
                ]]
                for cluster in clusters
            ],
            index = np.arange(len(self.clusters)),
            columns = ['Plasmids']
        )

class AffinityMatrix:

    def __init__(
        self
    ):
        pass

    def from_jaccard(
        self,
        proteins: pd.DataFrame
    ) -> pd.DataFrame:

        self.proteins = proteins.values[:,0]
        n_samples = len(proteins)
        self.affinity = np.zeros((n_samples, n_samples))

        for i in range(n_samples):
            for j in range(n_samples):
                self.affinity[i, j] = jaccard(self.proteins[i], self.proteins[j])

        self.affinity = pd.DataFrame(
            self.affinity,
            index = proteins.index,
            columns = proteins.index
        )

        return self.affinity

    def from_ani(
        self,
        plasmids: list,
        sequence_path: str,
        output_path: str,
        k: int = 16,
        min_ani_length: float = 0.2,
        load: bool = True
    ) -> pd.DataFrame:
        raise NotImplementedError

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

            self.clusters = self.__build_clusters().to_frame()
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

class CoreAggregate:

    def __init__(
        self,
        min_n_core: int,
        min_cluster_size: int,
        verbose: bool = True
    ):
        self.min_cluster_size = min_cluster_size
        self.min_n_core = min_n_core
        self.n_core = None
        self.core_genes = None
        self.verbose = verbose

    def __core_genes_in_cluster(
        self,
        plasmids: Iterable,
        proteins: pd.DataFrame
    ) -> Tuple[int, list]:

        proteins_in_cluster = proteins['Proteins'].loc[plasmids]
        core_genes = set(proteins_in_cluster.iloc[0])
        for plasmid in proteins_in_cluster.iloc[1:]:
            core_genes = core_genes.intersection(plasmid)

        return len(core_genes), list(core_genes)

    def __core_genes_in_all_clusters(
        self,
        clusters: pd.DataFrame,
        proteins: pd.DataFrame
    ) -> pd.Series:

        return pd.Series(
            [
                self.__core_genes_in_cluster(
                    clusters.loc[cluster]['Plasmids'],
                    proteins
                )[1]
                for cluster in clusters.index
            ],
            index = clusters.index
        )

    def get_core(
        self,
        clusters: pd.DataFrame,
        proteins: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:

        self.initial_core_genes = self.__core_genes_in_all_clusters(
            clusters,
            proteins
        )

        self.core_genes = np.zeros(
            (len(self.initial_core_genes), len(self.initial_core_genes)),
            dtype = object
        )
        self.n_core = np.zeros(
            (len(self.initial_core_genes), len(self.initial_core_genes))
        )

        for n, i in enumerate(self.initial_core_genes.index):
            for m, j in enumerate(self.initial_core_genes.index):
                self.core_genes[n, m] = set(
                    self.initial_core_genes.loc[i]
                ).intersection(
                    self.initial_core_genes.loc[j]
                )

                self.n_core[n, m] = len(self.core_genes[n, m])
        
        for n in range(self.n_core.shape[0]):
            self.n_core[n,n] = 0

        self.n_core = pd.DataFrame(
            self.n_core,
            index = self.initial_core_genes.index,
            columns = self.initial_core_genes.index
        )

        self.core_genes = pd.DataFrame(
            self.core_genes,
            index = self.initial_core_genes.index,
            columns = self.initial_core_genes.index
        )
        
        return self.n_core, self.core_genes

    def cluster(
        self,
        clusters: pd.DataFrame,
        proteins: pd.DataFrame,
        max_iter: int = 20
    ) -> pd.DataFrame:
        
        new_clusters = []
        # Clusters that do not fulfill the requirements for joining
        unjoined = []
        any_join = True
        n_iter = 0
        
        while any_join and n_iter < max_iter:
            n_core, core = self.get_core(
                clusters,
                proteins
            )
            n_plasmids = clusters['Plasmids'].apply(lambda x: len(x))
            any_join = False
            for cluster in clusters.index:
                if n_plasmids.loc[cluster] < self.min_cluster_size and not cluster in unjoined and cluster in clusters.index:

                    best_match_n = n_core.loc[cluster].max()
                    if best_match_n > self.min_n_core:
                        best_match = n_core.columns[np.argmax(n_core.loc[cluster])]
                        clusters = pd.concat(
                            [
                                clusters,
                                pd.DataFrame(
                                    [[clusters['Plasmids'].loc[cluster] + clusters['Plasmids'].loc[best_match]]],
                                    index = ['+'.join([str(cluster), str(best_match)])],
                                    columns = ['Plasmids']
                                )
                            ]
                        ).drop([cluster, best_match])

                        n_core, core = self.get_core(
                            clusters,
                            proteins
                        )
                        any_join = True
                        if self.verbose:
                            print(f'CoreAggregate | Joining cluster {cluster} with cluster {best_match} | No. core genes: {best_match_n}')
                    else:
                        unjoined.append(cluster)
                else:
                    if not cluster in unjoined:
                        unjoined.append(cluster)

            n_iter += 1
        
        return clusters

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
        print(f'SPLITTER [utils.clustering.py Line 311]: Split cluster {cluster_number} into {len(np.unique(new_clusters))} new clusters.')
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

class JaccardClusters(PlasmidClusters):

    def __init__(
        self,
        path_to_proteins: str,
        similarity: int = 9,
        k: int = 5,
        inflation: Union[Iterable, float] = 1.4,
        reg_factor: float = 1,
        aggregate: bool = True,
        min_n_core: int = 10,
        min_cluster_size: int = 10,
        drop_isolates: bool = True,
        resolution: float = 0.8,        
        score: str = 'core',
        core_threshold: float = 0.7,
        n_clusters: Union[int, None] = None,
        distance_threshold: float = 0.8,
        max_size: int = 80,
        linkage: str = 'complete'
    ):
        self.path = path_to_proteins
        self.similarity = similarity
        self.k = k
        super().__init__(
            inflation,
            reg_factor,
            aggregate,
            min_cluster_size,
            min_n_core,
            drop_isolates,
            score,
            core_threshold,
            resolution,
            n_clusters,
            distance_threshold,
            max_size,
            linkage
        )

        self.__get_proteins()

    def __get_proteins(self):
        self.proteins = joblib.load(
            os.path.join(
                self.path,
                f'Plasmids with Clustered Proteins_s{self.similarity}_k{self.k}.pkl'
            )
        )
    
    def build_affinity_matrix(self) -> pd.DataFrame:

        self.__affinity_builder = AffinityMatrix()
        self.affinity = self.__affinity_builder.from_jaccard(
            proteins = self.proteins
        )

        return self.affinity


def drop_isolates(
    clusters,
    return_idx = False
):
    counts = clusters['Plasmids'].apply(lambda x: len(x))
    clusters = clusters['Plasmids'][counts > 1].to_frame()
    print(f'Drop Isolates | Dropping {sum(counts < 2)} isolates.')
    clusters.columns = ['Plasmids']
    if return_idx:
        original_idx = clusters.index.to_list()
    clusters.index = np.arange(len(clusters))

    if return_idx:
        return clusters, original_idx
    else:
        return clusters

def get_membership(
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

class AdaptativeAgglomerativeClustering(AgglomerativeClustering):

    def __init__(
        self, 
        n_clusters=2, 
        *, 
        affinity="euclidean", 
        memory=None, 
        connectivity=None, 
        compute_full_tree="auto", 
        linkage="ward", 
        distance_threshold=None,
        distance_ratio = 0.6, 
        compute_distances=False
    ):
        self.distance_ratio = distance_ratio
        super().__init__(
            n_clusters, 
            affinity=affinity, 
            memory=memory, 
            connectivity = connectivity, 
            compute_full_tree = compute_full_tree, 
            linkage = linkage, 
            distance_threshold = distance_threshold, 
            compute_distances = compute_distances
        )

    def fit(
        self,
        X,
        y = None
    ):
        if self.distance_ratio is not None:
            self.distance_threshold = self.distance_ratio * np.max(X)
        
        return super().fit(X)
# %%
