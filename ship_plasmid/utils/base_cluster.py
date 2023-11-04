'''
Base class for clustering.
'''

from typing import Union, Tuple, Any
from pandas import DataFrame, Series
from numpy import array
import numpy as np
from scipy import sparse
from abc import abstractmethod, ABC
from itertools import product

def check_type_set(
    a,
    b
) -> Tuple[set, set]:
    '''
    Check input type and try to cast to set.
    '''

    if type(a) == set and type(b) == set:
        return a, b
    else:
        try:
            return set(a), set(b)
        except:
            raise ValueError('Invalid input type for vectors.') 

class Clustering(ABC):

    def __init__(
        self,
        return_type: str = 'numpy'
    ):
        '''
        Base clustering class for plasmid clustering through gene annotations and
        sequences.

        Parameters
        ----------

        return_type: 'numpy' or 'pandas'
            If 'numpy', the output of predict will be a numpy array. If 'pandas', the
            output will be a Pandas Series.
        '''
        self.return_type = return_type
        super().__init__()

    def __check_pandas_order(
        self,
        X: DataFrame
    ) -> DataFrame:
        '''
        Reorders a Pandas DataFrame so that the index and columns are aligned.
        Keeps column order.
        '''

        if not len(
            set(X.index).intersection(X.columns)
        ) == len (
            set(X.index).union(X.columns)
        ):
            raise ValueError('The passed index and column names do not match.')
        else:
            return X.loc[X.columns]

    def __to_numpy(
        self,
        X: Union[array, DataFrame, Series]
    ) -> array:
        '''
        Converts a Pandas DataFrame or Series input to a 1D Numpy array, if needed. Stores
        index and column information.
        '''

        if type(X) == DataFrame:
            self.columns = X.columns.to_list()
            self.index = X.index.to_list()
            return X.values[:,0]

        elif type(X) == Series:
            self.index = X.index.to_list()
            self.columns = X.index.to_list()
            return X.values

        else:
            if len(X.shape) > 1:
                if X.shape[1] > 1 and len(X.shape) > 2:
                    raise ValueError(f'Input vector should be 1D or 2D with dimension 1 of size 1. Received input with shape {X.shape}.')
                else:
                    X = X[:,0]
            self.columns = None
            self.index = None
            return X

    def __build_affinity_matrix(
        self,
        X: Union[array, DataFrame, Series],
    ):
        X = self.__to_numpy(X)
        n_samples = len(X)
        self.affinity = np.zeros((n_samples, n_samples))

        for i in range(n_samples):
            for j in range(n_samples):
                self.affinity[i, j] = self.sim_f(X[i], X[j])
    
    def fit(
        self,
        X: Union[array, Series, DataFrame],
        affinity_matrix = None
    ):

        if affinity_matrix is not None:
            self.affinity = affinity_matrix
        else:
            self.X = X
            self.__build_affinity_matrix(X)
        self.clusters = self.__cluster(self.affinity)
        if type(X) == DataFrame or type(X) == Series:
            self.index = X.index.to_list()

    def fit_predict(
        self,
        X: Union[array, Series, DataFrame],
        affinity_matrix = None
    ) -> Union[array, Series]:

        self.fit(X, affinity_matrix)
        
        if self.return_type == 'numpy':
            return self.clusters
        elif self.return_type == 'pandas':
            return self.__predictions_to_pandas(self.clusters)
        else:
            raise ValueError(f'Invalid cluster output type. \'return_type\' set to {self.return_type}, but only \'numpy\' and \'pandas\' are valid.')

    def build_adjacency_matrix(
        self,
        threshold: float = .5,
        as_sparse: bool = True
    ) -> array:
        '''
        Builds an adjacency matrix by thresholding the affinity matrix.
        Can return as SciPy sparse matrix, if 'as_sparse' is set to
        True.
        '''
        self.adjacency = self.affinity > threshold
        self.adjacency = self.adjacency.astype(int)

        if as_sparse:
            self.adjacency = sparse.csr_matrix(self.adjacency)
        return self.adjacency

    def grid_search(
        self,
        X: Union[array, Series, DataFrame],
        params: dict,
        affinity_matrix = None,
        refit: bool = True
    ) -> dict:
        '''
        Repeats clustering with several parameters. Selects the parameters yielding
        the best (highest) fitness score. The parameters to search should be passed 
        as a dictionary with lists as values and parameter names as keys. The class
        methods score() and set_params() should be implemented in derived classes.
        '''
        self.search_scores = []
        best_score = -np.inf
        self.best_params = None
        self.X =X
        for iparams in product(*params.values()):
            iparams = {
                k: v
                for k, v in zip(
                    params.keys(),
                    iparams
                )
            }
            self.set_params(iparams)
            self.fit(X, affinity_matrix)
            self.search_scores.append((iparams, self.score(X)))
            if self.search_scores[-1][1] > best_score:
                best_score = self.search_scores[-1][1]
                self.best_params = iparams
            print(f'Finished clustering with parameters {iparams} | Score: {self.search_scores[-1][1]} | Best Score: {best_score}')
        
        if refit:
            self.set_params(self.best_params)
            return self.fit_predict(X, affinity_matrix)

    @abstractmethod
    def score(
        self,
        X
    ) -> float:
        pass

    @abstractmethod
    def set_params(
        self,
        params
    ):
        pass

    @abstractmethod      
    def __cluster(
        self,
        affinity: array
    ) -> array:
        '''
        Template class for clustering implementation, from the affinity matrix.
        Should be implemented in inherited classes.
        '''
        pass
    
    @abstractmethod
    def __predictions_to_pandas(
        self,
        y: Any
    ) -> DataFrame:
        '''
        Template class to convert the output of clustering inference to a Pandas
        DataFrame.
        '''
        pass


# %%
