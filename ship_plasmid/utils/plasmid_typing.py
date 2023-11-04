# Class to display plasmid networks colored by given labels (i.e., plasmid types)
import pandas as pd
import os
import numpy as np
from pyvis.network import Network
from typing import Tuple, Union, Iterable
from numpy.typing import ArrayLike
import matplotlib.cm as cm
import webbrowser
from time import sleep
from copy import deepcopy
from sklearn.metrics import silhouette_score, calinski_harabasz_score, adjusted_mutual_info_score

class LabeledNetwork:
    '''
    PyVis Graph with nodes colored according to a user-specified mapping.
    '''

    def __init__(
        self,
        figsize: Tuple[str, str] = ('700px', '700px'),
        distance_threshold: float = .5,
        colors = None
    ):
        self.figsize = figsize
        self.network = Network(
            figsize[1],
            figsize[0],
            bgcolor = '#FFFFFF'
        )
        self.__threshold = distance_threshold
        self.colors = colors

    def __rbga_to_hex(
        self,
        rgba_: Iterable
    ) -> str:

        rgba = rgba_
        # Convert to [0, 255]
        rgba *= 255
        rgba = rgba[:-1].astype(int)
        rgba = np.minimum(255, rgba)

        return "#{0:02x}{1:02x}{2:02x}".format(*rgba)

    def __join_colors(
        self,
        colors: Iterable
    ):
        return np.mean(np.vstack(colors), axis=0)

    def __process_labels(
        self,
        labels: pd.Series,
    ):
        '''
        Calculates the number of labels and assigns colors to each one of them.
        Converts the labels Pandas Series into a one-hot encoding.
        
        If more than a label is assigned to the same sample, it interpolates the
        respective colors and displays a color halfway between them in the RGB
        space.

        Parameters
        ----------

        labels: Pandas Series
            Series with the label for each sample. May contain more than one entry
            for each sample, but needs at least one entry per sample.
        '''

        unique_labels = np.unique(labels.values)
        samples = np.unique(labels.index)
        self.n_samples, self.n_labels = len(samples), len(unique_labels)
        self.encoded = pd.DataFrame(
            np.zeros(
                (self.n_samples, self.n_labels),
                dtype = int
            ),
            columns = unique_labels,
            index = samples
        )

        for label in unique_labels:
            for sample in samples:
                if np.any(
                    labels.loc[sample] == label
                ):
                    self.encoded.loc[sample][label] = 1

        if self.colors is None:
            self.colormap = pd.Series({
                label_: color_ 
                for label_, color_ in zip(
                    unique_labels,
                    cm.rainbow(np.linspace(0,1,self.n_labels))
                )
            })
        else:
            self.colormap = pd.Series({
                label_: color_ 
                for label_, color_ in zip(
                    unique_labels,
                    self.colors
                )
            })
            
        self.colors = self.__assign_colors()

    def __assign_colors(
        self
    ):
        self.colors = pd.Series(
            np.zeros(self.n_samples),
            index = self.encoded.index
        )
        label_names = self.encoded.columns.to_numpy()
        for sample in self.encoded.index:
            # Check if the sample has only one label
            sample_labels = self.encoded.loc[sample]
            if sample_labels.sum() == 1:
                self.colors.loc[sample] = self.__rbga_to_hex(
                    deepcopy(self.colormap[
                        label_names[sample_labels.values.astype(bool)]
                    ].values[0])
                )
            else:
                self.colors.loc[sample] = self.__rbga_to_hex(
                    self.__join_colors(
                        [
                            deepcopy(self.colormap[x])
                            for x in label_names[sample_labels.values.astype(bool)]
                        ]
                    )
                )
        
        return self.colors

    def from_distance(
        self,
        M_: pd.DataFrame,
        labels: pd.Series
    ):
        '''
        Builds a network from a distance matrix.

        Parameters
        ----------

        M: Pandas DataFrame.
            Distance matrix. Index and columns will be used as node names.

        labels: Pandas Series
            Series with the label for each sample. May contain more than one entry
            for each sample, but needs at least one entry per sample.
        '''
        
        self.__process_labels(labels)
        M = 1-M_/np.max(M_.values)
        self.affinity = M
        self.distance = M_
        self.threshold = self.__threshold
        for sample in M.index:

            if type(labels.loc[sample]) == pd.Series:
                title = f"{sample} | Type {','.join(labels.loc[sample].values.flatten())}"
            else:
                title = f"{sample} | Type {labels.loc[sample]}"
            self.network.add_node(
                sample,
                label = sample,
                color = self.colors[sample],
                title = title
            )

        # Reduce number of iterations by dropping nodes with less than 
        # two entries above the threshold
        filtered_M = M[
            np.sum(
                M.values>self.threshold,
                axis = -1
            ) > 1
        ]
        filtered_M = filtered_M[filtered_M.index]
        n_samples = len(filtered_M.index)
        idx = filtered_M.index.to_list()
        for i in range(n_samples):
            for j in range(i+1, n_samples):
                value = filtered_M.iloc[i,j]
                if value > self.threshold:
                    self.network.add_edge(
                        idx[i],
                        idx[j],
                        value = value
                    )
        self.network.toggle_physics(False)
        self.network.show_buttons(['physics', 'nodes'])

    def show(
        self,
        path: Union[None, str] = None
    ):
        if path is None:
            path_ = os.path.join(
                os.getcwd(),
                'Typing-Graph-tmp.html'
            )
        else:
            path_ = path

        self.network.show(path_)

        webbrowser.open(
            'file://///wsl.localhost/Ubuntu'+path_,
            new = 1
        )

        sleep(5)

        if path is None:
            os.remove(path_)

    def get_metrics(
        self,
        labels: pd.Series
    ):
        metrics_f = {
            'Silhouette Score': silhouette_score,
            'Adjusted Mutual Information': adjusted_mutual_info_score
        }
        self.metrics = {k: {} for k in metrics_f}


        dist = deepcopy(self.distance.values)
        dist[self.affinity.values < self.threshold] = 1e20
        for k in self.encoded:
            y_true = self.encoded[k].loc[labels.index].values.flatten()
            y_pred = labels.values
            self.metrics['Adjusted Mutual Information'][k] = adjusted_mutual_info_score(
                y_true,
                y_pred
            )
            if len(np.unique(y_true)) > 1:
                self.metrics['Silhouette Score'][k] = silhouette_score(
                    dist,
                    labels = y_true,
                    metric = 'precomputed',
                    random_state=23
                )
            else:
                self.metrics['Silhouette Score'][k] = np.nan

        return pd.DataFrame(self.metrics)

# %%
