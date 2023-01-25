'''
Uses KDE and MLE to estimate the thresholds for recombination/mutation/IS synteny blocks.
'''


from abc import ABC, abstractmethod
from typing import Callable, Iterable, Tuple
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from utils.files import find_annotation_paths
import yaml
from BCBio.GFF import parse

class BasePDF(ABC):

    def __init__(self) -> None:
        super().__init__()
        self.values = np.array(self.get_values())
        self.kde = self.build_pdf()
        self.XLIM = 10000
    
    @abstractmethod
    def get_values(
        self
    ):
        pass

    def build_pdf(
        self
    ):
        self.kde = gaussian_kde(self.values)
        return self.kde

    def max_prob(self):
        xs = np.arange(max(self.values))
        probs = self.kde(xs)
        best_len = np.argmax(probs)
        return xs[best_len], np.max(probs)

    def pdf(self, length):
        return self.kde(length)

    def __call__(
        self,
        length
    ):
        try:
            return self.pdf(length)
        except:

            return [self.pdf(x) for x in length]

    def plot(self):
        
        with plt.style.context('ggplot'):
            f, ax = plt.subplots(
                1,2,
                gridspec_kw={'width_ratios': [3, 1]},
                constrained_layout = True
            )
            ax[0].hist(self.values[self.values<self.XLIM], bins = 50)
            ax01 = ax[0].twinx()
            ax01.plot(self(np.arange(self.XLIM)), color = 'dodgerblue')
            ax[1].boxplot(
                self.values,
            )
            plt.show()


class MutationSizeProbability(BasePDF):

    def __init__(
        self,
        path_to_annotations: str,
        plasmid_ids: Iterable
    ) -> None:
        
        self.__paths = find_annotation_paths(
            plasmid_ids,
            path_to_annotations,
            format = '.gff'
        )
        super().__init__()

    def get_values(
        self
    ) -> list:
        values = []
        for plasmid_path in self.__paths:
            content = next(parse(plasmid_path)).features
            values += [
                x.location.end.position - x.location.start.position
                for x in content
            ]
        return values

class MobileElementSizeProbability(BasePDF):

    def __init__(
        self,
        path_to_isfinder: str
    ):

        self.organisms = [
            'Escherichia coli',
            'Staphylococcus aureus',
            'Klebsiella pneumoniae',
            'Enterococcus faecalis',
            'Acinetobacter baumannii',
            'Pseudomonas aeruginosa'
        ]
        self.path_to_isfinder = path_to_isfinder
        super().__init__()

    def get_values(
        self
    ) -> list:

        values = []
        self.sizes_by_organism = {}
        for organism in self.organisms:

            organism_code = organism[0]+'_'+organism.split()[1]
            with open(
                os.path.join(
                    self.path_to_isfinder,
                    f'{organism_code}_ISFinder.txt'
                ),
                'r'
            ) as instream:
                content = instream.read().split(organism)

            self.sizes_by_organism[organism] = np.array([int(x.split('\t')[1].split()[0]) for x in content[1::]])
            values = np.hstack([self.sizes_by_organism[organism], values])

        return values

class RecombinationSizeProbabiliy(MutationSizeProbability):

    def get_values(self) -> list:
        self.SAMPLE_SIZE = 10000
        values = super().get_values()
        # Sample 20 genes at random to build a selection
        # of total HR section sizes
        rng = np.random.default_rng(seed = 23)
        values = np.sum(
            rng.choice(
                values,
                size = (self.SAMPLE_SIZE,20)
            ),
            axis=-1
        )
        return values

def func_intersection(
    x: Iterable,
    f1: Callable,
    f2: Callable,
    left: float,
    right: float,
    tol: float = 1e-10
) -> Tuple[float, float]:

    x = np.array(x)
    x = x[x<right]
    x = x[x>left]

    y1 = np.array([f1(xx) for xx in x])
    y2 = np.array([f2(xx) for xx in x])

    diff = np.abs(y2 - y1)
    intercept_x = x[np.argmin(diff)]
    intercept_y = np.min(diff)

    if intercept_y > tol:
        return None
    else:
        return (intercept_x, intercept_y)

if __name__ == '__main__':
    with open('configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/clustering_config.yaml', 'r') as config_file:
        cluster_config = yaml.load(config_file, Loader=yaml.Loader)
    pIS = MobileElementSizeProbability(
        'data/ISFinder Databases'
    )
    pMut = MutationSizeProbability(
        data_config['paths']['annotations'],
        cluster_config['accessions']
    )
    pHR = RecombinationSizeProbabiliy(
        data_config['paths']['annotations'],
        cluster_config['accessions']
    )

    with plt.style.context('ggplot'):
        fig, ax = plt.subplots(
            1,1,
            constrained_layout = True,
            figsize = (10,7)
        )
        xs = np.arange(1, 20000, 2)
        mut_is = func_intersection(
            xs,
            pMut,
            pIS,
            10, 5000,
            tol = 1e-5
        )
        is_hr = func_intersection(
            xs,
            pHR,
            pIS,
            2000, 20000,
            tol = 1e-5
        )

        print(f'Mutation/IS threshold: {int(mut_is[0])} bp.')
        print(f'IS/Recombination threshold: {int(is_hr[0])} bp.')

        for label, object_ in zip(
            ['Mutation', 'IS', 'Recombination'],
            [pMut, pIS, pHR]
        ):
            ax.fill_between(
                xs,
                np.zeros_like(xs),
                object_(xs),
                label = label,
                alpha=.7
            )

        ax.vlines(
            [i[0] for i in [mut_is, is_hr]],
            0,
            0.0015,
            ls = 'dashed',
            alpha = .7,
            color = 'black'
        )

        ax.set_xlabel('Sequence Length (bp)')
        ax.set_ylabel('Probability')
        ax.legend()
        ax.grid(False)
        ax.set_ylim(0, 0.0015)
        ax.set_xlim(0, 18000)
    plt.show()

# %%
