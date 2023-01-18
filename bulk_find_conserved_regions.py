'''
Finds all conserved regions containing AMR genes in a panplasmidome, without manual selection
of query nodes and paths.
'''

import os
import yaml
import joblib
import numpy as np
from utils.files import get_amr
from utils.motifs import BulkMotifFinder
import pandas as pd
from matplotlib import pyplot as plt
import argparse

parser = argparse.ArgumentParser(
    description = '''
Finds all conserved regions containing AMR genes in a panplasmidome, without manual selection
of query nodes and paths.
'''
)
parser.add_argument(
    '--min-dist',
    default = 0.1,
    nargs = 1,
    help = 'Minimum average plasmid distance for region inclusion. Default is 0.1.'
)
parser.add_argument(
    '--min_len',
    default = 5,
    nargs = 1,
    help = 'Minimum fragment length in CDS. Default is 5.'
)
parser.add_argument(
    '--max_len',
    default = 9,
    nargs = 1,
    help = 'Minimum fragment length in CDS. Default is 9.'
)
parser.add_argument(
    '--min_n',
    default = 3,
    nargs = 1,
    help = 'Minimum number of plasmids containing a region for inclusion. Default is 3'
)
parser.add_argument(
    '--out',
    default = './',
    nargs = 1,
    help = 'Output directory. Default is CWD.'
)
args = parser.parse_args() 

if __name__ == '__main__':

    with open('./configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('./configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)
    #with open('./configs/motif_config.yaml', 'r') as config_file:
    #    motif_config = yaml.load(config_file, Loader=yaml.Loader)

    motif_config = {
        'min-distance': args.min_dist[0],
        'min-length': args.min_len[0],
        'max-length': args.max_len[0],
        'min-n-plasmids': args.min_n[0],
        'motif-finder-output-dir': args.out[0]
    }

    for k in phylo_config['output-paths']:
        phylo_config['output-paths'][k] = os.path.join(
            phylo_config['results-dir'],
            phylo_config['output-paths'][k]
        )

    phylos = joblib.load(phylo_config['output-paths']['phylogenies'])

    with open(
        os.path.join(
            motif_config['motif-finder-output-dir'],
            f'Motif_Finder_Parameters.txt'
        ),
        'w'
    ) as outstream:
        outstream.write(
            '\n'.join(
                [
                    'Phylo Config:',
                    str(phylo_config),
                    'Motif Config:',
                    str(motif_config)
                ]
            )
        )

    for n, (cluster, phylo) in enumerate(phylos.items()):
        if cluster == 37:
            print(f'Finding motifs for cluster {cluster} ({n+1}/{len(phylos)})...')
            plasmid_ids = phylo.accessions
            amr_df = get_amr(
                plasmid_ids,
                data_config['paths']['amr_hits']
            )
            resistant_ids = np.unique(amr_df[amr_df['Gene symbol']!= 'Susceptible'].index)

            bmf = BulkMotifFinder(
                resistant_ids,
                phylo,
                data_config,
                phylo_config
            )
            result = bmf.search(
                motif_config['min-distance'],
                motif_config['min-length'],
                motif_config['max-length'],
                motif_config['min-n-plasmids']
            )
            bmf.save_report(
                os.path.join(
                    motif_config['motif-finder-output-dir'],
                    f'Motif_Finder_Results_Cluster_{cluster}.tsv'
                )
            )
            joblib.dump(
                result,
                os.path.join(
                    motif_config['motif-finder-output-dir'],
                    f'Motif_Finder_Results_Cluster_{cluster}.pkl'
                ),
                compress = 3
            )
    # Concatenate data for all clusters
    results = pd.concat(
        [
            joblib.load(
                os.path.join(
                    motif_config['motif-finder-output-dir'],
                    f'Motif_Finder_Results_Cluster_{cluster}.pkl'
                )
            )
            for cluster in phylos.keys()
        ]
    )
    results.to_csv(
        os.path.join(
            motif_config['motif-finder-output-dir'],
            'Motif_Finder_Results.tsv'
        ),
        sep = '\t'
    )
    joblib.dump(
        results,
        os.path.join(
            motif_config['motif-finder-output-dir'],
            'Motif_Finder_Results.pkl'
        ),
        compress = 3
    )

    # Add Nested column: true if a motif is a subset of another motif, with the same query gene

    results = joblib.load(
        os.path.join(
            motif_config['motif-finder-output-dir'],
            'Motif_Finder_Results.pkl'
        )
    )

    results = results.reset_index()
    results['Nested'] = [False] * len(results)

    for gene in np.unique(results['Query Gene ID']):
        
        all_plasmid_sets = results[results['Query Gene ID']==gene]['Genes'].values
        results['Nested'].loc[
            results[results['Query Gene ID']==gene].index
        ] = results['Genes'][results['Query Gene ID']==gene].apply(
            lambda x: np.any(
                [
                    len(set(pset).intersection(x)) == len(x)
                    for pset in all_plasmid_sets
                    if pset != x
                ]
            )
        )

    # Discard entries with the same gene set
    idx_keep = results['Genes'].astype('str').drop_duplicates(keep='first').index
    results = results.loc[idx_keep]

    results.to_csv(
        os.path.join(
            motif_config['motif-finder-output-dir'],
            'Motif_Finder_Results_Filtered.tsv'
        ),
        sep = '\t'
    )
    joblib.dump(
        results,
        os.path.join(
            motif_config['motif-finder-output-dir'],
            'Motif_Finder_Results_Filtered.pkl'
        ),
        compress = 3
    )

    def boxplot_column(
        results: pd.DataFrame,
        column: str
    ):
        with plt.style.context('ggplot'):
            fig, ax = plt.subplots(
                1,1,
                figsize = (7,7),
                constrained_layout = True
            )

            colors = [
                '#2176DD',
                '#DB840B',
                '#49BEE3',
                '#7F7F7F',
                '#64DB16'
            ]

            ax.set_facecolor('#FFFFFF')
            ax.boxplot(
                results[column],
                patch_artist = True,
                boxprops = {
                    'facecolor': colors[0],
                    'edgecolor': 'white'
                },
                capprops = {
                    'color': colors[0],
                    'lw': 5
                },
                whiskerprops = {
                    'color': colors[0],
                    'lw': 3
                },
                medianprops = {
                    'color': 'white',
                    'lw': 5
                },
                flierprops = {
                    'color': colors[0],
                },
                widths = .5
            )
            ax.set_xticks([])
            ax.set_title(
                column,
                fontsize = 20
            )
            ax.set_yticks(
                np.arange(
                    np.floor(ax.get_ylim()[0]),
                    np.ceil(ax.get_ylim()[1]),
                    10,
                    dtype = int
                )
            )
            ax.set_yticklabels(
                np.arange(
                    np.floor(ax.get_ylim()[0]),
                    np.ceil(ax.get_ylim()[1]),
                    10,
                    dtype = int
                ),
                fontsize = 17
            )
            ax.grid(
                True,
                axis = 'y',
                color = colors[3],
                lw = 1
            )
            ax.spines['left'].set_color(colors[3])
            ax.spines['left'].set_linewidth(3)
            ax.spines['bottom'].set_color(colors[3])
            ax.spines['bottom'].set_linewidth(3)
        
            plt.show()

    boxplot_column(results, 'Number of Plasmids')
    colors = [
                '#2176DD',
                '#DB840B',
                '#49BEE3',
                '#7F7F7F',
                '#64DB16'
            ]
    with plt.style.context('ggplot'):
        fig, ax = plt.subplots(
            1,1,
            figsize = (7,7),
            constrained_layout = True
        )

        ax.barh(
            [
                'Integrases',
                'Transposases'
            ],
            [
                results['Integrases'].apply(lambda x: len(x)>0).sum(),
                results['Transposases'].apply(lambda x: len(x)>0).sum()
            ],
            label = 'With',
            color = colors[1]
        )
        ax.barh(
            [
                'Integrases',
                'Transposases'
            ],
            [
                results['Integrases'].apply(lambda x: len(x)==0).sum(),
                results['Transposases'].apply(lambda x: len(x)==0).sum()
            ],
            left = [
                results['Integrases'].apply(lambda x: len(x)>0).sum(),
                results['Transposases'].apply(lambda x: len(x)>0).sum()
            ],
            label = 'Without',
            color = colors[0]
        )

        ax.set_facecolor('#FFFFFF')
        ax.spines['left'].set_color(colors[3])
        ax.spines['left'].set_linewidth(3)
        ax.spines['bottom'].set_color(colors[3])
        ax.spines['bottom'].set_linewidth(3)

        plt.plot()
