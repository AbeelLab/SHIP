'''
Script for extracting statistics about tables with conserved regions.
'''
import warnings
import pandas as pd
import numpy as np
import joblib
import yaml
from utils.files import get_amr
from matplotlib import pyplot as plt
import os
from copy import deepcopy

def line_splitter(
    line: pd.Series,
    on: str = 'Plasmids'
):
    line = line.to_frame().transpose()
    out_df = pd.concat(
        [line for _ in range(len(line.iloc[0][on]))]
    )
    for i, _ in enumerate(out_df.index):
        out_df.iloc[i][on] = line.iloc[0][on][i]

    return out_df

def is_nested(
    table
):
    all_genes = deepcopy(table['Genes'].values)
    all_plasmids = deepcopy(table['Plasmids'].values)
    is_nested = []

    nested_genes = []
    for genes, plasmids in zip(table['Genes'], table['Plasmids']):
        is_nested_it = np.any(
            [
                set(genes).issubset(gg)
                and set(pp).issubset(plasmids)
                for gg, pp in zip(all_genes, all_plasmids)
                if (gg != genes and gg not in nested_genes)
            ]
        )
        is_nested.append(is_nested_it)
        if is_nested_it:
            nested_genes.append(genes)

    return is_nested

if __name__ == '__main__':
    with open('./configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('./configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('./configs/motif_config.yaml', 'r') as config_file:
        motif_config = yaml.load(config_file, Loader=yaml.Loader)

    for k in phylo_config['output-paths']:
        phylo_config['output-paths'][k] = os.path.join(
            phylo_config['results-dir'],
            phylo_config['output-paths'][k]
        )

    info = pd.read_csv(
        data_config['paths']['info_df'],
        sep = '\t',
        index_col = 0
    )
    info.index = [x.split('.')[0] for x in info.index]
    amr = get_amr(
        info.index.to_list(),
        data_config['paths']['amr_hits']
    )

    motifs = joblib.load(
        os.path.join(
            motif_config['motif-finder-output-dir'],
            f'Motif_Finder_Results_Filtered.pkl'
        )
    ).drop('index', axis='columns')
    motifs['Nested'] = is_nested(motifs)

    for n, motifs_ in enumerate([
        deepcopy(motifs),
        deepcopy(motifs[~motifs['Nested']])
    ]):

        if n: print('\n:::::: NON-NESTED ::::::')
        else: print('\n:::::: NESTED ::::::')
        # Split lines in motifs table into 1 entry per plasmid and motif
        split_motifs = line_splitter(
            motifs_.iloc[0]
        ) 
        for line in motifs_.index.to_list()[1:]:
            split_motifs = pd.concat(
                [
                    split_motifs,
                    line_splitter(motifs_.loc[line])
                ]
            )

        split_motifs['Genes'] = split_motifs['Genes'].apply(
            lambda x: np.array(x)
        )
        split_motifs.reset_index(drop=False, inplace=True)
        split_motifs = split_motifs.join(info, on='Plasmids')

        #% How many cross-species transfers?
        query = split_motifs
        # Species per fragment
        n_species = query.groupby(
                ['index', 'Organism']
            ).size().groupby('index').size()

        total = len(motifs_)
        print(
            f'{len(n_species[n_species==2])} ({len(n_species[n_species==2])/total*100:.3}%) motifs are present in two species.'
        )
        # Plot species distribution
        fragments_per_species = split_motifs.groupby(
            ['index', 'Organism']
        ).size().groupby('Organism').size()
        print('Number of fragments per species:')
        print(fragments_per_species)

        # Species co-occurrences
        multispecies = query.groupby(
                ['index', 'Organism']
            ).size().loc[n_species[n_species>1].index]
        species = np.unique(query['Organism'])
        cooccurrence = pd.DataFrame(
            np.zeros((len(species), len(species)), int),
            columns = species,
            index = species
        )
        for i in np.unique(multispecies.index.get_level_values(0)):
            # Check species with motif
            species_w_motif = multispecies.loc[i].index.to_numpy()
            cooccurrence[species_w_motif[0]][species_w_motif[1]] += 1
            cooccurrence[species_w_motif[1]][species_w_motif[0]] += 1

        print(f'Motif species co-occurrence:'.capitalize())
        print(cooccurrence)

        # % Multiresistant fragments
        # Split lines in motifs table into 1 entry per gene and motif
        amr_splits = line_splitter(
            motifs_.iloc[0],
            on = 'Gene Annotations'
        )
        for line in motifs_.index.to_list()[1:]:
            amr_splits = pd.concat(
                [
                    amr_splits,
                    line_splitter(motifs_.loc[line], on='Gene Annotations')
                ]
            )

        # Get unique annotations of AMR genes
        amr_annotations = np.unique(amr['Sequence name'])
        query = amr_splits
        query['Is AMR'] = query['Gene Annotations'].apply(
            lambda x: x in amr_annotations
        )
        query.reset_index(drop=False, inplace=True)

        for att in ['Class', 'Subclass']:
            def find_class(x):
                class_ = np.unique(amr[att][amr['Sequence name']==x])
                if len(class_): return class_[0]
                else: return 'None'
            query[att] = query['Gene Annotations'].apply(
                find_class
            )
        multiresistant = query[query['Is AMR']].groupby('index').size()
        multiresistant = multiresistant[multiresistant>1]
        n_multiresistant = len(multiresistant)
        
        print(
            f'{n_multiresistant} ({n_multiresistant/total*100:.3}%) motifs are multiresistant.'
        )

        res_genes = np.unique(query['Gene Annotations'])
        annotations = query[query['Is AMR']].groupby(['index', 'Gene Annotations']).size()[multiresistant.index]
        classes = query[query['Is AMR']].groupby(['index', 'Class']).size()[multiresistant.index]
        subclass = query[query['Is AMR']].groupby(['index', 'Subclass']).size()[multiresistant.index]

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            for table, table_name in zip(
                [annotations, classes, subclass],
                ['AMR gene', 'AMR class', 'AMR subclass']
            ):
                entries = np.unique(table.index.get_level_values(1))
                amr_cooccurrence = pd.DataFrame(
                    np.zeros((len(entries), len(entries)), int),
                    columns = entries,
                    index = entries
                )
                for i in np.unique(table.index.get_level_values(0)):
                    # Check species with motif
                    amr_genes = table.loc[i].index.to_numpy()
                    if len(amr_genes) == 1:
                        amr_cooccurrence[amr_genes[0]][amr_genes[0]] += 1
                    else:
                        for gene_a in amr_genes:
                            for gene_b in amr_genes:
                                if gene_a != gene_b:
                                    amr_cooccurrence[gene_a][gene_b] += 1
                    

                amr_cooccurrence = amr_cooccurrence.loc[
                    np.any(
                        amr_cooccurrence,
                        axis=0
                    )
                ]
                amr_cooccurrence = amr_cooccurrence[
                    amr_cooccurrence.index
                ]            
                print(f'Motif {table_name} co-occurrence:'.capitalize())
                #print(amr_cooccurrence) 

                with plt.style.context('ggplot'):

                    fontsize = 22*13/len(amr_cooccurrence)+4

                    f, ax = plt.subplots(
                        1,1,
                        constrained_layout = True,
                        figsize = (20,20)
                    )
                    ax.matshow(
                        amr_cooccurrence.values,
                        cmap='inferno'
                    )
                    ax.set_xticks(
                        np.arange(len(amr_cooccurrence))
                    )
                    ax.set_xticklabels(
                        amr_cooccurrence.index.to_list(),
                        rotation = -90,
                        fontsize = fontsize
                    )
                    ax.set_yticks(
                        np.arange(len(amr_cooccurrence))
                    )
                    ax.set_yticklabels(
                        amr_cooccurrence.index.to_list(),
                        fontsize = fontsize
                    )
                    ax.grid(False, which='major')
                    # Set grid
                    ax.set_xticks(
                        np.arange(len(amr_cooccurrence)) + .5,
                        minor = True,
                        lw = 1
                    )
                    ax.set_yticks(
                        np.arange(len(amr_cooccurrence)) + .5,
                        minor = True,
                        lw = 1
                    )
                    ax.grid(True, which='minor')

                    for i in np.arange(len(amr_cooccurrence)):
                        for j in np.arange(len(amr_cooccurrence)):
                            value_ = amr_cooccurrence.iloc[i,j]
                            if value_ > 0:
                                if value_ > 0.6*np.max(amr_cooccurrence.values):
                                    color_ = 'black'
                                else:
                                    color_ = 'white'
                                ax.text(
                                    i, j,
                                    str(value_),
                                    color = color_,
                                    va = 'center',
                                    ha = 'center',
                                    fontsize = fontsize
                                )

                    plt.show()

    #%% Length of motifs and number of plasmids
    lens = []
    n_plasmids = []
    with plt.style.context('ggplot'):
        f, ax = plt.subplots(
            2,2,
            constrained_layout = True,
            figsize = (20,20)
        )
        len_ax = ax[0,:]
        plm_ax = ax[1,:]
        for n, motifs_ in enumerate([
            motifs,
            motifs[~motifs['Nested']]
        ]):  
            lens.append(motifs_['Length'])
            n_plasmids.append(motifs_['Number of Plasmids'])
            color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

        for i, name_ in enumerate(['All', 'Non-nested']):
            len_ax[i].hist(
                lens[i],
                bins = np.arange(5, 10),
                lw = 1,
                edgecolor = 'white'
            )
            len_ax[i].set_title(
                ' '.join(
                    [
                        'Number of genes -',
                        name_, 
                        'regions'
                    ]
                ),
                fontsize = 30
            )
            len_ax[i].set_xticks(
                np.arange(5, 10)
            )
            len_ax[i].set_xticklabels(
                np.arange(5, 10),
                fontsize = 20
            )
            len_ax[i].set_yticklabels(
                len_ax[i].get_yticklabels(),
                fontsize = 20
            )
            len_ax[i].set_xlabel(
                'Number of genes',
                fontsize = 20
            )
            len_ax[i].set_ylabel(
                'Number of regions',
                fontsize = 20
            )

        for i, name_ in enumerate(['All', 'Non-nested']):
            plm_ax[i].hist(
                n_plasmids[i],
                bins = np.arange(3, 34, 1),
                lw = 1,
                edgecolor = 'white'
            )
            plm_ax[i].set_title(
                ' '.join(
                    [
                        'Number of plasmids with\nregion -',
                        name_, 
                        'regions'
                    ]
                ),
                fontsize = 30
            )
            plm_ax[i].set_xticks(
                np.arange(3, 34, 2)
            )
            plm_ax[i].set_xticklabels(
                np.arange(3, 34, 2),
                fontsize = 20
            )
            plm_ax[i].set_yticklabels(
                plm_ax[i].get_yticklabels(),
                fontsize = 20
            )
            plm_ax[i].set_xlabel(
                'Number of plasmids',
                fontsize = 20
            )
            plm_ax[i].set_ylabel(
                'Number of regions',
                fontsize = 20
            )

    plt.show()

# %%
