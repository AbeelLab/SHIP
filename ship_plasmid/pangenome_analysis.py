'''
Extracts stats about the pangenome, including local variability.
'''
import warnings
import networkx as nx
from pyvis.network import Network
from copy import deepcopy
import os
import yaml
import pandas as pd
import numpy as np
from ship_plasmid.utils.phylogeny import plot_pangenome
from ship_plasmid.utils.motifs import GeneMotifs
from ship_plasmid.utils.pangenome import PangenomeAnalysis
import joblib
from scipy.stats import mannwhitneyu
from ship_plasmid.utils.interpretation import get_protein_function, linebreak
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy.stats import linregress, logser, geom, poisson, zipf, chisquare, norm
from scipy.optimize import differential_evolution
import argparse

parser = argparse.ArgumentParser(
    description = '''
Extracts stats about the pangenome, including local variability.
'''
)
parser.add_argument(
    '--cluster',
    default = None,
    nargs = 1,
    help = 'Jaccard cluster for which to build a pangenome. If not specified, assumes that there is only one cluster.'
)
parser.add_argument(
    '--paper',
    action = 'store_const',
    const = True, 
    default = False, 
    help = 'Loads the Jaccard-based clusters used in the original paper.'
)
args = parser.parse_args() 

if __name__ == '__main__':
    with open('configs/data_config.yaml', 'r') as config_file:
        data_config = yaml.load(config_file, Loader=yaml.Loader)
    with open('configs/phylo_config.yaml', 'r') as config_file:
        phylo_config = yaml.load(config_file, Loader=yaml.Loader)

    if args.paper:
        phylo_config['results-dir'] = os.path.join(
            phylo_config['results-dir'],
            'Paper'
        )
    for k in phylo_config['output-paths']:
        phylo_config['output-paths'][k] = os.path.join(
            phylo_config['results-dir'],
            phylo_config['output-paths'][k]
        )

    phylos = joblib.load(phylo_config['output-paths']['phylogenies'])

    if args.cluster is not None:
        phylo = phylos[int(args.cluster[0])]
    else:
        phylo = phylos


    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pangenome = plot_pangenome(
            phylo.accessions,
            data_config,
            phylo_config
        )
    m = GeneMotifs(
        phylo, data_config, phylo_config, min_n_edges=0
    )
    m.from_ids(pangenome, multigraph = True)
    m.show_motifs()
    joblib.dump(
        deepcopy(m),
        os.path.join(
            phylo_config['results-dir'],
            f'GeneMotifs_cluster_{args.cluster}.pkl'
        ),
        compress=3
    )

    graph_analysis = PangenomeAnalysis(
        m.get_motif_graph_as_nx(),
        m.motif_graph,
        phylo_config
    )
    # Get degree and annotations in each community
    mean_degrees, annotations = graph_analysis.degree_report()

    mean_degrees_df = pd.DataFrame.from_dict(
        mean_degrees, 
        orient = 'index', 
        columns = ['Degree']
    )
    n_communities = max(annotations['Community'])
    # Add missing communities (< 5 nodes) to the mean_degrees_df
    missing = pd.DataFrame(
        [[None]] * (n_communities - len(mean_degrees) +1),
        index = [x for x in range(n_communities) if x not in mean_degrees_df.index],
        columns = ['Degree']
    )
    mean_degrees_df = pd.concat([mean_degrees_df, missing]).sort_index()
    annotations_df = deepcopy(annotations)
    annotations_df['Degree'] = mean_degrees_df.loc[annotations['Community']].values

    # Set annoatations_df index as (Annotation, Community) and group by annotation
    annotations_df.index = pd.MultiIndex.from_frame(annotations_df[['Annotation', 'Community']])
    annotations_df = annotations_df.drop(['Annotation', 'Community'], axis = 1).sort_index()

    # Perform Mann-Whitney U test on the average degree of the communities in which gene
    # is present. Tests the hypothesis that the gene is included in regions of higher degree
    # than all other genes in the panplasmidome

    # Filter for annotations appearing >= 5 times in the pangenome
    n_instances = annotations_df['Count'].groupby(
        'Annotation'
    ).sum().sort_values(ascending=False)
    n_instances = n_instances[n_instances>4]
    filtered_annotations = n_instances.index.to_numpy()
    filtered_annotations = filtered_annotations[filtered_annotations != '']
    filtered_annotations = filtered_annotations[filtered_annotations != 'hypothetical protein']
    # Get the degree and count for each selected gene
    filtered_annotations = annotations_df.loc[filtered_annotations]
    all_prompts = []
    for gene in np.unique(filtered_annotations.index.get_level_values('Annotation')):
        # List of mean degree of the community in which the gene appears
        mean_gene_degs = []
        for x in filtered_annotations.loc[gene].values:
            mean_gene_degs += [x[1]]*int(x[0])
        # List of mean degrees of communities in which all other genes appear
        all_other_mean_degs = []
        for x in annotations_df.drop(gene).values:
            all_other_mean_degs += [x[1]]*int(x[0])
        u_gene, pval = mannwhitneyu(mean_gene_degs, all_other_mean_degs, method='asymptotic')
        # Look for 95% confidence
        if pval <= 0.05:
            prompt = f'The distribution for average degree in the community for {gene} is different than that of other genes, with AUC of {u_gene*100/(len(all_other_mean_degs)*len(mean_gene_degs)):.3}% (p-value: {pval:.5}, Mann-Whitney U Test).\n'
            print(prompt)
            all_prompts += [prompt]
        
    # What are the nodes with higher degree?

    degrees = pd.DataFrame.from_dict(
        {x[0]: x[1] for x in nx.degree(m.get_motif_graph_as_nx())},
        orient = 'index',
        columns = ['Degree']
    )
    pf = get_protein_function(
        degrees.index.to_list(),
        phylo_config['paths']['protein-names']
    )
    # Remove duplicate entries and delete inverted commas and trailing accession from
    # the gene annotations
    pf = pf.loc[~pf.index.duplicated(keep='first')].apply(lambda x: ' '.join(x[1:-1].split()[:-1]))
    degrees['Annotation'] = pf.loc[degrees.index].values
    degrees.drop(
        degrees[degrees['Annotation']==''].index,
        inplace = True
    )
    degrees.drop(
        degrees[degrees['Annotation']=='hypothetical protein'].index,
        inplace = True
    )
    degrees = degrees.sort_values('Degree', ascending = False)

    # Keep the top 10 and plot
    top_degrees = degrees.iloc[:10]
    with plt.style.context('ggplot'):
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        fig, ax = plt.subplots(
            1,1,
            figsize = (10,7),
            constrained_layout = True
        )
        bar_ax = ax

        bar_ax.bar(
            np.arange(len(top_degrees['Annotation'].values)),
            top_degrees['Degree'].values,
            color = colors[5]
        )
        
        bar_ax.set_xlim(
            bar_ax.get_xlim()[0],
            bar_ax.get_xlim()[1]+1
        )
        avg_color = colors[1]
        
        bar_ax.plot(
            bar_ax.get_xlim(),
            [degrees['Degree'].mean()]*2,
            color = avg_color,
            lw = 3,
            ls = 'dashed'
        )
        bar_ax.text(
            bar_ax.get_xlim()[1]-.1, degrees['Degree'].mean()+.3,
            'Average',
            ha = 'right',
            fontweight = 'bold',
            fontsize = 13,
            color = avg_color
        )
        
        bar_ax.set_ylabel('Degree')
        n_words = 2
        bar_ax.set_xticks(
            np.arange(len(top_degrees))
        )
        bar_ax.set_xticklabels(
            [
                linebreak(x, max_char = 30)
                for x in top_degrees['Annotation']
            ],
            rotation = -45,
            ha = 'left',
            fontsize = 12
        )
        
        bar_ax.set_yscale('log')
        bar_ax.grid(False, axis='x')

        plt.show()