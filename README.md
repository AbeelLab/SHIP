# plasmidHGT

Code developed for _Identifying Antimicrobial Resistance Gene Transfer Between Plasmids_. 

By Marco Teixeira, Stephanie Pillay, Aysun Urhan and Thomas Abeel,
Delft Bioinformatics Lab, Delft University of Technology, Delft, 2628 XE, Netherlands

## Dependencies

Requires:
- Python 3.8
- AMRFinder+ version 3.10.45
- MOBSuite version 3.0.3
- Biopython version 1.79
- Prokka version 1.12
- CD-HIT version 4.8.1
- Pyvis version 0.2.1
- NetworkX version 2.8.6
- markov-clustering 0.0.6
- pyMC 4.2.2
- Numpy
- Pandas
- Joblib
- Matplotlib
- Pyyaml
- Scipy
- scikit-learn

Developed for Ubuntu 20.04.

# Running

The FASTA files for the plasmids to analyze should be in data/Plasmid FASTA. The
GenBank files for these plasmids, with features, should be placed in data/Plasmid_Features. The subdirectory data/Concatenated_Plasmid_FASTA should contain the plasmid sequence FASTA files concatenated per group. In both directories, plasmids should be placed inside one or multiple folders, corresponding to groups (e.g. host species).

`pipeline_clustering.sh` implements an initial Jaccard-based clustering of the plasmids in data/Plasmid FASTA and creates all necessary input files for `pipeline_networks.py`.

`pipeline_networks.py` builds and analyzes detailed plasmid networks for some of the clusters found with `pipeline_clustering.sh`. Alternatively, the data used in the manuscript can be loaded.

`pipeline_mobsuites.py` runs MOBSuite on the plasmid sequences. Required for further network analysis.

To run each script, `bash pipeline_networks.py` or `bash pipeline_clustering.sh`

# Scripts

## `annotate`

Requires Prokka version 1.12

Runs Prokka on plasmid FASTA files to transfer the original annotations. Outputs the annotations to an Annotations folder. Requires a FASTA file with all products concatenated per group (FASTA folder) in data/Plasmid_Features.

If run for the first time, --build should be set.

### Options

- `--build`: Builds a protein FASTA database with all plasmid features. Must be set at least once.

## `make_plasmid_info_df`

Writes plasmid information into a TSV file. This information includes host species, strain and BioSample ID.

## `find_amr`

Requires AMRFinder+ 3.10.45

Calls AMRFinder+ to find AMR genes in all plasmids. Uses the concatenated FASTA files in data/Concatenated_Plasmid_FASTA. AMRFinder+ sometimes fails when using individual FASTA files. Creates one AMRFinder file per group.

## `make_amrfinder_df`

Processes the output files from AMRFinder+ and concatenates it into a single TSV file.

## `homolog_clustering`

Requires CD-HIT version 4.8.1

Defines homolog gene clusters with 90% amino acid similarity using CD-HIT.

## `replace_with_homologs`

Processes the CD-HIT output to replace the CDS annotations in the plasmid sequences with the representative CDS.

## `cluster_plasmids`

Clusters plasmids according to the Jaccard similarity on gene content. Plots the panplasmidome of the plasmids in each cluster. The computed clusters are saved in data/Jaccard Clusters.

### Options

- `--load`: Loads pre-computed clusters and plots stats and the panplasmidome graphs.

## `cluster_analysis`

Provides statistics for the clusters computed by `cluster_plasmids.py`. Can use either the clusters more recently obtained with `cluster_plasmids.py` or those used in the paper.

### Options

- `--paper`: If set, uses the clusters computed and described in the paper, instead of those lastly obtained with `cluster_plasmids.py`.

## `prepare_mobsuite_input`

Copies the plasmid sequence FASTA files into one folder, to act as input for MOBSuite. Uses the files in data/Plasmid Networks or the clusters defined in the paper. Includes all plasmids in the Jaccard-Subclusters.tsv file.

*Note:* if `--paper` is not set, `build_networks.py` must be executed before calling `prepare_mobsuite_input`.

### Options

- `--paper`: If set, uses the clusters computed and described in the paper.

## `make_mobsuite_db`

Concatenates the MOBSuite outputs and creates one table per plasmid typing scheme.

## `build_networks`

Creates detailed plasmid similarity networks based on the proposed method for the specified clusters. Stores a SimpleDendrogram object that can be used for further analysis. Can use the Jaccard-based clusters used in the paper.

### Options

- `--paper`: If set, uses the Jaccard-based clusters computed and described in the paper.
- `--clusters`: One or more clusters for which to build the more detailed networks. Must be integers corresponding to Jaccard-based cluster numbers. Redundant and not required if `--paper` is set.

### Example

`python build_networks.py --clusters 0 5 10 16 37`

## `bulk_find_conserved_regions`

Finds conserved regions containing AMR genes on the clusters and networks built using `build_networks.py`. Alternatively, it can use the networks used in the paper. Requires setting search criteria. Creates a table with all found regions.

### Options

- `--paper`: If set, uses the networks computed in the paper.
- `--min_dist`: Minimum average plasmid distance for region inclusion. Default is 0.1.
- `--min_len`: Minimum fragment length in CDS. Default is 5.
- `--max_len`: Minimum fragment length in CDS. Default is 9.
- `--min_n`: Minimum number of plasmids containing a region for inclusion. Default is 3.
- `--out`: Output directory. Default is CWD.

### Example

`python bulk_find_conserved_regions.py --min_dist 0.1 --min_len 5 --max_len 9 --min_n 3 --out data/Conserved\ AMR\ Regions`

## `conserved_region_stats`

Gets statistics of the conserved AMR regions found by bulk_find_conserved_regions.py.

### Options

- `--input`: Directory containing the output files from bulk_find_conserved_regions.py.

## `find_conserved_regions`

Allows manual search of conserved regions containing AMR genes. Requires a pre-computed phylogeny.

### Options

- `--paper`: If set, uses the networks computed in the paper.
- `--cluster`: Jaccard cluster in which to find AMR genes. If not specified, assumes that there is only one cluster.

## `motif_analysis`

Extracts the regions described in the paper (recombination in E. faecalis and complex class 1 integron in _E. coli_ and _K. pneumoniae_) into FASTA files. Performs alignment with Biopython's pairwise2 module and prepares files for BLAST alignment.

### Options

- `--region`: Conserved region to analyse. Either "integron" (default) or "recombination".

## `pangenome_analysis`

Extracts stats about the panplasmidome, including local variability. Uses the networks in data/Plasmid Networks or those used in the paper.

### Options

- `--paper`: If set, uses the networks computed in the paper.
- `--cluster`: Jaccard cluster to analyze. If not specified, assumes that there is only one cluster.

## `plot_networks`

Plots  the detailed plasmid similarity networks obtained with `build_networks.py`. Shows the estimated frequencies for evolutionary events. Optionally, shows networks colored on plasmid types. MOBSuite must be called first (using `pipeline_mobsuite.sh`) if `--types` is set.

### Options

- `--paper`: If set, uses the networks computed in the paper.
- `--panplasmidomes`: If set, shows the panplasmidome graphs for each subcluster in the networks.
- `--types`: If set, plots plasmid networks colored on plasmid types from MOBSuite.

## `length_thresholds`

Uses KDE and MLE to estimate the thresholds for recombination/mutation/IS synteny blocks.
Length thresholds are defined using maximum likelihood estimation. The probability distributions for synteny block lengths are approximated using Gaussian kernel density
estimation. The length of mutation regions is modeled using the size of all CDS in the dataset; for integrons/transposons, uses the length of all ESKAPE organisms entries in ISFinder as of Nov 2022; the sizes of recombination blocks are modeled by sampling seven CDS with replacement and uniform probability from the dataset, and adding their lengths. As a result, extensive mutation regions were defined as having less than 915 bp and connected to
the same node as another gene present only in the other plasmid. Regions shorter than 915 bp that did not fulfill other requirements as outlined in the paper, and those shorter than 2752 bp were considered integrons/transposons.

Not required to run other scripts.
