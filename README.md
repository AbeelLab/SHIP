# plasmidHGT

Code developed for _Identifying Antimicrobial Resistance Gene Transfer Between Plasmids_. 

## Dependencies

Requires:
- AMRFinder+ 3.10.45
- MOBSuite 3.0.3
- Biopython 
- Prokka version 1.12
- CD-HIT version 4.8.1
- Pyvis version 0.2.1
- NetworkX version 2.8.6

# Running

The FASTA files for the plasmids to analyze should be in data/Plasmid FASTA. The
GenBank files for these plasmids, with features, should be placed in data/Plasmid_Features. In both directories, plasmids should be placed inside one or multiple folders, corresponding to groups (e.g. host species).

`pipeline_clustering.sh` implements an initial Jaccard-based clustering of the plasmids in data/Plasmid FASTA and creates all necessary input files for `pipeline_networks.py`.

`pipeline_networks.py` builds and analyzes detailed plasmid networks for some of the clusters found with `pipeline_clustering.sh`. Alternatively, the data used in the manuscript can be loaded.

To run each script, `bash pipeline_networks.py` or `bash pipeline_clustering.sh`

# Scripts

## `annotate.py`

Requires Prokka version 1.12

Runs Prokka on plasmid FASTA files to transfer the original annotations. Outputs the annotations to an Annotations folder. Requires a FASTA file with all products concatenated per group (FASTA folder) in data/Plasmid_Features.

If run for the first time, --build should be set.

### Options

- `--build`: Builds a protein FASTA database with all plasmid features. Must be set at least once.


