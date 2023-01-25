# plasmidHGT

Code developed for _Identifying Antimicrobial Resistance Gene Transfer Between Plasmids_. 

By Marco Teixeira, Stephanie Pillay, Aysun Urhan and Thomas Abeel,
Delft Bioinformatics Lab, Delft University of Technology, Delft, 2628 XE, Netherlands

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
GenBank files for these plasmids, with features, should be placed in data/Plasmid_Features. The subdirectory data/Concatenated_Plasmid_FASTA should contain the plasmid sequence FASTA files concatenated per group. In both directories, plasmids should be placed inside one or multiple folders, corresponding to groups (e.g. host species).

`pipeline_clustering.sh` implements an initial Jaccard-based clustering of the plasmids in data/Plasmid FASTA and creates all necessary input files for `pipeline_networks.py`.

`pipeline_networks.py` builds and analyzes detailed plasmid networks for some of the clusters found with `pipeline_clustering.sh`. Alternatively, the data used in the manuscript can be loaded.

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

