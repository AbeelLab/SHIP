#!/bin/bash

# Create concatenated GBFF file for all plasmid features
# Annotate plasmids
python annotate.py --build
# Create DataFrame with plasmid host species, strain and NCBI accession information
python make_plasmid_info_df.py
# Find AMR genes with AMRFinder+
python find_amr.py
# Join outputs into a DataFrame
python make_amrfinder_df.py
# Cluster homolog proteins
python homolog_clustering.py
python replace_with_homologs.py
# Initial Jaccard-based clustering
python cluster_plasmids.py
# Get stats and plots from the clustering results
python cluster_analysis.py