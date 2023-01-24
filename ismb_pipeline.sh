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
