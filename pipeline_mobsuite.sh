#!/bin/bash


# Prepare FASTA files for MOBSuite
python prepare_mobsuite_input.py
# Run MOBSuite
for file in data/MOBSuite\ Input/*
do
    mob_typer -i "$file" -o "data/MOB\ Types/$file".txt
    echo $file
done
python make_mobsuite_df.py