#!/bin/bash
tmp_dir="tmp"
mkdir $tmp_dir

# Run Prokka. Disable annotations or use a curated database of relevant plasmid products
# as Prokka annotations are likely bogus with the default databases
annotations_dir="prokka_annotations"
mkdir $annotations_dir

#echo "Annotating plasmids in $annotations_dir."
#for FILE in *.fa
#    do
#    filename=${FILE%.fa*}
#    prokka --outdir "$annotations_dir/$filename" --proteins "plasmid_blast_db/proteins.faa" $FILE
#done

# Run AMRFinder+
amrfinder_dir="amrfinder_results"
#mkdir $amrfinder_dir
#echo "Finding AMR genes with AMRFinder+."
#cat *.fa > "$tmp_dir/concat_plasmids.fa"
#amrfinder -n "$tmp_dir/concat_plasmids.fa" -o "$amrfinder_dir/amrfinder_out.txt" --log "$amrfinder_dir/amrfinder.log"

# Run CD-HIT
cdhit_dir="cdhit_results"
mkdir $cdhit_dir
#echo "Clustering homologs with CD-HIT."
#touch $tmp_dir"/concat_cds.faa"
#for FILE in *.fa
#    do
#    filename=${FILE%.fa*}
#    echo "$(cat $annotations_dir/$filename/*.faa)\n" >> $tmp_dir"/concat_cds.faa"
#done
#cdhit -i $tmp_dir"/concat_cds.faa" -c 0.9 -n 5 -o $cdhit_dir"/cdhit_output"

# Run SHIP
cd ..
mkdir "SHIP_results"
echo "Finding HGT events with SHIP."
ship -a "$annotations_dir" -o "SHIP_results" \
    -c "$cdhit_dir/cdhit_output.clstr" -r "$amrfinder_dir" -m 0.1 -l 5 -L 9 -n 3 -i

rmtree "$tmp_dir"

echo "Done!"