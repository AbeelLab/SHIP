#!/bin/bash
test_dir=./test
cd $test_dir

tmp_dir="tmp"

# Run Prokka. Disable annotations as these are likely bogus without a specific plasmid CDS database
annotations_dir="prokka_annotations"
mkdir $annotations_dir

echo "Annotating plasmids in $annotations_dir."
for FILE in *.fa
    do
    filename=${FILE%.fa*}
    prokka --outdir "$annotations_dir/$filename" --force --plasmid $filename --noanno $FILE
done

# Run AMRFinder+
amrfinder_dir="amrfinder_results"
mkdir $amrfinder_dir
echo "Finding AMR genes with AMRFinder+."
cat *.fa > "$tmp_dir/concat_plasmids.fa"
amrfinder -n "$tmp_dir/concat_plasmids.fa" -o "$amrfinder_dir/amrfinder_out.txt" --log "$amrfinder_dir/amrfinder.log"

# Run CD-HIT
cdhit_dir="cdhit_results"
mkdir $cdhit_dir
echo "Clustering homologs with CD-HIT."
touch $tmp_dir"/concat_cds.faa"
for FILE in *.fa
    do
    filename=${FILE%.fa*}
    echo -e "$(cat $annotations_dir/$filename/*.faa)\n" >> $tmp_dir"/concat_cds.faa"
done
cdhit -i $tmp_dir"/concat_cds.faa" -c 0.9 -n 5 -o $cdhit_dir

# Run SHIP
cd ..
mkdir "$test_dir/SHIP_results"
echo "Finding HGT events with SHIP."
python ship.py -a "$test_dir/$annotations_dir" -o "$test_dir/SHIP_results" \
    -c "$test_dir/$cdhit_dir" -r "$test_dir/$amrfinder_dir" -m 0.1 -l 5 -L 9 -n 3

rmtree "$test_dir/$tmp_dir"

echo "Done!"