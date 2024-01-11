# SHIP
## Synteny-aware HGT Identification in Plasmids

SHIP is a tool to find antimicrobial resistant (AMR) DNA regions with evidence of lateral transfer between plasmids. To do this,
SHIP starts by computing pairwise distances between all plasmids provided as input. These distances
reflect divergence times between the query plasmids, by considering the number of structural variants
(SV) needed to obtain one plasmid sequence from another. SHIP also infers the frequencies of SVs
from the input data. Then, SHIP searches for regions containing AMR genes present in plasmids with
a large average pairwise distance.

**Note:** currently, SHIP only looks for regions containing AMR genes and complete plasmid assemblies.

## Installation

SHIP is available in PyPI, and its latest version can be installed by running:

```pip install ship_plasmid```

## Dependencies

Installing SHIP through pip or conda should also install all required dependencies. However, SHIP
uses as input files created by other bioinformatics tools, namely:
- Prokka or Bakta
- CD-HIT
- AMRFinderPlus

## Usage

For basic usage, SHIP can be executed with the command

```ship -a path/to/annotations/ -c path/to/orthogroups/clusters.clstr -r path/to/amr.tsv -o path/to/output/```

### Inputs

SHIP requires gene annotations from Prokka or Bakta. However, product prediction is not required to run SHIP
(altough it may make results more dificult to analyze), so feel free to set the ```--noanno``` flag when
running Prokka/Bakta to speed computation. Furthermore, note that the standard Prokka databases are not
adequate for plasmid annotation, so consider creating your own if you need functional annotations.
The directory containing the output from Prokka/Bakta should be provided to SHIP using the ```--annotations``` or
```-a``` argument. SHIP assumes that the plasmid assemblies are complete, and it may fail if provided draft 
assemblies. SHIP will take all plasmid annotations inside this directory as the set of plasmids in which to
find AMR regions with evidence of HGT.

After predicting ORF in your plasmid sequences with Prokka/Bakta, you should cluster them into ortholog groups
with CD-HIT. SHIP uses SVs to estimate plasmid distances, and disregards SNVs. When developing SHIP, and in the
experiments outlined in its publication, a value of 90% amino acid similarity was used when finding ortholog groups.
The .clstr output file from CD-HIT should be provided to SHIP in the ```--cdhit``` or ```-c``` argument.

Finally, you should provide SHIP with the report from AMRFinderPlus through the ```--amr``` or ```-r``` argument.

#### Optional arguments

```--plot-dendogram```, ```-d```: If set, plots a dendrogram of the plasmid phylogeny as estimated by SHIP.

```--min_dist```, ```-m```:  Minimum average plasmid distance between plasmids containing a region for it to be considered as having evidence for HGT. Default is 0.1.

```--min_len```, ```-l```:  Minimum number of CDS in a fragment with evidence for HGT for it to be included. Default is 5.

```--max_len```, ```-L```:  Maximum searched number of CDS in fragments. Longer fragments with evidence for HGT will be split in fragments of ```--max_len```. Default is 9.

```--min_n```, ```-n```:  Minimum number of plasmids containing a region with evidence for HGT for it to be included in the output. Default is 3.

```--keep-intermediate```, ```-i```: If set, keeps all intermediate files in path/to/SHIP/output/tmp.

### Output files

SHIP outputs a report table (HGT_AMR_regions.tsv) in TSV format, in which each row represents a region with evidence of HGT
between plasmids. This table contains the representative ID of the ortholog groups present in the region, the IDs of the
plasmid sequences containing it, the average plasmid distance, and the transposases and integrases present, if gene annotations are provided.

SHIP also creates two pickle files - one holds the network of plasmid distances (plasmid_dissimilarity.pkl), while another is the RegionFinder object used by SHIP to search for the plasmid regions with evidence of HGT.

## Citation

Marco Teixeira, Stephanie Pillay, Aysun Urhan, Thomas Abeel, SHIP: identifying antimicrobial resistance gene transfer between plasmids, 
Bioinformatics, Volume 39, Issue 10, October 2023, btad612, https://doi.org/10.1093/bioinformatics/btad612