# Make detailed plasmid networks for the specified clusters
python build_networks.py --clusters 0 5 10 16 37
# Find conserved AMR regions
python bulk_find_conserved_regions.py --min_dist 0.1 --min_len 5 --max_len 9 --min_n 3 --out data/Conserved\ AMR\ Regions
# Get stats on conserved regions
python conserved_region_stats.py --input data/Conserved\ AMR\ Regions
