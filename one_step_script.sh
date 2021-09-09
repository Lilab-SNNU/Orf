#!user/bin/bash


python make_index_genomes.py -y config.yaml
python filter_coding_reads.py -y config.yaml
python filter_coding_circ.py -y config.yaml

python filter_coding_orf.py -y config.yaml
python visual_circ_orf.py -y config.yaml
