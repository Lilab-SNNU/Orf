#!user/bin/bash


python filter_coding_reads_circ.py -y config.yaml
python filter_coding_orf.py -y config.yaml
python visual_circ_orf.py -y config.yaml
