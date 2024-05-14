#!/bin/bash

# run on all CPU threads
parallel --dry-run mkdir data/mitos_output/{/.} ';' runmitos.py --input {} -c 2 -o data/mitos_output/{/.} -R data/mitos_db/ -r refseq89m --noplots --prot 0 --trna 0 --rrna 0 --intron 0 ::: data/genomes_to_annotate/*.fasta

