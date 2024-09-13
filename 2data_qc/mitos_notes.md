## Build env for mitos execution

```bash
sudo apt install libcairo2-dev

pip3 install mitos reportlab pycairo rlpycairo
```

## Annotate mt genomes

First of all run analysis and data preparation in `scripts/cytb_location.ipynb`

Next:

```bash
#single run
#runmitos.py --input data/sample_genome.fasta -c 2 -o mitos_test -R data/mitos_db/ -r refseq89m --noplots --prot 0 --trna 0 --rrna 0 --intron 0

# run prepared genome dataset on all CPUs
parallel  mkdir data/mitos_output/{/.} ';' runmitos.py --input {} -c 2 -o data/mitos_output/{/.} -R data/mitos_db/ -r refseq89m --noplots --prot 0 --trna 0 --rrna 0 --intron 0 ::: data/genomes_to_annotate/*.fasta

# aggregate information
cat data/mitos_output/*/result.mitos | grep OL > data/mitos_oriL_labels.txt
```

Return to analyses in `scripts/cytb_location.ipynb`
