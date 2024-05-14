## How to build env for mitos execution

```bash
sudo apt install libcairo2-dev

pip3 install mitos reportlab pycairo rlpycairo

runmitos.py --input data/sample_genome.fasta -c 2 -o mitos_test -R data/mitos_db/ -r refseq89m --noplots --prot 0 --trna 0 --rrna 0 --intron 0

bash ./scripts/annotate_mt_genomes_with_mitos.sh

cat data/mitos_output/*/result.mitos | grep OL > data/mitos_oriL_labels.txt
```
