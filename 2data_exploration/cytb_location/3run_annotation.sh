# run incomplete annotation (only Ori) on prepared genomes dataset on all CPUs
parallel  mkdir mitos_output/{/.} ';' runmitos.py --input {} -c 2 -o mitos_output/{/.} -R mitos_db/ -r refseq89m --noplots --prot 0 --trna 0 --rrna 0 --intron 0 ::: genomes_to_annotate/*.fasta

# aggregate separate annotations to single table
cat mitos_output/*/result.mitos | grep OL > mitos_oriL_labels.txt


#single run for tests
#runmitos.py --input data/sample_genome.fasta -c 2 -o mitos_test -R data/mitos_db/ -r refseq89m --noplots --prot 0 --trna 0 --rrna 0 --intron 0
