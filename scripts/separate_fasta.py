from Bio import SeqIO
from typing import Set
import pandas as pd


taxa_df = pd.read_csv('../data/taxa_gb.csv', sep='\t')

taxa_names: Set[str] = {'Actinopterygii', 'Amphibia', 'Lepidosauria', 'Mammalia', 'Aves'}
files_to_write = {taxa: open(f'../data/CytB_Classes/cytb_{taxa}.fasta', 'w') for taxa in taxa_names}

with open('../data/fasta/CYTB.fasta', 'r') as fasta_input:
    for rec in SeqIO.parse(fasta_input, 'fasta'):
        name = rec.description.replace(' ', '_')
        if name in taxa_df['Species'].values:
            ind_taxa = taxa_df[taxa_df['Species'].isin([name])].index[0]
            class_in_taxa = taxa_df.iloc[ind_taxa, 0]
            header = '>' + name
            files_to_write[class_in_taxa].write(header + '\n')
            files_to_write[class_in_taxa].write(str(rec.seq) + '\n')


for fout in files_to_write.values():
    fout.close()
