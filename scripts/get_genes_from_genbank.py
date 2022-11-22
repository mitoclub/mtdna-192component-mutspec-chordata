from typing import Set
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
import pandas as pd

# Defining 13 mtDNA genes
genenames: Set[str] = {'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3',
                       'ND4', 'ND4L', 'ND5', 'ND6'}
# File to write our fasta files
files_to_write_nc = {gene: open(f"../data/fasta/{gene}.fasta", "w") for gene in genenames}
taxes = pd.read_csv('../data/taxa_gb.csv', sep='\t')

# Open GenBank file
with open('../data/sequence.gb') as file:
    for rec in SeqIO.parse(file, 'genbank'):
        try:
            feature: SeqFeature = None
            for feature in rec.features:
                # Check if feature is gene
                if hasattr(feature, 'type') and 'gene' in feature.qualifiers and feature.type == 'gene':
                    cur_genename = feature.qualifiers['gene'][0]
                    # Take only mtDNA genes
                    if cur_genename not in genenames:
                        continue
                    gene: SeqRecord = feature.extract(rec)
                    seq = str(gene.seq)
                    header = '>' + rec.annotations['organism']  # organism name

                    # sq = SeqRecord(seq, name=header)
                    files_to_write_nc[cur_genename].write(header + '\n')
                    files_to_write_nc[cur_genename].write(seq + '\n')

        except:
            print(feature)

for fout in files_to_write_nc.values():
    fout.close()
