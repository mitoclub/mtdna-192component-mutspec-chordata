from collections import Counter

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

## Prepare fasta genomes for custom OriL annotation using Mitos

df_loc = pd.read_csv('./cytb_and_oriL_location_genbank_annot.csv')

acc_to_annotate = set(df_loc['acc'].unique())
len(acc_to_annotate)
path_to_vertebrates_mtdna = './vertebrates_mtdna.gb'
for rec in SeqIO.parse(path_to_vertebrates_mtdna, 'genbank'):
    if rec.id in acc_to_annotate:
        organism = rec.annotations.get('organism').replace(' ', '_')
        if not organism:
            continue
        try:
            SeqIO.write(rec, f'./genomes_to_annotate/{organism}.fasta', 'fasta')
        except:
            pass