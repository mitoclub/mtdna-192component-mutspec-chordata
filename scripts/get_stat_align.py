from Bio import SeqIO
from typing import List

taxes_names: List[str] = ['Actinopterygii', 'Amphibia', 'Aves', 'Lepidosauria', 'Mammalia']

files_to_read = {taxa: open(f'../data/align_cytb/cytb_{taxa}.fna', 'r') for taxa in taxes_names}

with open('../data/statistics_algn.csv', 'w') as stat_file:
    stat_file.write('\t'.join(['Species', 'Class', 'PosOfErr', 'LenOfSeq']))
    stat_file.write('\n')
    indx = list(range(0, 5))
    for file, tax_id in zip(files_to_read.values(), indx):
        for rec in SeqIO.parse(file, 'fasta'):
            positions = [pos for pos, char in enumerate(rec.seq) if char == '!']
            if len(positions) == 0:
                pos_cor = '0'
            else:
                pos_cor = ';'.join(str(e) for e in positions)
            cls = file
            len_of_seq = str(len(rec.seq))
            out_line = rec.name + '\t' + taxes_names[tax_id] + '\t' + pos_cor + '\t' + len_of_seq
            stat_file.write(out_line)
            stat_file.write('\n')
