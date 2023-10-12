import os
import pandas as pd
import re
from Bio import SeqIO
from Bio.Data import CodonTable


def define_type_of_mutation(codon_to_change: str, changed_codon: str, sbs_pos):
    """ Define type of substitution, where:
     0 - NonSynonymous mutation
     1 - Synonymous mutation
     2 - Mutation in 4Fold site
     -1 - Stop codon, deletion or N in codon
    Args:
        codon_to_change: Ancestral Codon
        changed_codon: Descendant codon
        sbs_pos: position of substitution
    Returns:
        Type of substitution
    """
    # Define 4Fold sites
    syn_fourfold_codons = ['CTA', 'GTA', 'TCA', 'CCA',
                           'ACA', 'GCA', 'CGA', 'GGA',
                           'CTT', 'GTT', 'TCT', 'CCT',
                           'ACT', 'GCT', 'CGT', 'GGT',
                           'CTG', 'GTG', 'TCG', 'CCG',
                           'ACG', 'GCG', 'CGG', 'GGG',
                           'CTC', 'GTC', 'TCC', 'CCC',
                           'ACC', 'GCC', 'CGC', 'GGC']

    mito_table = CodonTable.unambiguous_dna_by_name['Vertebrate Mitochondrial']
    aa1 = mito_table.forward_table.get(codon_to_change, '*')
    aa2 = mito_table.forward_table.get(changed_codon, '*')
    if aa1 == '*' or aa2 == '*':
        return -1
    elif aa1 != aa2:
        return 0
    elif codon_to_change in syn_fourfold_codons and sbs_pos == 3:
        return 2
    elif aa1 == aa2:
        return 1


path = '../data/fasta'
fasta_files = os.listdir(path)
obs_mut = pd.read_csv('../data/ObsMutSpec.csv')
obs_sps = list(obs_mut['Species'])

with open('../data/exp_mut_spec.csv', 'w') as wr_file:
    col = ['Species', 'Gene', 'Mut3', 'MutType', 'Pos', '3Pos', 'Codon']  # define colnames of future csv
    wr_file.write('\t'.join(col))
    wr_file.write('\n')
    for fasta_index in fasta_files:
        fasta_path = os.path.join(path, fasta_index)
        # Go through all fasta records in gene file
        for record in SeqIO.parse(fasta_path, "fasta"):
            species_name = record.description
            species_name = species_name.replace(' ', '_')
            if species_name in obs_sps:
                seq = list(record.seq)
                gene_in_file = fasta_index.replace('.fasta', '')
                # Get all expected mutations
                for nucleotide_position in range(1, len(seq) - 1):
                    # Define position in codon previous and next nucleotides
                    # also define structure of codon
                    position_in_codon = nucleotide_position % 3
                    prev_nuc = nucleotide_position - 1
                    next_nuc = nucleotide_position + 1
                    start_of_codon = nucleotide_position - position_in_codon
                    codon = ''.join((seq[start_of_codon:start_of_codon + 3]))
                    subs = ['A', 'T', 'G', 'C']
                    all_nuc = ''.join(seq[prev_nuc] + seq[nucleotide_position] + seq[next_nuc])
                    if len(re.findall(r'[ATGC]', all_nuc)) == 3 and len(re.findall(r'[ATGC]', codon)) == 3:  # all three nucleotides must be A,T,G or C
                        subs.remove(seq[nucleotide_position])  # take only non-repeated nucleotides
                        for nucleotide_change in subs:
                            changed_codon = list(codon)
                            changed_codon[position_in_codon] = nucleotide_change
                            changed_codon = ''.join(changed_codon)
                            type_of_subs = define_type_of_mutation(codon, changed_codon, position_in_codon+1)
                            if position_in_codon == 2:
                                pos3 = 1
                            else:
                                pos3 = 0
                            # MutSpec in format A[T>G]C
                            complete_exp_codon = seq[prev_nuc] + '[' + seq[nucleotide_position] + '>' + \
                                                    nucleotide_change + ']' + seq[next_nuc]
                            out = [species_name, gene_in_file, complete_exp_codon, str(type_of_subs), 
                                    str(position_in_codon), str(pos3), str(codon)]
                            out = '\t'.join(out)
                            wr_file.write(out)
                            wr_file.write('\n')

        print(f'Gene file {fasta_index} completed')
