from Bio import SeqIO

needed_taxa = ['Aves', 'Mammalia', 'Amphibia', 'Actinopterygii', 'Lepidosauria']
with open('../data/taxa_gb.csv', 'w') as wr_file:
    wr_file.write('\t'.join(['Class', 'Species']))
    wr_file.write('\n')
    with open('../data/sequence.gb') as file:
        for rec in SeqIO.parse(file, 'genbank'):
            taxes = rec.annotations.get('taxonomy')
            needed_class = None
            for taxa in taxes:
                if taxa in needed_taxa:
                    needed_class = taxa
            if needed_class is None:
                continue
            org_name = rec.annotations.get('organism').replace(' ', '_')
            wr_file.write('\t'.join([needed_class, org_name]))
            wr_file.write('\n')