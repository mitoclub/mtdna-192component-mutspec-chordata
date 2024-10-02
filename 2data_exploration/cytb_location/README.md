# Supplementary Fig. S4. Distribution of distances between Cytb gene and oriL in chordate species

## Analysis workflow

### Build env for mitos execution

```bash
sudo apt install libcairo2-dev

pip3 install mitos==2.1.9 reportlab pycairo rlpycairo
pip3 install requirements.mitos.txt
```

### Try to extract Ori annotation from genbank

[1cytb_location_from_genbank.ipynb](./1cytb_location_from_genbank.ipynb)

Out spectra dataset have a little intersection with genbank annotation data,
therefore we annotate MT genomes using Mitos software

### Prepare genomes 

[2prepare_genome_for_mitos.py](./2prepare_genome_for_mitos.py) will split genbank file to separate mtDNA genomes in fasta format 

### Annotate mt genomes using Mitos

[3run_annotation.sh](./3run_annotation.sh)

### Analyse annotations and get relative distance of Cytb based on OriL position

[4cytb_location_from_mitos.ipynb](./4cytb_location_from_mitos.ipynb)

## Results

- [Supplementary Fig. S4](./distance_between_oriL_start_and_Cytb_start.pdf)
- [Mitos OriL annotation](./mitos_oriL_labels.txt)
- [Genbank Cytb and OriL annotatation](./cytb_and_oriL_location_genbank_annot.csv)
