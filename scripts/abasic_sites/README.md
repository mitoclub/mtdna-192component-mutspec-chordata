# Scripts and notebooks for analysis of abasic sites in *Mus musculus* mtDNA

Fig. 6b in our paper

## Scripts

- contexter.py - derive Abasic Sites counts for heavy chain and Light chains. Control region excluded
- compic.py - compare Heavy and Light chains
- AbasicSites.Rmd - EDA: T-test, Wilcoxon test, correlations and plots
- AbasicSites_seq_logo.ipynb - plot triplet logos for both chains

## Raw data

- [Table](../../data/abasic_sites/41467_2022_33594_MOESM26_ESM.csv) from Source data (41467_2022_33594_MOESM26_ESM.xlsx Figure 9a,b) in the paper [1]
- [mm10 mtDNA](../../data/abasic_sites/mm10_ChrM_oneline.fasta)

## Processed data

- [Table](../../data/abasic_sites/AbasicSitesMtDNAcontext_HLcompare.csv) with abasic sites normalized counts in heavy and light chains

## Reference

1. Cai, Y., Cao, H., Wang, F. et al. Complex genomic patterns of abasic sites in mammalian DNA revealed by a high-resolution SSiNGLe-AP method. Nat Commun 13, 5868 (2022). https://doi.org/10.1038/s41467-022-33594-1
