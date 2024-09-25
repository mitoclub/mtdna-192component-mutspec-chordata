# Reanalyse human cancer mtDNA mutations

## Raw data

1. [raw mutations](./data/mtDNA_snv_Oct2016.txt)
2. [reference human mtDNA](./data/NC_012920.1.gb)

## Scripts

1. [1cancer_analysis.ipynb](./1cancer_analysis.ipynb) - analyse mutations and calculate spectra
2. [2calculate_asymmetry.py](./2calculate_asymmetry.py) - calculate mtDNA asymmetry based on different genome regions

## Figures

- [Supplementary Fig. S8](./figures/syn_spectrum.pdf)
- [Supplementary Fig. S11c](./figures/cancer_ChTh.pdf)
- [Supplementary Fig. S11d](./figures/cancer_AhGh.pdf)
- [Supplementary Fig. S20a](./figures/patients_mut_num.pdf)
- [Supplementary Fig. S20b](./figures/cancer_samples_nmut.png)
- [Supplementary Fig. S20c](./figures/human_cancer_spectra_syn_samples_umap.png)
- [Supplementary Fig. S21b](./figures/nmuts_classes_cancer.png)

## Calculated mutational spectra

Mutaional spectra of human cancer mtDNA is availavle [here](./data/cancer_mutspec.csv)

We used it to compare with chordata classes spectra.

We calculated spectra based on different mutational samples, but used primarily only SYN spectra:

- syn (based on synonymous mutations)
- syn-fourfold (based on synonymous mutations at fourfold sites)
- all (syn+nonsyn)
- D-loop (mutational spectra of the control region)
- ExDloop (syn+nonsyn spectra on the full genome excluding control region mutaions)

## References

1. Yuan, Y., Ju, Y.S., Kim, Y. et al. Comprehensive molecular characterization of mitochondrial genomes in human cancers. Nat Genet 52, 342–352 (2020). https://doi.org/10.1038/s41588-019-0557-x
2. Used Data - https://ibl.mdanderson.org/tcma/download.html
