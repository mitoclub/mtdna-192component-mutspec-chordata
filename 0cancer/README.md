# Reanalyse human cancer mtDNA mutations

## Raw data

1. [cancer mutations](./data/mtDNA_snv_Oct2016.txt)
2. [reference human mtDNA](./data/NC_012920.1.gb)

## Scripts

1. [1cancer_analysis.ipynb](./1cancer_analysis.ipynb) - analyse mutations and calculate spectra
2. [2calculate_tbss_specific_spectra.py](./2calculate_tbss_specific_spectra.py) - calculate mtDNA asymmetry based on different genome regions

## Figures

See `1cancer_analysis.ipynb`

- [Supplementary Fig. S8](./figures/syn_spectrum.pdf)
- Supplementary Fig. S11 [c](./figures/cancer_ChTh.pdf) [d](./figures/cancer_AhGh.pdf)
- Supplementary Fig. S19 [a](./figures/patients_mut_num.pdf) [b](./figures/cancer_samples_nmut.png) [c](./figures/human_cancer_spectra_syn_samples_umap.png)
- [Supplementary Fig. S20b](./figures/nmuts_classes_cancer.png)

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
2. Link to the cancer mutations - https://ibl.mdanderson.org/tcma/download.html
