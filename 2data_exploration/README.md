# Exploration of Mutational Spectra Dataset

This part of the project focuses on **checking the quality** and **exploration** of the obtained mutational spectrum for vertebrates.

The initial data used in this analysis can be found in the [192-comp mutational spectra data](../1data_derivation/dataset/MutSpecVertebrates192.csv.gz).

## Content

- [Plotting mutaional spectra](./spectra_visualization/) - drawing the mutational spectra for all vertebrate species for the genes CO1, CO3, ND2, and CytB, as well as for drawing spectra for 5 vertebrata classes separately
- [Mutational Spectrum exploration](./explore_mutations/) - explorary data analysis (EDA) of the mutational spectrum dataset
- [Location of CytB gene in mtDNA](./cytb_location/) - checking the CytB gene location in available chordates species to prove the similar mutagenic condition of Cytb during the replication (different genes stay different time in single-strand condition; is this feature common among chordates species mtDNA?)
- [Comparison of species-specific spectra](./genes_cls_exploration/) - compare Cytb mutaitonal spectra of chordates species and compare their expected mutations counts
- [Check spectrum bias based on genetic code](./syn_spectrum_bias_based_on_gencode.ipynb) - some mutations contexts cannot be observed in synonymous mutational spectrum. We discover this gencode feature in this notebook.
- [Estimation of synonymous sites neutrality](./estimate_selection_effect/) - To assess the neutrality of synonymous sites in mtDNA, we compared expected substitutions counts (192-component vectors) in different sites samples according to estimated site variability rate for each species.

## Main Files 

- [Class-specific spectra](./spectra_visualization/ClassesSpectraCytb.csv) - Supplementary Table S5, class specific mtDNA mutational spectra for each class and gene 
- [Average expected substitutions counts](./average_class_exp_freqs.csv) - Supplementary Table S4
## Main Figures

- [mtDNA 12-comp Mutational Spectrum for CytB](./spectra_visualization/MutSpec12VertCytb.pdf) - 12-component mutational spectra fig 1a.
- mtDNA Mutational Spectrum for CytB in COSMIC format [upper part](./spectra_visualization/SBS_96_plots_vert_1.pdf) and [bottom part](./spectra_visualization/SBS_96_plots_vert_2.pdf) - Two parts of mtDNA mutational spectrum in COSMIC format for figure 1b 
