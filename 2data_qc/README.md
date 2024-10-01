# Quality Control and Plotting of Mutational Spectra

This part of the project focuses on checking the quality of the obtained mutational spectrum for vertebrates. It includes plotting mutational spectra for all vertebrate species for the CytB, CO1, CO3, and ND2 genes. 

The initial data used in this analysis can be found in the [192-comp mutational spectra data](../1data_derivation/dataset/MutSpecVertebrates192.csv.gz).

## Content

### Notebooks

- [Plotting mutaional spectra](./spectra_visualization/plot_mutspec.ipynb) is used for drawing mutational spectra for all vertebrate species for the genes CO1, CO3, ND2, and CytB, as well as for drawing 192 component mutational spectra for different vertebrata classes and cancers
- [Location in CytB](./cytb_location/cytb_location.ipynb) is used to check the location of the CytB gene in different species
- [Check selection effect](./estimate_selection_effect/check_selection_effect.ipynb) TODO Bogdan
- [Mutational Spectrum EDA](./explore_mutations/dataset_EDA.ipynb) is used to look at EDA of the mutational spectrum
- [Gene Code Bias](./gencode_bias.ipynb) TODO BOGDAN !!!

### Scripts

- [Plotting mutational spectrum in COSMIC format](./spectra_visualization/plot_fig1b.py) is used to plot two pieces of Figure 1b, which contain mutational spectra in COSMIC format for all vertebrate species


## Used Data

- [192-comp mutational spectra data](../1data_derivation/dataset/MutSpecVertebrates192.csv.gz) contains species-specific mutational spectra for all available vertebrate species, covering the four genes: CO1, CO3, ND2, and CytB
- [] TODO Bogdan