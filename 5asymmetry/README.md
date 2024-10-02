# Analysis of asymmetry in mtDNA of vertebrate

## Workflow

All steps are performed in [Asymmetry analysis notebook](./asymmetry_mutspec.csv)

1. Transformation of a 192-component mutational spectrum of verterbrate into a 96-component spectrum with asymmetry calculation
2. Correlation of mtDNA asymmetry and T- and R- asymmetry\
3. Plotting mtDNA asymmetry and T- R- relationships (fig. 4a and 4c)
4. Transformation of a 192-component mutational spectrum of birds (warm), fishes (cold) and cancer [high](../0cancer/data/for_asymmetry/ms_high_tsss_all.csv) and [low](../0cancer/data/for_asymmetry/ms_low_tsss_all.csv) into a 96-component spectrum
5. Plotting boxplots with asymmetry ratio in Cold vs Warm and High TSSS vs Low TSSS (fig. 4b)


## Used Data

- [192-comp mutational spectra data](../1data_derivation/dataset/MutSpecVertebrates192.csv.gz)
- Mutational spectra of cancer data devided into [low](../0cancer/data/for_asymmetry/ms_low_tsss_all.csv) and [high](../0cancer/data/for_asymmetry/ms_high_tsss_all.csv) zones

### External data

- [Table of T- and R- nuclear asymmetry](./T-R_plot_allCpG.txt) - taken from [Seplyarskiy et al. 2019](10.1038/s41588-018-0285-7)

## Final dataset

- [Asymmetry Table of mtDNA with T- and R- nuclear asymmetry](./asymmetry_mutspec.csv) - Supplementray Table S7

## Figures 

- [Scatter Plot of mtDNA asymmetry and T- and R- asymmetry](./figures/AsymmetryTRM.pdf) (Fig. 4a)
- [Error bars of mtDNA Mutational Bases and T- and R- asymmetry](./figures/AsymmetryErrorBars.pdf) (Fig. 4c)
- Scatter plots of mtDNA Mutatioanal bases separately and [T-](./figures/AsymmetryMutBaseT.pdf) and [R-](./figures/AsymmetryMutBaseR.pdf) asymmetry (Supl. Figures S13 and S14)
- [Boxplots of asymmetry ration in Cold vs Warm](./figures/ColdvsWarmBox.pdf) and [cancer LowTSSS and High TSSS](./figures/TBSSViolin.pdf) (Fig. 4b)

## Separate analysis

[Spectrum of POLG mutations](./pol_gamma/) -  symmetrical mutagenesis of mtDNA shaped by POLG (Supplementary Figure S15)