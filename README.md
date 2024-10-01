# Deciphering the Foundations of Mitochondrial Mutational Spectra: Replication-Driven and Damage-Induced Signatures Across Chordate Classes

In this project, we use 2,591 CytB sequences of chordata species from the [NeMu pipeline](https://nemu-pipeline.com). We derive and analyse a 192-component mutational spectra to identify the main factors driving mtDNA mutagenesis in chordata species.

## Repository details

We splitted our data collection and analyses steps to subdirectories:

- [Cancer data preparation](./0cancer/)
- [Mutational spectra derivation](./1data_derivation/)
- [Data exploration](./2data_exploration/)
- [Chordates classes comparison](./3compare_classes/)
- [Signatures assigment (COSMIC)](./4signatures/)
- [Spectrum asymmetry analysis](./5asymmetry/)
- [AP sites context analysis](./6abasic_sites/)

[Main figures for the article](./main_figures/)

See more details in READMEs placed in subfolders

## Used software

- Python 3.9+
- R 4.4.1
- [NeMu-pipeline](https://nemu-pipeline.com/)

## Environment for analyses

### Python packages

see [requirements.txt](./requirements.txt)

- PyMutSpec-0.0.8
- statsmodels~=0.14.0
- scipy~=1.10.1
- biopython~=1.81
- matplotlib~=3.5.3
- numpy~=1.22.4
- pandas~=1.4.4
- seaborn~=0.12.2
- scikit-learn~=1.2.2
- SigProfilerPlotting~=1.3.23
- SigProfilerAssignment~=0.1.6
- umap_learn~=0.5.6

### R packages

- dplyr_1.1.4
- cosmicsig_1.1.1
- mSigAct_3.0.1
- ICAMS_3.0.5
- Biostrings_2.72.1

## Article

Mitochondrial mutation spectrum in Chordates: damage versus replication signatures, causes, and dynamics
Dmitrii Iliushchenko, Bogdan Efimenko, Alina G. Mikhailova, Victor Shamanskiy, Murat K. Saparbaev, Ilya Mazunin, Dmitrii Knorre, Wolfram S. Kunz, Philipp Kapranov, Stepan Denisov, Jacques Fellay, Konstantin Khrapko, Konstantin Gunbin, Konstantin Popadin

**bioRxiv** 2023.12.08.570826; doi: https://doi.org/10.1101/2023.12.08.570826 

