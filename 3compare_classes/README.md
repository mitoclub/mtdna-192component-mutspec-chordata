# Comparison of the mutational spectra

## Workflow

All steps are performed in [Similarity analysis notebook](./similarity_jack.ipynb)

1. We conducted pairwise comparisons of mutational spectra at the species-specific level within and between classes, considering all possible combinations

We randomly take 50 species of each class in the pair, calculate MS for each class in pair and calculate cosine similarity in tyhis pair. We have done this process 1000 times for each unique pair. In cancer case, we take 50% of all pacient cancer data 50 instead of 50 species due to low number of mutations in each patient

2. Draw matrix plots of similarity analysis, where in each cell we have three numbers: Q1, Q2, Q3 from 1000 pairwise comparisons in combination 


## Used Data

- [192-comp mutational spectra data](../1data_derivation/dataset/MutSpecVertebrates192.csv.gz)
- [Mutational spectra based on synonymous mutations of cancer patient data](../0cancer/data/human_cancer_spectra_patient_specific_syn.csv)

## Output Data 

All possible pairwise comparisons using cosine similarities:

- [For all mutaitons](./data/CosSimAll.csv)
- [For transitions](./data/CosSimTS.csv)
- [For transversions](./data/CosSimTV.csv)
- For base substitutions separately in transitions: [A>G](./data/CosSimA%3EG.csv); [C>T](./data/CosSimC%3ET.csv); [G>A](./data/CosSimG%3EA.csv); [T>C](./data/CosSimT%3EC.csv)

## Figures 

Cosine similarity matrix plots for:

- [All mutations](./figures/Jackknife_SimilarityAll.pdf) - Supplementary Fig. S9a
- [Transitions](./figures/Jackknife_SimilarityTS.pdf) - Figure 2a
- [Transversions](./figures/Jackknife_SimilarityTV.pdf) - Supplementary Fig. S9b
- [Transitions separately](./figures/Jackknife_SimilarityTS_Separate.pdf) - Figure 2b
