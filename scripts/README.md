# Processing and Analysis

**[AP sites analysis](./abasic_sites/) - (fig. 6b)**

## Notebooks

- [asymmetry_TR.ipynb](asymmetry_TR.ipynb) - counts mitochondrial asymmetry and correlation with T and R assymetry, also compare mitochondrial asymmetry of cold (fishes) and warm (aves) species ; creates **AsymmetryTRM.pdf, AsymmetryErrorBars.pdf, TBSS_and_Temperature.pdf** (fig. 7a, fig. 7b, fig. 7c)
- [cancer_analysis.ipynb](cancer_analysis.ipynb) - prepare data for TSSS asymmetry analysis
- [class_homogenity.ipynb](class_homogenity.ipynb) - UMAP, clustermaps for species spectra in classes (fig. S1)
- [count_codons.ipynb](count_codons.ipynb) - count how many mutations with context in each species for observed and expected mutations in cytB and merge it; creates **counted_codons_cytb.csv**
- [get_mutspec.ipynb](get_mutspec.ipynb) - spectra derivation and vizualization (fig. 2, fig. 5ab, fig. S2, fig. S4, fig. S5) ; creates file with 192 component mutational spectrum for each of 974 species (SupplementaryTable1.csv)
- [observed_mutspec.ipynb](observed_mutspec.ipynb) - gives 192 component format for observed mutations and annotation; creates **ObsMutSpec.csv**
- [signatures_analysis.ipynb](signatures_analysis.ipynb) - signatures analysis (fig. 6a)
- [similarity_jack.ipynb](similarity_jack.ipynb) - performs comparison using cosine similarity between classes and then mutbases using jackknife method; creates **fig3, fig4**
- [TsTv_and_little_asymmetry.ipynb](TsTv_and_little_asymmetry.ipynb) - Ts and Tv precence distribution (fig. S3)
- [PolGspectrum.Rmd](PolGspectrum.Rmd) & [PolGspectrum.html](PolGspectrum.html) -

## Python functions

- [utils.py](utils.py) - helper code used in several analyses
- [plotActivity.py](plotActivity.py) - code for visuzlization of signatures decomp results

## Scripts

- [EDA_mutspec.R](EDA_mutspec.R) - creates **EDA_mutspec.pdf** with description of observed mutations and species destribution
- [get_genes_from_genbank.py](get_genes_from_genbank.py) - takes mitochondrial sequences of all genes from genbank; creates list of files in **../data/fasta/{gene}.fasta**
- [from_gb_to_exp_subs.py](from_gb_to_exp_subs.py) - performs *in silico* saturation mutagenesis procedure on RefSeq fasta files ; creates **exp_mut_spec.csv** with expected substiotutions and context for each species
- [get_taxa_from_gb.py](get_taxa_from_gb.py) - takes taxonomy from genbank file; creates **taxa_gb.csv**

