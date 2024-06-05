# mtdna-192component-mutspec-chordata

## Deciphering the Foundations of Mitochondrial Mutational Spectra: Replication-Driven and Damage-Induced Signatures Across Chordate Classes

In this project, we use 2,591 chordata species with the CytB mitochondrial gene from the [NeMu pipeline](https://nemu-pipeline.com). We plot a 192-component mutational spectrum to identify the main factors driving mutagenesis in mitochondrial DNA in chordata species.

## Main authors

* Dmitrii Iliushchenko
* Bogdan Efimenko
* Konstantin Popadin

## Environment

- python 3.9+
- [PyMutSpec-0.0.8](https://pypi.org/project/PyMutSpec/)
- statsmodels~=0.14.0
- scipy~=1.10.1
- biopython~=1.81
- matplotlib~=3.5.3
- numpy~=1.22.4
- pandas~=1.4.4
- seaborn~=0.12.2
- scikit-learn~=1.2.2
- SigProfilerPlotting~=1.3.23
- umap_learn~=0.5.6

## Data

Here is the description of all relatable data that we use in paper:...

## Workflow

### 1. NeMu pipeline 


### 2. Collection of a mutation spectrum


### 3. Plot mutational spectrum for CytB, COX1, COX3, ND2

Notebook `notebooks/plot_mutspec.ipynb` plots 12- and 192- component mutational spectrum for 4 most offen genes (COX1, COX3, ND2, CytB). 
    
    Input:
        - data/new_dataset/MutSpecVertebrates192.csv.gz  - calculated mutation spectra for all vertebrate classes for COX1, COX3, ND2, CytB
        
    Output: 
        - pictures/Mut12Vert.pdf - 12-component mutaitonal spectrum for all available chordates in Cytb
        - pictures/MutSpec192Vert.pdf - 192-component mutaitonal spectrum for all available chordates in Cytb
        - pictures/fig5a.pdf - C>T substitution with all contexts
        - pictures/fig5b.pdf - A>G substitution with all contexts
        - ToPaper/SuplFiles/SupplementaryTable1.csv - Table in which each raw represents species with all possible 192 substitutions

        additional pictures of mutational spectra in COX1, COX3, ND2 

### 4. Comparison of the mutational spectrum

Notebook `notebooks/similarity_jack.ipynb` performs JackKnife Similarity Analysis. In this process we use cosine similarity to compare spectra. We estimated differences between classes using jackknife resampling of species spectra. Next, we randomly selected **20** species from each pair of classes, calculated the 192-component mutational spectrum for both classes, and computed the cosine similarity of either the overall mutational spectrum or its parts

    Output:
        -
        -
        -
        -

