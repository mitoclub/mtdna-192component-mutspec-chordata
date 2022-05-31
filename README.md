# mtdna-192component-mutspec-chordata

## Authors:
- Dmitrii Iliushchenko
- Konstantin Popadin (**Supervisor**)

## Environment

- python 3.9+

## Workflow

### 1. Download sequences from RefSeq.

Using web-site [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) with options `((complete mitochondrion RefSeq)) AND "chordates"[porgn:__txid7711] `, download GenBank file. In our case, we got 'sequence.gb'

### 2. Get genes from GenBank

To obtain all 13 mitochondrial genes use script `scripts/get_genes_from_genbank.py`. It takes GenBank file and parses the data to 13 separated fasta files with all chordates included.

### 3. Obtain observed mutspec

Obtaining the 3rd and 5th component mutation spectrum is carried out using the `scripts/observed_mutspec.ipynb` script and `data/VertebratePolymorphisms.MutSpecDataWithNonSyn.txt` table. At the output, we get a table in `.csv` format, you can see an description of this table below.

* Species - Name of species

* MutType- Type of a substitution, where:

    a. `0` - Non Synonymous
    
    b. `1` - Synonymous
    
    c. `2` - FourFold Synonymous

* 3Pos - `1` if the substitution is in the third position of codon, else `0` (second or first position)
* Mut3 - 3 component format of substitution
* Mut5 - 5 component format of substitution
* Pos - position of substitution in codon

![](https://github.com/DIliushchenko/mtdna-192component-mutspec-chordata/blob/main/pictures/Example_Obs_mutspec.png?raw=true)
