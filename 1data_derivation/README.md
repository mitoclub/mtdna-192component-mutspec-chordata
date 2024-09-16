# Dataset preparation details

[NeMu pipeline](nemu-pipeline.com) was used to derive mutational spectra

## Content

### Data generation

- [Jupyter notebook](./prepare_beasts_datasets.ipynb) used to generate input protein sequences for  pipeline execution
- [Used version of NeMu pipeline](./nemu.nf) for spectra calculation using protein sequences and nucleotide blast database (MIDORI2) as input
- [Config for the pipeline](./nemu_chordata.config)
- [Bash script](./run_chordata.sh) for parallel pipeline execution on around 10,000 mithochondrial genes of many different species

### Data aggregatrion and analysis

- [collect_spectra.ipynb](./collect_spectra.ipynb) - Jupyter notebook used to aggregate NeMu-pipeline outputs to simple tables with mutational spectra for the gene-species

### Final dataset

- [Mutational spectra (192-comp)](./dataset/MutSpecVertebrates192.csv.gz)
- [Input sequences info](./dataset/info.csv)