# Dataset preparation details

[NeMu pipeline](https://github.com/mitoclub/nemu-webpage)

## Content

- [Used version of NeMu pipeline](./nemu.nf) for spectra calculation using protein sequences and nucleotide blast database (MIDORI2) as input
- [Config for pipeline](./nemu_chordata.config)
- [Jupyter notebook](./prepare_beasts_datasets.ipynb) used to generate input protein sequences for  pipeline execution
- [Bash script](./run_chordata.sh) for parallel pipeline execution on around 10,000 mithochondrial genes of many different species
- [Input sequences info](../../data/new_dataset/info.csv)