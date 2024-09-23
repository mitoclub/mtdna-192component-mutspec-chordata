# Decomposition of mtDNA spectra to COSMIC signatures

1. Need to renormalize chordates classes spectra to human genome (see details in article methods) - [py script](./0count_human_triplets_freqs.py)
2. Run SigProfilerAssignment cosmic_fit function to run signatures assignment, aggregate results and plot figures - [ipy notebook](./1signatures_analysis_sigpro.ipynb)

We splitted 192-component spectra to parts based on transitions freqs: Low High Diff (see details in article methods). Besides we renormalized spectra to emulate human genome mutation counts

3. To run another tool mSigAct we prepared priors for signatures distribution in the data based on SigProfilerAssignment results - [py script](./2prepare_priors_for_mSigAct.py)
4. Run mSigAct MAPAssignActivity func to decompose spectra by another approach - [R script](./3mSigAct_analysis.R)
5. Aggregation and little EDA of mSigAct results - [ipy notebook](./4aggregate_mSigAct_outputs.ipynb)
6. Plot figures for mSigAct output - [py script](./5plot_output_for_mSigAct.py)


- [plotActivity.py](./plotActivity.py) - visualization of barplots in the SigProfiler way
- [utils.py](./utils.py) - some functions for data processing