## Plot output stacked barplots for mSigAct

import pandas as pd
from plotActivity import plotActivity

custom_colors = [
    '#acf2d0', '#63d69e', '#f8b6b3', '#c4abc4', '#f2aeae', '#d9f7b0', 
    '#faf1dc', '#fcebc2', '#fae4af', '#fae1a5', '#fcde97', '#fad682', 
    'tab:pink', 'tab:orange', 'tab:purple', 'tab:olive', 'tab:brown', 
    'tab:red', 'tab:green', 'tab:cyan', 'deeppink', 'lightgray', 'blueviolet', 
    'chocolate', 'darkgreen', 'dodgerblue', 'gray', 'salmon']



# total
inpath  = "./data/mSigAct/output/all_relatable_sbs_prop1_Activities.txt"
outpath = "./data/mSigAct/output/all_relatable_sbs_prop1_Activities.pdf"

df = pd.read_csv(inpath, sep='\t')
print(df.shape)
df['Set'] = df.Samples.apply(lambda x: x.split('_')[0])
df['WithTv'] = df.Samples.str.contains('Tv')
df['Samples'] = df['Samples'].str.replace('.only', '').str.replace('\.\.\.', ' & ')
df = df.sort_values(['Set', 'WithTv', 'Samples']).set_index('Set').loc[['high', 'low', 'diff']]
df = df.set_index('Samples').astype(int).drop('WithTv', 1)
# drop signatures with low fraction (< 0.5%)
df = df.loc[:, (df.T / df.sum(1)).T.mean() >= 0.005]
print(df.shape)
df.to_csv(inpath+'.sorted', sep='\t')

plotActivity(
    inpath+'.sorted', 
    outpath, 
    bin_size=30, 
    # custom_colors=custom_colors,
    delimiter_step=10, delimiter_size=3,
    # rename=True,
)



# only Ts
inpath  = "./data/mSigAct/output/all_relatable_sbs_prop1_Activities.txt"
outpath = "./data/mSigAct/output/all_relatable_sbs_prop1_Activities_Ts.pdf"

df = pd.read_csv(inpath, sep='\t')
print(df.shape)
df['Set'] = df.Samples.apply(lambda x: x.split('_')[0])
df['WithTv'] = df.Samples.str.contains('Tv')
df['Samples'] = df['Samples'].str.replace('__Ts.only', '').str.replace('_', '-')
df = df[~df['Samples'].str.contains('\.\.\.')]
df = df.sort_values(['Set', 'WithTv', 'Samples']).set_index('Set').loc[['high', 'low', 'diff']]
df = df.set_index('Samples').astype(int).drop('WithTv', 1)
# drop signatures with low fraction (< 0.5%)
df = df.loc[:, (df.T / df.sum(1)).T.mean() >= 0.005]
print(df.shape)
df.to_csv(inpath+'.sorted.ts', sep='\t')

plotActivity(
    inpath+'.sorted.ts', 
    outpath, 
    bin_size=50, 
    scale=True,
    # custom_colors=custom_colors,
    delimiter_step=5, delimiter_size=2,
    rename=True, figure_width=4,
)



### When using some priors
custom_colors = ['#63d69e', '#fad682', '#fad682', 'tab:purple', 
                 'tab:pink', 'tab:orange', 'lightgray']

inpath  = "./data/mSigAct/output/custom_prop_Activities.txt"
outpath = "./data/mSigAct/output/custom_prop_Activities.pdf"

df = pd.read_csv(inpath, sep='\t')
df['Set'] = df.Samples.apply(lambda x: x.split('_')[0])
df['Samples'] = df['Samples'].str.replace('.only', '').str.replace('\.\.\.', ' & ')
df = df.sort_values(['Set', 'Samples']).set_index('Set').loc[['high', 'low', 'diff']]
df = df.set_index('Samples', append=True).astype(int).reset_index(1)
df.to_csv(inpath+'.sorted', index=False, sep='\t')

plotActivity(
    inpath+'.sorted', 
    outpath, 
    bin_size=30, 
    # custom_colors=custom_colors,
    delimiter_step=10, delimiter_size=3,
    # rename=True,
)



### When using some priors
custom_colors = ['#63d69e', '#fad682', '#fad682', 'tab:purple', 
                 'tab:pink', 'tab:orange', 'lightgray']

inpath  = "./data/mSigAct/output/custom_from1_Activities.txt"
outpath = "./data/mSigAct/output/custom_from1_Activities.pdf"

df = pd.read_csv(inpath, sep='\t')
df['Set'] = df.Samples.apply(lambda x: x.split('_')[0])
df['Samples'] = df['Samples'].str.replace('.only', '').str.replace('\.\.\.', ' & ')
df = df.sort_values(['Set', 'Samples']).set_index('Set').loc[['high', 'low', 'diff']]
df = df.set_index('Samples', append=True).astype(int).reset_index(1)
df.to_csv(inpath+'.sorted', index=False, sep='\t')

plotActivity(
    inpath+'.sorted', 
    outpath, 
    bin_size=30, 
    # custom_colors=custom_colors,
    delimiter_step=10, delimiter_size=3,
    # rename=True,
)