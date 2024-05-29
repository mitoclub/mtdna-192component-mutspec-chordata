library(ICAMS)
library(mSigAct)
library(cosmicsig)
library(dplyr)

path = '../data/decomp'
samples_names = list.files(path, pattern='_samples.txt', full.names=T)
samples_names

diff_df = read.table("../data/decomp/high_minus_low_Ts_samples.txt", sep='\t', header =T)
high_df = read.table("../data/decomp/high_Ts_samples.txt" , sep='\t', header =T)
low_df = read.table("../data/decomp/low_Ts_samples.txt", sep='\t', header =T)

diff_df = diff_df %>% 
          rename_at(2:ncol(diff_df),~ paste0('diff', '_', .))

high_df = high_df %>% 
  rename_at(2:ncol(high_df),~ paste0('high', '_', .)) 

low_df = low_df %>% 
  rename_at(2:ncol(low_df),~ paste0('low', '_', .)) 

low_df = low_df[,-1]
high_df = high_df[,-1]

final_catalog = cbind(diff_df, high_df, low_df)

write.table(final_catalog, '../data/decomp/samples_mSigAct.txt', sep='\t', row.names = F)

path_full_catalog = '../data/decomp/samples_mSigAct.txt'
input_catalog <- ICAMS::ReadCatalog(file = path_full_catalog)


sigs <- cosmicsig::COSMIC_v3.3$signature$GRCh37$SBS96

sig_use= c('SBS2', 'SBS23', 'SBS44', 'SBS5', 'SBS30', 'SBS12')

sigs_to_use <- sigs[, colnames(sigs) %in% sig_use]


### Proportions are taken from SigProfiler

output_home <- "../data/decomp/mSigAct/custom_prop"

sig_prop = c(0.01, 0.03, 0.06, 0.18, 0.33, 0.38)
names(sig_prop) = colnames(sigs_to_use)

retval <-
  mSigAct::MAPAssignActivity(spectra = input_catalog, 
                             sigs = sigs_to_use,
                             sigs.presence.prop = sig_prop,
                             output.dir = output_home,
                             max.level = ncol(sigs_to_use) - 1,
                             p.thresh = 0.05 / ncol(sigs_to_use), 
                             num.parallel.samples = 2, 
                             mc.cores.per.sample = 5)

### All proportions are the same == 1

output_home <- "../data/decomp/mSigAct/prop1"

sig_prop = c(1, 1, 1, 1, 1, 1)
names(sig_prop) = colnames(sigs_to_use)

retval <-
  mSigAct::MAPAssignActivity(spectra = input_catalog, 
                             sigs = sigs_to_use,
                             sigs.presence.prop = sig_prop,
                             output.dir = output_home,
                             max.level = ncol(sigs_to_use) - 1,
                             p.thresh = 0.05 / ncol(sigs_to_use), 
                             num.parallel.samples = 2, 
                             mc.cores.per.sample = 5)


### All proportions are the same == 0.5

output_home <- "../data/decomp/mSigAct/prop05"

sig_prop = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
names(sig_prop) = colnames(sigs_to_use)

retval <-
  mSigAct::MAPAssignActivity(spectra = input_catalog, 
                             sigs = sigs_to_use,
                             sigs.presence.prop = sig_prop,
                             output.dir = output_home,
                             max.level = ncol(sigs_to_use) - 1,
                             p.thresh = 0.05 / ncol(sigs_to_use), 
                             num.parallel.samples = 2, 
                             mc.cores.per.sample = 5)


### Proportions from mSigAct with 1 prop

output_home <- "../data/decomp/mSigAct/custom_from1"


sig_prop = c(0.047459, 0.236982, 0.390326, 0.224110, 0.225671, 0.018647)
names(sig_prop) = colnames(sigs_to_use)

retval <-
  mSigAct::MAPAssignActivity(spectra = input_catalog, 
                             sigs = sigs_to_use,
                             sigs.presence.prop = sig_prop,
                             output.dir = output_home,
                             max.level = ncol(sigs_to_use) - 1,
                             p.thresh = 0.05 / ncol(sigs_to_use), 
                             num.parallel.samples = 2, 
                             mc.cores.per.sample = 5)


### Take all relatable SBS with 1 prop 
output_home <- "../data/decomp/mSigAct/all_relatable_sbs"


sigs <- cosmicsig::COSMIC_v3.3$signature$GRCh37$SBS96

sig_to_delete= c('SBS32', 'SBS11', 'SBS25', 'SBS31', 'SBS32', 'SBS35', 'SBS86',
                 'SBS87', 'SBS90', 'SBS22', 'SBS88', 'SBS27', 'SBS43', 'SBS45',
                 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50', 'SBS51', 'SBS52',
                 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58', 'SBS59',
                 'SBS60', 'SBS95', 'SBS9', 'SBS84', 'SBS85')

sigs_to_use = sigs[, !(colnames(sigs) %in% sig_to_delete)] 

sig_prop = rep(1, ncol(sigs_to_use))
names(sig_prop) = colnames(sigs_to_use)

retval <-
  mSigAct::MAPAssignActivity(spectra = input_catalog, 
                             sigs = sigs_to_use,
                             sigs.presence.prop = sig_prop,
                             output.dir = output_home,
                             max.level = ncol(sigs_to_use) - 1,
                             p.thresh = 0.05 / ncol(sigs_to_use), 
                             num.parallel.samples = 2, 
                             mc.cores.per.sample = 8, 
                             max.subsets = 100)
