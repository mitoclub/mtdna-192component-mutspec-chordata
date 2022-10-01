rm(list=ls(all=TRUE))

library(dplyr)
library(ggplot2)

mutspec = read.csv('../data/MutSpecVertebratescytb.csv')
mutspec = mutspec %>% select(-c(RawMutSpecSum, MutSpec))

ms_cls = mutspec %>% group_by(Class, Mut) %>%
  summarise(RawMutSpecSum = sum(RawMutSpec)) %>% 
  
ms_cls = merge(ms_cls, mutspec %>% select(Class, Mut, RawMutSpec, MutBase),
               by = c('Class', 'Mut'))

ms_CT = ms_cls %>% filter(MutBase == 'C>T')
ms_AG = ms_cls %>% filter(MutBase == 'A>G')
