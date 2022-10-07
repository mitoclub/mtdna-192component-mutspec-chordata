rm(list=ls(all=TRUE))

library(dplyr)
library(ggplot2)

mutspec = read.csv('../data/MutSpecVertebratescytb.csv')

ms_CT = mutspec %>% filter(MutBase == 'C>T')
ms_AG = mutspec %>% filter(MutBase == 'A>G')

ggplot(data = ms_CT, aes(x=Mut, y=MutSpec,color=Class))+
  geom_violin()
