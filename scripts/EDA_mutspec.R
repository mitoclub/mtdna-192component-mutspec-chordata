rm(list = ls(all=TRUE))

library(ggplot2)
library(dplyr)
library(tidyr)

obs_mut = read.csv('../data/ObsMutSpec.csv')
counted_codons = read.csv('../data/counted_codons_cytb.csv')
taxonomy = read.csv('../data/taxa_gb.csv', sep='\t')

com_sps = counted_codons$Species

obs_mut = obs_mut %>% filter(Species %in% com_sps)
obs_mut = merge(obs_mut, taxonomy, by = 'Species')
colnames(obs_mut)[4] = 'Pos3'

tr_vec = c('A>G','G>A', 'C>T', 'T>C')
tv_vec = c('C>A','A>C','G>T','T>G','T>A', 'A>T', 'C>G', 'G>C')

tr_or_tv = function(x)
{
  x = substr(x, 3, 5)
  if (x %in% tr_vec)
  {
    return('transition')
  }
  else if (x %in% tv_vec)
  {
    return('transversion')
  }
}

obs_mut$TsTv = apply(as.matrix(obs_mut$Mut3), 1, FUN = tr_or_tv)  

pdf(file='../pictures/EDA_mutspec.pdf')
to_bar = obs_mut %>% distinct(Species, .keep_all = TRUE) 

ggplot(to_bar, aes(x=Class, fill=Class))+
  geom_bar()+
  theme_classic()+
  ggtitle('Destribution of Species by Classes')


from_num_to_type = function(x)
{
  if (x == 0){return('NonSyn')}
  else if (x==1){return('Syn')}
  else if(x==2){return('4Fold')}
}

obs_mut$MutTypeFull = apply(as.matrix(obs_mut$MutType),1,from_num_to_type)
obs_mut$Pos3 = as.factor(obs_mut$Pos3)

df_pos3 = obs_mut[obs_mut$Pos3 == '1', ]
df_pos3$MutTypeFull = 'Pos3'

gr2 = rbind(obs_mut, df_pos3)
gr2 = gr2 %>% group_by(Class, MutTypeFull) %>%summarise(n=n())
gr2$MutTypeFull = as.factor(gr2$MutTypeFull)

ggplot(gr2, aes(x=Class, y=n, fill=MutTypeFull))+
  geom_bar(stat='identity', position=position_dodge())+
  theme_classic()+
  #scale_fill_manual(values = c('#FFB266','#CC99FF', '#66B2FF'))+
  ggtitle('Destribution of Mutations by Classes and Types of Mutaitons')


gr3 = obs_mut %>% group_by(Class) %>% summarise(n=n())

ggplot(gr3, aes(x=Class, y=n, fill=Class))+
  geom_bar(stat='identity')+
  theme_classic()+
  ggtitle('Destribution of all Mutations by Classes')


gr4 = obs_mut %>% group_by(Species) %>% summarise(n=n())

ggplot(gr4, aes(x=n))+
  geom_histogram(fill="lightblue", col='black', size=0.5)+
  theme_classic()+
  ggtitle('Destribution of all Mutations by Species')+
  scale_x_continuous(breaks = seq(0, 2800, 100), lim = c(0, 2800))+ 
  theme(axis.text.x = element_text(angle = 90))

gr5 = obs_mut %>% group_by(Species, MutTypeFull) %>% summarise(n=n())

ggplot(gr5, aes(x=MutTypeFull, y=n, fill=MutTypeFull))+
  geom_violin()+
  theme_classic()+
  ggtitle('Destribution of Mutations by Species and Mutation Type \n(with sps that have less then 500 mutations)')+
  scale_y_continuous(breaks = seq(0, 500, 20), lim = c(0, 500))
  #theme(axis.text.y = element_text(angle = 90))

gr6 = obs_mut %>% group_by(Class, Species,TsTv) %>% summarise(n=n())

ggplot(gr6, aes(x=Class, y=n, fill=TsTv))+
  geom_bar(stat='identity', position=position_dodge())+
  theme_classic()+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  ggtitle('Destribution of TsTv by Classes')


gr7 = obs_mut %>% group_by(Class,Species, TsTv) %>% summarise(n=n()) %>%
  pivot_wider(names_from=TsTv, values_from = n)

gr7[is.na(gr7)] <- 0
gr7$TsTv= gr7$transition/gr7$transversion
is.na(gr7)<-sapply(gr7, is.infinite)
gr7[is.na(gr7)]<-0

ggplot(gr7, aes(x=Class, y=TsTv, fill=Class))+
  geom_violin()+
  theme_classic()+
  ggtitle('Destribution of TsTv metric by Species \n(with sps that have less then 500 mutations)')+
  scale_y_continuous(breaks = seq(0, 60, 5), lim = c(0, 60))
  
dev.off()

