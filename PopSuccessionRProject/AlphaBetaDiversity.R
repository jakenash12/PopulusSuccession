library(vegan)
library(tidyverse)
library(car)

#make sure the data import + initial processiong 
#has been done prior to running this script
source("LSU_DataImport.R")

#calculates alpha diversity metrics
AlphaDiv=
  data.frame(Shannon=diversity(otu_mat_rare_t, index="shannon"),
             Simpson=diversity(otu_mat_rare_t, index="simpson"),
             Invsimpson=diversity(otu_mat_rare_t, index="invsimpson"),
             Richness=specnumber(otu_mat_rare_t)) %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df)

AlphaDiv %>%
  filter(Species!="Neg", Stand!="C") %>%
  group_by(Time, Species) %>%
  summarise(mean_rich=mean(Richness),
            sd_rich=sd(Richness))

Anova(lm(Richness~Species*Time, AlphaDiv))

ggplot(filter(AlphaDiv, Species!="Neg", Stand!="C"), 
       aes(x=Time, y=Richness, color=Species)) +
  geom_boxplot()

AlphaDiv %>%
  filter(Species!="Neg", Stand!="C") %>%
  lm(Richness~Species*Time,.) %>%
  Anova


LSU_bray_dm=
  otu_mat_rare_t %>%
  filter(rownames(.) %in% FilterList) %>%
  metaMDSdist(distance="bray")

LSU_bray_pcoa=cmdscale(LSU_bray_dm,k=3)

LSU_bray_pcoa_df=
  LSU_bray_pcoa %>%
  as.data.frame %>% 
  mutate(sample=rownames(.)) %>%
  left_join(samples_df)

ggplot(LSU_bray_pcoa_df, 
       aes(x=V1, y=V2, color=Species, shape=Species))+
  geom_point(size=2) +
  facet_wrap(.~Time)



