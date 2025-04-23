#This script tests for an effect of field soil
#inoculum levels on alpha and beta diversity
#This was done because some pots had lower 
#amounts of soil inoculum used because we ran out

library(tidyverse)
library(car)

#make sure the data import + initial processiong 
#has been done prior to running this script
source("LSU_DataImport.R")

AlphaDiv %>%
  filter(Species!="Neg", Stand!="C") %>%
  lm(Richness~Time+Species+FieldSoilFraction,.) %>%
  Anova

#scatter plot to look at richness vs field soil fraction
ggplot(filter(AlphaDiv, Species!="Neg", Stand!="C"), 
       aes(x=FieldSoilFraction, y=Shannon)) +
  geom_point() +
  theme_test() +
  facet_grid(Species~Time)

#pcoa to look at field soil fraction
ggplot(LSU_bray_pcoa_df, 
       aes(x=V1, y=V2, color=FieldSoilFraction))+
  geom_point(size=2) +
  facet_grid(Species~Time)

