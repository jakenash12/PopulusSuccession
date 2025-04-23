#Tests for correlations between abundance of different fungal guilds
#in the sequencing and microscopy dataset

library(tidyverse)
library(car)

GuildSumm_wide %>%
  filter(Species=="T") %>% 
  ggplot(aes(x=Ectomycorrhizal, y=TotalEndophyte, color=Time)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  theme_test() +
  facet_grid(.~Time, scales="free_x")


GuildSumm_wide %>%
  filter(Species=="T", Time=="12Months") %>% 
  lm((Ectomycorrhizal)~ `Arbuscular Mycorrhizal`, .) %>% 
  Anova

Colonization_bysample %>%
  filter(!is.na(Ectomycorrhizal))
