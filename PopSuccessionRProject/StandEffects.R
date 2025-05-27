#Tests for effects of the stand that the soil inoculum was 
#sourced from on various parameters
library(tidyverse)
library(car)

View(GuildSumm_wide)

Colonization_bysample %>%
  filter(Species=="T", Time=="7Months") %>%
  ggplot(aes(x=Stand, y=TotalEndophyte)) +
  geom_boxplot()

Colonization_bysample %>%
  filter(Species=="T", Time=="7Months") %>%
  lm(sqrt(`Arbuscular Mycorrhizal`)~Stand,.) %>%
  Anova
