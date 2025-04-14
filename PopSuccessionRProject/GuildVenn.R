#Calculates number of samples classified as ECM, AM, or dual-mycorrhizal
#based on sequencing and microscopy. For sequencing, a cutoff of 
#0.5% sequence abundance is required for a positive detection of
#a mycorrhizal guild

library(tidyverse)

###################Sequencing Based##########################
#sets guild abundance threshold
guild_cut=0.001

Guild_Venn=
  GuildSumm_wide %>%
  mutate(MycStat=
           case_when(
             `Arbuscular Mycorrhizal` > guild_cut & 
               Ectomycorrhizal < guild_cut ~ "AM",
             `Arbuscular Mycorrhizal` > guild_cut & 
               Ectomycorrhizal > guild_cut ~ "DM",
             `Arbuscular Mycorrhizal` < guild_cut & 
               Ectomycorrhizal > guild_cut ~ "EM",
             `Arbuscular Mycorrhizal` < guild_cut & 
               Ectomycorrhizal < guild_cut ~ "NM",
             .default="Invalid"
           ))

Guild_Venn_sum=
  Guild_Venn %>%
  group_by(Species, Time, MycStat) %>%
  summarise(MycStat_abs = n(), .groups = 'drop') %>% 
  ungroup() %>% 
  group_by(Species, Time) %>%
  mutate(MycStat_rel = MycStat_abs / sum(MycStat_abs))

ggplot(Guild_Venn_sum, aes(x=Time, y =MycStat_rel, fill=MycStat)) +
  geom_bar(stat="identity") +
  facet_wrap(.~Species)

#plotting reveals that there are a few Trichocarpa samples that
#were non-EM at 7 months that became EM at 12 months.
#Identifies samples that went from Non-EM to EM to look for evidence
#of contamination
Trich_trans_EM=
  Guild_Venn %>%
  filter(Species=="T") %>%
  select(MycStat, Pot, Time) %>% 
  spread("Time", "MycStat")  

  
###################Microscopy Based##########################
View(Colonization_bysample)

Guild_Venn_mic=
  Colonization_bysample %>%
  filter(!is.na(Ectomycorrhizal)) %>%
  mutate(MycStat=
           case_when(
             `Arbuscular Mycorrhizal` > 0 & 
               Ectomycorrhizal == 0  ~ "AM",
             `Arbuscular Mycorrhizal` > 0 & 
               Ectomycorrhizal > 0 ~ "DM",
             `Arbuscular Mycorrhizal` == 0 & 
               Ectomycorrhizal > 0 ~ "EM",
             `Arbuscular Mycorrhizal` == 0 & 
               Ectomycorrhizal == 0 ~ "NM",
             .default="Invalid"
           )) %>%
  mutate(Time=factor(Time, levels=c("7Months", "12Months")))

Guild_Venn_mic_sum=
  Guild_Venn_mic %>%
  group_by(Species, Time, MycStat) %>%
  summarise(MycStat_abs = n(), .groups = 'drop') %>% 
  ungroup() %>% 
  group_by(Species, Time) %>%
  mutate(MycStat_rel = MycStat_abs / sum(MycStat_abs))

ggplot(Guild_Venn_mic_sum, aes(x=Time, y =MycStat_rel, fill=MycStat)) +
  geom_bar(stat="identity") +
  facet_wrap(.~Species)


