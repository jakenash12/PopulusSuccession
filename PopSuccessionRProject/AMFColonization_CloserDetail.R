#This script looks at the type of AM colonization (microscopy)
#in more detail. Examines whether arbuscular, coils, or vesicles
#are the predominant form of colonization
#
#Also looks for effect of mycorrhizal colonization on plant growth
#(i.e. height)

library(tidyverse)
library(magrittr)
library(car)

#make sure the data import + initial processing 
#has been done prior to running this script
source("LSU_DataImport.R")

#sets input path
Input_path="./InputFiles/"

#reads in detailed dataframe of fungal colonization
#(i.e. has info on specific structures like arbuscules, vesicles, coils) 
#in addition to total colonization by mycorrhizal type
Colonization_df_detailed= read.delim(paste(Input_path, "FungalColonization_detailed.txt", sep=""))

#calculates what proportion of vesicles, coils, arbuscules, and 
#AMF hyphae contribute to total AMF colonization
AMF_detail=
  Colonization_df_detailed %>%
  filter(TotalAMColonization>0, Species=="T") %>%
  mutate(Arb_prop=ArbusculeColonization/TotalAMColonization,
         Ves_prop=VesicleColonization/TotalAMColonization,
         Coil_prop=CoilColonization/TotalAMColonization,
         AMH_prop=AMHyphaeColonization/TotalAMColonization) %>%
  select(SampleID, HarvestDate, Species,
         Arb_prop, Ves_prop, Coil_prop, AMH_prop, TotalAMColonization) %>%
  mutate(HarvestDate=factor(HarvestDate, levels=c("Ju2020", "Ja2021")))

#calculates means of above df
colMeans(AMF_detail[sapply(AMF_detail, is.numeric)], na.rm = TRUE)

Colonization_df_detailed_summ=
  Colonization_df_detailed %>%
  filter(Species=="T")  %>%
  summarise(across(where(is.numeric), list(mean = mean, sd = sd), na.rm = TRUE))


###########Testing for differences in plant height between
##AM, EM, dual-mycorrhizal, and non-mycorrhizal plants
Colonization_df_detailed_growth <- Colonization_df_detailed %>%
  mutate(MycStatus = case_when(
    grepl("Control", SampleID) ~ "NM",
    TotalAMColonization > 0 & ECMColonization > 0 ~ "DM",
    TotalAMColonization > 0 & ECMColonization == 0 ~ "AM",
    TotalAMColonization == 0 & ECMColonization > 0 ~ "EM",
    TotalAMColonization == 0 & ECMColonization == 0 ~ "NM",
    .default = NA
  )) %>%
  mutate(
    AM = case_when(
      grepl("Control", SampleID) ~ "AMF",
      TotalAMColonization > 0 ~ "AMT",
      TRUE ~ "AMF"
    ),
    EM = case_when(
      grepl("Control", SampleID) ~ "EMF",
      ECMColonization > 0 ~ "EMT",
      TRUE ~ "EMF"
    )
  ) %>%
  mutate(Cont = case_when(
    grepl("Control", SampleID) ~ "Control",
    TRUE ~ "NotControl"
  )) %>%
  mutate(MycStatus = factor(MycStatus, levels = c("NM", "AM", "EM", "DM")))

Colonization_df_detailed_growth %>%
  filter(Species == "T", Cont == "NotControl") %>%
  ggplot(aes(x=MycStatus, y=Height)) +
  geom_boxplot() +
  facet_grid(.~HarvestDate)

Colonization_df_detailed_growth %>%
  filter(Species=="T", 
         HarvestDate=="Ja2021") %>%
  lm(Height~MycStatus,.) %>%
  Anova

#pairwise t-test for growth differences between myco types
#in p trichocarpa
Colonization_df_detailed_growth %>%
  filter(Species == "T", HarvestDate == "Ja2021", MycStatus!="NM") %>%
  with(pairwise.t.test(Height, MycStatus, p.adjust.method = "fdr"))

#calculates the frequency of mycorrhizal status in trichocarpa
Trich_mycstatus_summary=
  Colonization_df_detailed_growth %>%
  filter(Species == "T", Cont == "NotControl") %>%
  group_by(MycStatus) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(RelFreq = Count / sum(Count)) %>%
  ungroup()
  
  