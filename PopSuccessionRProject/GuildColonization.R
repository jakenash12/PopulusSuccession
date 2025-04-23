#This script calculates percent colonization of three
#fungal guilds (endophytes, ectomycorrhizal, and
#arbuscular mycorrhizal) based on sequencing and
#microscopy data, then makes a stacked bar chart

library(tidyverse)

#make sure the data import + initial processing 
#has been done prior to running this script
source("LSU_DataImport.R")

#sets input path
Input_path="./InputFiles/"

#reads in microscopic fungal colonization data
Colonization_df= read.delim(paste(Input_path, "FungalColonization.txt", sep=""))

#calcualates guild colonization from sequence abundance per sample
#in long format
GuildSumm=
  otu_mat_rare_rel %>%
  mutate(otu=rownames(.)) %>%
  gather(sample, abundance, -otu) %>%
  left_join(select(tax_mat_fg, otu, primary_guild)) %>%
  mutate(primary_guild=
           case_when(
             primary_guild %in% c("Ectomycorrhizal",
                                  "DSE",
                                  "Arbuscular Mycorrhizal") ~ primary_guild,
             primary_guild=="Endophyte" ~ "Non-DSE Endophyte",
             .default="Other")) %>% 
  mutate(primary_guild = case_when(
    primary_guild %in% c("DSE", "Non-DSE Endophyte") ~ "TotalEndophyte",
    TRUE ~ primary_guild
  )) %>%
  group_by(primary_guild, sample) %>%
  summarise(abundance=sum(abundance))

#converts guild summary to wide format
GuildSumm_wide=
  GuildSumm %>%
  spread(primary_guild, abundance) %>%
  left_join(samples_df) %>%
  mutate(Time = dplyr::recode(Time, "June20" = "7Months", "Jan21" = "12Months")) %>%
  mutate(SpeciesTime=paste(Species, Time, sep="_")) %>%
  filter(sample %in% FilterList) 
  

#calculates mean sequence based guild colonization
#per treatment (i.e. time point and tree species)
GuildMeans=
  GuildSumm %>%
  left_join(samples_df) %>%
  filter(sample %in% FilterList) %>% 
  group_by(Time, Species, primary_guild) %>%
  summarise(abundance=mean(abundance)) %>%
  mutate(Time=factor(Time, levels=c("June20","Jan21")),
         primary_guild=factor(primary_guild, levels=c(
           "Arbuscular Mycorrhizal",
           "Ectomycorrhizal",
           "TotalEndophyte",
           "Other"
         )))  %>%
  mutate(Time = dplyr::recode(Time, "June20" = "7Months", "Jan21" = "12Months")) %>%
  mutate(Method="Sequencing")

GuildMeans_species=
  GuildSumm %>%
  left_join(samples_df) %>%
  filter(sample %in% FilterList) %>% 
  group_by(Species, primary_guild) %>%
  summarise(meanabundance=mean(abundance),
            sd_abundance=sd(abundance)) %>%
  mutate(primary_guild=factor(primary_guild, levels=c(
           "Arbuscular Mycorrhizal",
           "Ectomycorrhizal",
           "TotalEndophyte",
           "Other"
         )))  %>%
  mutate(Method="Sequencing")

#creates mean fungal colonization based on microscopy
#per treatment
Colonization_bysample=
  Colonization_df %>%
  filter(Stand!="Control") %>%
  rename("Ectomycorrhizal"="ECM",
         "Arbuscular Mycorrhizal"="AM") %>%
  mutate(TotalEndophyte=Endophyte,
         SpeciesTime=paste(Species, Time, sep="_"))
  
Colonization_means=
  Colonization_bysample %>% 
  filter(!is.na(Uncolonized)) %>% View
  select(-c("Microsclerotia","Endophyte","Stand", "SpeciesTime")) %>% 
  gather("primary_guild","abundance", -Species, -Time, -SampleID) %>% View
  group_by(Species,Time, primary_guild) %>% 
  summarise(abundance=mean(abundance)) %>%
  select(Time,Species,primary_guild,abundance) %>%
  mutate(Method="Microscopy",
         abundance=abundance/100,
         primary_guild=factor(primary_guild, levels=c(
           "Arbuscular Mycorrhizal",
           "Ectomycorrhizal",
           "TotalEndophyte",
           "Uncolonized"
         )))

#calculates average and sd of guild colonization by species (not by time)
Colonization_means_byspecies=
  Colonization_bysample %>%
  select(-c("Microsclerotia","Endophyte","Stand", "SpeciesTime")) %>% 
  gather("primary_guild","abundance", -Species, -Time, -SampleID) %>%
  filter(!is.na(abundance)) %>%
  group_by(Species, primary_guild) %>%
  summarise(meanabundance=mean(abundance),
            sd_abundance=sd(abundance)) %>%
  select(Species,primary_guild,meanabundance, sd_abundance) %>%
  mutate(Method="Microscopy",
         primary_guild=factor(primary_guild, levels=c(
           "Arbuscular Mycorrhizal",
           "Ectomycorrhizal",
           "TotalEndophyte",
           "Uncolonized"
         )))

#binds together sequence and microscopy based colonization
#assessment dataframes so that they can be plotted together
sequencing_microscopy_df=
  rbind(Colonization_means, GuildMeans) %>%
  filter(Species!="Neg") %>%
  mutate(Time=factor(Time, levels=c("7Months","12Months")))

#custom color palette for stacked barplot
guild_color_map=color_map <- c("Ectomycorrhizal" = "#603540",
                               "Arbuscular Mycorrhizal" = "#c4462b",
                               "TotalEndophyte"="#19424e",
                               "Other" = "#827e73",
                               "Uncolonized" = "burlywood3")

GuildBarplot=
  ggplot(sequencing_microscopy_df, aes(x=Species, 
                                     y=abundance, 
                                     fill=primary_guild)) +
  geom_bar(stat="identity", position="stack") +
  facet_wrap(Method~Time, nrow=1) +
  scale_fill_manual(values = guild_color_map) +
  theme_test() +
  scale_y_continuous(limits = c(0,1.01), expand = expansion(mult = c(0, .05)))

#outputs figure as PDF so that it can be 
#formatted in Adobe Illustrator and cleaned up
pdf(file="./OutputFigures/GuildBarplotRaw.pdf", 
    width = 8, height = 4.5)
GuildBarplot
dev.off()


###Guild colonization ANOVA on sequencing and microscopy data####


###################Microscopy#######################
#creates subsets of data by timepoint
Colonization_bysample_7months=
  filter(Colonization_bysample, Time=="7Months")

Colonization_bysample_12months=
  filter(Colonization_bysample, Time=="12Months")

#ECM tests
mod_col_ECM=lm(sqrt(Ectomycorrhizal)~ Species*Time, Colonization_bysample)
car::Anova(mod_col_ECM)

pairwise.wilcox.test(Colonization_bysample$Ectomycorrhizal,
                     Colonization_bysample$SpeciesTime, p.adjust.method = "fdr")

hist(resid(mod_col_ECM))
shapiro.test(resid(mod_col_ECM))

mod_col_ECM7=lm(sqrt(Ectomycorrhizal)~ Species, Colonization_bysample_7months)
car::Anova(mod_col_ECM7)

mod_col_ECM12=lm(sqrt(Ectomycorrhizal)~ Species, Colonization_bysample_12months)
car::Anova(mod_col_ECM12)

#AM tests
mod_col_AM=lm(sqrt(`Arbuscular Mycorrhizal`)~ Species*Time, Colonization_bysample)
car::Anova(mod_col_AM)

pairwise.wilcox.test(Colonization_bysample$`Arbuscular Mycorrhizal`,
                     Colonization_bysample$SpeciesTime, p.adjust.method = "fdr")


hist(resid(mod_col_AM))
shapiro.test(resid(mod_col_AM))

mod_col_AM7=lm(sqrt(`Arbuscular Mycorrhizal`)~ Species, Colonization_bysample_7months)
car::Anova(mod_col_AM7)

mod_col_AM12=lm(sqrt(`Arbuscular Mycorrhizal`)~ Species, Colonization_bysample_12months)
car::Anova(mod_col_AM12)

#Endophyte tests
mod_col_En=lm(sqrt(TotalEndophyte)~ Species*Time, Colonization_bysample)
car::Anova(mod_col_En)

pairwise.wilcox.test(Colonization_bysample$TotalEndophyte,
                     Colonization_bysample$SpeciesTime, p.adjust.method = "fdr")


hist(resid(mod_col_En))
shapiro.test(resid(mod_col_En))

mod_col_En7=lm(sqrt(TotalEndophyte)~ Species, Colonization_bysample_7months)
car::Anova(mod_col_En7)

mod_col_En12=lm(sqrt(TotalEndophyte)~ Species, Colonization_bysample_12months)
car::Anova(mod_col_En12)


########################Sequencing Data###################
#creates subsets of data by timepoint
GuildSumm_wide_7months=
  filter(GuildSumm_wide, Time=="7Months")

GuildSumm_wide_12months=
  filter(GuildSumm_wide, Time=="12Months")

#ECM tests
mod_LSU_ECM=lm(sqrt(Ectomycorrhizal)~ Species*Time, GuildSumm_wide)
car::Anova(mod_LSU_ECM)

pairwise.wilcox.test(GuildSumm_wide$Ectomycorrhizal,
                     GuildSumm_wide$SpeciesTime, p.adjust.method = "fdr")

pairwise.t.test(GuildSumm_wide$Ectomycorrhizal,
                GuildSumm_wide$SpeciesTime, p.adjust.method = "fdr")


hist(resid(mod_LSU_ECM))
shapiro.test(resid(mod_LSU_ECM))

mod_LSU_ECM7=lm(sqrt(Ectomycorrhizal)~ Species, GuildSumm_wide_7months)
car::Anova(mod_LSU_ECM7)

mod_LSU_ECM12=lm(sqrt(Ectomycorrhizal)~ Species, GuildSumm_wide_12months)
car::Anova(mod_LSU_ECM12)

#AM tests
mod_LSU_AM=lm(sqrt(`Arbuscular Mycorrhizal`)~ Species * Time, GuildSumm_wide)
car::Anova(mod_LSU_AM)

pairwise.wilcox.test(GuildSumm_wide$`Arbuscular Mycorrhizal`,
                     GuildSumm_wide$SpeciesTime, p.adjust.method = "fdr")

hist(resid(mod_LSU_AM))
shapiro.test(resid(mod_LSU_AM))

mod_LSU_AM7=lm(sqrt(`Arbuscular Mycorrhizal`)~ Species, GuildSumm_wide_7months)
car::Anova(mod_LSU_AM7)

mod_LSU_AM12=lm(sqrt(`Arbuscular Mycorrhizal`)~ Species, GuildSumm_wide_12months)
car::Anova(mod_LSU_AM12)

#Endophyte tests
mod_LSU_En=lm(sqrt(TotalEndophyte)~ Species*Time, GuildSumm_wide)
car::Anova(mod_LSU_En)

hist(resid(mod_LSU_En))
shapiro.test(resid(mod_LSU_En))

mod_LSU_En7=lm(sqrt(TotalEndophyte)~ Species, GuildSumm_wide_7months)
car::Anova(mod_LSU_En7)

mod_LSU_En12=lm(sqrt(TotalEndophyte)~ Species, GuildSumm_wide_12months)
car::Anova(mod_LSU_En12)


