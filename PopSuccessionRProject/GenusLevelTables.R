#This script summarizes OTU table at the genus level
#i.e. sums all OTUs that have same genus assignment

library(tidyverse)

#make sure the data import + initial processing 
#has been done prior to running this script
source("LSU_DataImport.R")


Genus_mat_rare_rel=
  otu_mat_rare_rel %>%
  mutate(otu=rownames(.)) %>%
  left_join(select(tax_mat, otu, Genus)) %>% 
  group_by(Genus) %>%
  summarise(across(starts_with("bc"), sum, na.rm = TRUE), .groups = "drop") %>%
  mutate(Genus=case_when(is.na(Genus)~"Unassigned",
                         .default=Genus)) %>% 
  column_to_rownames("Genus") %>%
  t %>%
  as.data.frame %>% 
  mutate(sample=rownames(.)) %>%
  left_join(samples_df)

#creates occurrence counts of each genus by treatment
#and then makes columns indicating which genera have greater than
#or equal to 5 occurrences in one treatment and are absent in
#the other (treatment refers to species or time in this case)
Genus_occurrences=
  Genus_mat_rare_rel %>%
  filter(sample %in% FilterList) %>%
  mutate(TimeSpecies=paste(Time, Species, sep="_")) %>%
  group_by(TimeSpecies) %>%
  summarise(across(Agaricus:Wilcoxina, ~ sum(. != 0, na.rm = TRUE))) %>%
  ungroup %>%
  column_to_rownames("TimeSpecies") %>%
  t %>%
  as.data.frame %>% 
  mutate(June20=June20_A + June20_T,
         Jan21=Jan21_A + Jan21_T,
         A=Jan21_A+June20_A,
         T=Jan21_T+June20_T) %>% 
  mutate(TimeFilter=
           case_when(June20==0 & Jan21>=5 ~ "Jan21",
                     June20>5 & Jan21==0 ~ "June20",
                     .default=NA),
         SpeciesFilter=
           case_when(A==0 & T>=5 ~ "Ptrich",
                     A>5 & T==0 ~ "Ptrem",
                     .default=NA)) %>%
  mutate(FilterSum=
           case_when(!is.na(TimeFilter) & is.na(SpeciesFilter) ~"Time",
                     is.na(TimeFilter) & !is.na(SpeciesFilter) ~"Species",
                     !is.na(TimeFilter) & !is.na(SpeciesFilter) ~"TimeSpecies",
                     .default="NonSig"))

#creates list of genera that differ by species or time based
#on presence/absence. Then merges this with list of ANCOMBC2
#differentially abundant taxa. The only additional genus found
#to be differentially abundant by ANCOMBC2 (compared to 
#presence/absence) is Tulasnella
#
#Make sure that the script ANCOMBC2.R has been run to ensure 
#this list exists if needed:
#source("ANCOMBC2.R")
DiffGenera =
  Genus_occurrences %>%
  filter(!is.na(TimeFilter) | !is.na(SpeciesFilter)) %>% View
  rownames() %>%
  c(., ancom_diff_gen_list) %>%
  unique

#filter genus-level otu table to only include significant
#genera and converts to long format for plotting
Genus_mat_rare_rel_sig=
  Genus_mat_rare_rel %>%
  select(DiffGenera, sample) %>%
  gather("Genus", "Abundance", -sample) %>%
  left_join(samples_df) %>%
  filter(sample %in% FilterList)


ggplot(Genus_mat_rare_rel_sig, 
       aes(x = Time, y = Abundance, color = Species)) +
  #geom_boxplot(outlier.shape = NA) +  # Hide outliers to avoid overlap with points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 3, alpha = 0.6) + 
  facet_wrap(. ~ Genus, scales = "free_y") +
  theme_test()

#performs hierarchical clustering to determine plotting order
Genus_mat_rare_rel_sig_hclust=
  Genus_mat_rare_rel %>%
  select(DiffGenera, sample) %>%
  mutate(across(-sample, ~ (.-min(.)) / (max(.) - min(.)))) %>% 
  drop_na() %>%
  column_to_rownames("sample") %>% 
  t %>%
  dist %>% 
  hclust

#creates dataframe of genus plotting order -
#this somewhat follows the hierarchical clustering, but I made
#some adjustments for viewability
Genus_hclust_order=
  data.frame(Genus=c("Glomus",
                     "Rhizophagus",
                     "Sphaerosporella",
                     "Tulasnella",
                     "Exophiala",
                     "Aquilomyces",
                     "Cenococcum",
                     "Trichocladium",
                     "Solicoccozyma",
                     "Polyphilus",
                     "Diaporthe",
                     "Mycena"),
             GenusOrder=seq(1:12))

#filter genus-level otu table to only include significant
#genera and converts to long format for plotting
#data are scaled by Genus to make pretty heatmap
Genus_mat_rare_rel_sig_scale = Genus_mat_rare_rel %>%
  select(DiffGenera, sample) %>%
  mutate(across(-sample, ~ (.-min(.)) / (max(.) - min(.)))) %>% 
  gather("Genus", "Abundance", -sample) %>%
  left_join(samples_df) %>%
  filter(sample %in% FilterList) %>%
  mutate(sample_mod = paste(Time, Species, sample, sep="_"),
         Time = factor(Time, levels = c("June20", "Jan21"))) %>%
  left_join(Genus_hclust_order) %>%
  arrange(Time, Species, GenusOrder) %>%
  mutate(
    Time = factor(Time, levels = c("June20", "Jan21")), # Ensure order
    Species = factor(Species, levels = unique(Species)),
    sample_mod = factor(sample_mod, levels = unique(sample_mod)), # Ensure sample_mod respects the order
    Genus=factor(Genus, levels = unique(Genus))
    )

#makes heatmap of scaled abundance of differentially abundant genera
diffabun_heat=
  ggplot() +
  # Plot white tiles for zero Abundance values
  geom_tile(data = filter(Genus_mat_rare_rel_sig_scale, Abundance == 0),
            aes(x = sample_mod, y = Genus), fill = "white") +
  # Plot colored tiles for non-zero Abundance
  geom_tile(data = filter(Genus_mat_rare_rel_sig_scale, Abundance != 0),
            aes(x = sample_mod, y = Genus, fill = Abundance)) +
  # Add thin white lines to separate rows (ensures no gap)
  geom_hline(yintercept = seq(0.5, length(unique(Genus_mat_rare_rel_sig_scale$Genus)) + 0.5, by = 1), 
             color = "grey40", size = 0.25) +
  theme_test() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  scale_y_discrete(expand = c(0, 0)) + # Remove padding between tiles
  scale_x_discrete(expand = c(0, 0)) + # Remove padding between columns
  facet_wrap(Time ~ Species, scales = "free_x", nrow = 1) +
  theme(
    axis.text.x = element_blank(),        
    axis.title.x= element_blank(),
    axis.ticks.x = element_blank(),        
    strip.text = element_blank(),
    axis.text.y = element_text(size = 8)
  )

#outputs figure as PDF so that it can be 
#formatted in Adobe Illustrator and cleaned up
pdf(file="./OutputFigures/DiffAbun_heat.pdf", 
    width = 11, height = 2.2)
diffabun_heat
dev.off()


Genus_mat_rare_rel_sig %>%
  filter(Genus=="Aquilomyces") %>%
ggplot(aes(x=Stand, y=Abundance)) +
  geom_boxplot() +
  facet_wrap(.~Time, scales="free_y")



Family_mat_rare_rel=
  otu_mat_rare_rel %>%
  mutate(otu=rownames(.)) %>%
  left_join(select(tax_mat, otu, Family)) %>% 
  group_by(Family) %>%
  summarise(across(starts_with("bc"), sum, na.rm = TRUE), .groups = "drop") %>%
  mutate(Family=case_when(is.na(Family)~"Unassigned",
                         .default=Family)) %>% 
  column_to_rownames("Family") %>%
  t %>%
  as.data.frame %>% 
  mutate(sample=rownames(.)) %>%
  left_join(samples_df)

ggplot(filter(Family_mat_rare_rel, sample %in% FilterList), 
       aes(x=Time, y=Thelephoraceae, color=Species)) +
  geom_boxplot()
  
ggplot(filter(Genus_mat_rare_rel, sample %in% FilterList), 
       aes(x=Time, y=Sphaerosporella, color=Species)) +
  geom_boxplot()


filter(Genus_mat_rare_rel, sample %in% FilterList) %>%
  lm(Thelephora~Species*Time,.) %>%
  Anova

Colonization_df %>%
  filter(Stand!="Control") %>%
  mutate(TotalEndophye=Microsclerotia+Endophyte) %>%
  ggplot(aes(x=Time, y=AM, color=Species)) +
  geom_boxplot()

Colonization_df %>%
  filter(Stand!="Control") %>%
  mutate(TotalEndophye=Microsclerotia+Endophyte) %>%
  lm(AM~Species*Time,.) %>%
  Anova


