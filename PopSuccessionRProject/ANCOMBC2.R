library(ANCOMBC)
library(dplyr)
library(phyloseq)

#make sure the data import + initial processing 
#has been done prior to running this script
source("LSU_DataImport.R")

#####first we identify taxa that are present
#in one condition and absent in the other
#we define this by the genus occuring in at least
#five samples in one condition and being absent from
#the other

#first create genus level summary table
otu_mat_rare_genus=
  otu_mat_rare %>% 
  mutate(otu=rownames(.)) %>%
  left_join(select(tax_mat, otu, Genus)) %>%
  mutate(Genus=case_when(
    is.na(Genus)~"Unassigned",
    .default=Genus
  )) %>%
  column_to_rownames("otu") %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  column_to_rownames("Genus")

#creates occurrence counts of each genus
#creates columns calculating significance
otu_mat_rare_genus_counts=
  otu_mat_rare_genus %>%
  t %>%
  as.data.frame %>%
  filter(rownames(.) %in% FilterList) %>%
  mutate(sample=rownames(.)) %>%
  left_join(select(samples_df, sample, Species, Time)) %>%  
  select(-sample) %>%
  group_by(Species, Time) %>%
  summarise(across(where(is.numeric), ~ sum(. != 0))) %>% 
  ungroup %>%
  mutate(SpeciesTime=paste(Species, Time, sep="_")) %>%
  select(-Species, -Time) %>% 
  column_to_rownames("SpeciesTime") %>%
  t %>% 
  as.data.frame %>% 
  mutate(A=A_June20+A_Jan21,
         T=T_June20+T_Jan21,
         June20=A_June20+T_June20,
         Jan21=A_Jan21+T_Jan21) %>%
  mutate(TimeSig=case_when(
    June20>=5 & Jan21==0 ~ "June20",
    June20==0 & Jan21>=5 ~ "Jan21",
    .default="NA"
  ),
  SpeciesSig=case_when(
    A>=5 & T==0 ~ "Tremuloides",
    A==0 & T>=5 ~ "Trichocarpa",
    .default="NA"
  ))

#creates phyloseq object from rarefied otu table
otu_mat_rare_phyloseq=
  otu_mat_rare_t %>%
  filter(rownames(.) %in% FilterList) %>% 
  otu_table(taxa_are_rows = FALSE)

#creates phyloseq format metadata table
samples_phyloseq=
  samples_df %>% 
  filter(sample %in% FilterList) %>%
  column_to_rownames("sample") %>%
  sample_data

#creates phyloseq format taxonomy matrix
tax_mat_phyloseq =
  tax_mat %>%
  column_to_rownames("otu") %>%
  as.matrix %>%
  tax_table

#combines tax mat, otu table and sample df into 
#phyloseq object
Phyloseq_df=phyloseq(otu_mat_rare_phyloseq,
                     samples_phyloseq,
                     tax_mat_phyloseq)

Phyloseq_test = ancombc2(data = Phyloseq_df,
                                 fix_formula = "Species+Time",
                                 p_adj_method = "fdr", pseudo_sens = TRUE,
                                prv_cut = 0.02, lib_cut = 0, s0_perc = 0.01,
                                tax_level="Genus",
                                 group = "Species", struc_zero = FALSE, neg_lb = FALSE,
                                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                                 global = TRUE, pairwise = FALSE, 
                                 dunnet = FALSE, trend = FALSE,
                                 iter_control = list(tol = 1e-5, max_iter = 20, 
                                                     verbose = FALSE),
                                 em_control = list(tol = 1e-5, max_iter = 100),
                                 lme_control = lme4::lmerControl(), 
                                 mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100), 
                                 trend_control = NULL)

#converts Phyloseq results to dataframe
Phyloseq_test_df=
  Phyloseq_test$res %>% 
  as.data.frame 

#makes list of genera that differed by either Species or Time
#This is done based on FDR-corrected p-values
#this is filtered so that only genus-level results are considered
ancom_diff_gen=
  Phyloseq_test_df %>%
  filter(diff_SpeciesT | diff_TimeJan21) %>%
  filter(str_detect(taxon, "Genus")) %>%
  mutate(Genus = str_remove(taxon, "Genus:"))

ancom_diff_gen_list=
  ancom_diff_gen$Genus

  
  
  
  
