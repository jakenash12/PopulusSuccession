library(tidyverse)
library(fungarium)
library(vegan)

#Sets path for input files
Input_path="./InputFiles/"

#Reads in OTU table
otu_mat_raw <- 
  read.delim(paste(Input_path, "LSU_dada2_table.tsv", sep=""), skip=1) %>%
  rename_with(~ sub(".*bc", "bc", .x)) %>%
  rename("otu"="X.OTU.ID") %>%
  column_to_rownames("otu")

#reads in mapping file
samples_df <- read.delim(paste(Input_path, "PAMB1_metadata.txt",sep=""))  %>%
  mutate(sample = sub(".*bc", "bc", sample)) %>%
  mutate(Time=factor(Time, levels=c("June20", "Jan21")))

#generates list of samples to be filtered from 
#various analyses (uninoculated controls/negs)
FilterList=
  filter(samples_df,Stand!="C", Species!="Neg")$sample

#reads in custom taxonomy edits file
TaxEdits=read.delim(paste(Input_path,"TaxonomyEdits.txt", sep=""))

#reads in BLAST-NCBI taxonomy and formats the taxonomy 
#string into different columns using a helper function 
extract_text <- function(x) {
  sub(".+__", "", x)
}

tax_mat <- read.delim(paste(Input_path, "LSU_dada2_repseqs_BLAST_taxonomy.tsv", sep="")) %>%
  rename("otu"="Feature.ID") %>%
  left_join(TaxEdits) %>%
  mutate(Taxon=case_when(
    !is.na(EditedTaxonomy)~EditedTaxonomy,
    .default=Taxon
  )) %>%
  select(-EditedTaxonomy) %>%
  separate_wider_delim(too_few="align_start", cols = Taxon, delim = ";", names = c("Kingdom", "Phylum","Class", "Order","Family","Genus","Species")) %>%
  mutate(across(!contains("otu"), ~extract_text(.))) %>%
  mutate(Species=case_when(
    !is.na(Genus) & !is.na(Species) ~ paste(Genus, Species, sep="_"),
    .default=Species)
  ) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Kingdom = gsub("Unassigned", NA, Kingdom))

#generates list of fungal ASVs to filter table
#(excludes plant, nematode, etc)
fungal_ASVs=
  filter(tax_mat, Kingdom=="Fungi")$otu

otu_mat =
  otu_mat_raw %>%
  filter(rownames(.) %in% fungal_ASVs)

#reads in custom dark septate endophyte list
DSE_List=read.delim(paste(Input_path, "DSE_List.txt", sep=""))

#assigns guild with funguild
tax_mat_fg=
  tax_mat %>%
  as.data.frame %>%
  fungarium::fg_assign(., tax_cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) %>%
  as.data.frame %>%
  mutate(primary_guild = if_else(grepl("\\|", guild), 
                                 sub(".*\\|(.*?)\\|.*", "\\1", guild), 
                                 NA_character_)) %>%
  mutate(primary_guild=case_when(
    Genus %in% DSE_List$Genus ~ "DSE",
    .default=primary_guild
  ))

#calculates per sample sequencing depth
seq_depth=
  otu_mat %>%
  colSums %>%
  as.data.frame

#decided on a rarefaction depth of 816
#samples with less seqs than that are discarded
set.seed(67)
raredepth=816
otu_mat_rare=
  otu_mat %>%
  t %>%
  as.data.frame %>%
  rrarefy(raredepth) %>%
  as.data.frame %>%
  filter(rowSums(.) == raredepth) %>%
  t %>%
  as.data.frame
  
#relativizes rarefied dataframe
otu_mat_rare_rel=
  otu_mat_rare %>%
  mutate(across(everything(), ~ . / sum(.)))

#creates transposed dataframe that is needed
#for some functions in vegan
otu_mat_rare_t=as.data.frame(t(otu_mat_rare))

