#This script prepares a description of the fungi found on the
#control plants through LSU sequencing

#make sure the data import + initial processing 
#has been done prior to running this script
source("LSU_DataImport.R")

TopOTUs_Control <- otu_mat_rare_rel %>% 
  t() %>%
  as.data.frame() %>%
  mutate(sample = rownames(.)) %>%
  left_join(select(samples_df, sample, Species, Stand, Time), by = "sample") %>%
  filter(Stand == "C") %>%
  group_by(Species, Time) %>%
  summarise(across(
    .cols = where(is.numeric),
    .fns = mean,
    .names = "{.col}"
  ), .groups = "drop") %>%
  pivot_longer(
    cols = -c(Species, Time),
    names_to = "otu",
    values_to = "Abundance"
  ) %>%
  mutate(SpeciesTime = paste(Species, Time, sep = "_")) %>%
  select(otu, SpeciesTime, Abundance) %>%
  pivot_wider(
    names_from = SpeciesTime,
    values_from = Abundance
  ) %>%
  mutate(GlobalMean=(A_June20+T_June20+A_Jan21+T_Jan21)/4) %>%
  left_join(tax_mat)




TopOTUs_Control <- otu_mat_rare_rel %>% 
  t() %>%
  as.data.frame() %>%
  mutate(sample = rownames(.)) %>%
  left_join(select(samples_df, sample, Species, Stand, Time), by = "sample") %>%
  filter(Stand == "C") %>% select(Species, Time) %>% View
  
  # Reshape to long format to compute both mean and non-zero count
  pivot_longer(
    cols = where(is.numeric),
    names_to = "otu",
    values_to = "Abundance"
  ) %>%
  group_by(Species, Time, otu) %>%
  summarise(
    Abundance = mean(Abundance),
    NonZeroCount = sum(Abundance > 0),
    .groups = "drop"
  ) %>%
  mutate(SpeciesTime = paste(Species, Time, sep = "_")) %>%
  
  # Pivot means to wide format
  select(otu, SpeciesTime, Abundance) %>%
  pivot_wider(
    names_from = SpeciesTime,
    values_from = Abundance,
    names_prefix = "mean_"
  ) %>%
  
  # Add non-zero count columns
  left_join(
    otu_mat_rare_rel %>%
      t() %>%
      as.data.frame() %>%
      mutate(sample = rownames(.)) %>%
      left_join(select(samples_df, sample, Species, Stand, Time), by = "sample") %>%
      filter(Stand == "C") %>%
      pivot_longer(
        cols = where(is.numeric),
        names_to = "otu",
        values_to = "Abundance"
      ) %>%
      group_by(Species, Time, otu) %>%
      summarise(
        NonZeroCount = sum(Abundance > 0),
        .groups = "drop"
      ) %>%
      mutate(SpeciesTime = paste(Species, Time, sep = "_")) %>%
      select(otu, SpeciesTime, NonZeroCount) %>%
      pivot_wider(
        names_from = SpeciesTime,
        values_from = NonZeroCount,
        names_prefix = "count_"
      ),
    by = "otu"
  ) %>%
  
  # Add global mean of means
  mutate(GlobalMean = rowMeans(select(., starts_with("mean_")), na.rm = TRUE)) %>%
  
  # Join taxonomy info
  left_join(tax_mat, by = "otu")
