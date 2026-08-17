#===============================================================================
# phylANOVA
#===============================================================================
# Run phylogenetic ANOVA (phylANOVA) to test associations between
# bee traits (sociality, nesting) and climatic variables.
#===============================================================================
# Setup
#-------------------------------------------------------------------------------
rm(list = ls())

data_dir <- "/home/lenarh/data/bee_nesting_sociality"
curated_data_dir <- file.path(data_dir, "curated_data")
results_dir <- file.path(data_dir, "results", "phylANOVA")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(phytools)

#-------------------------------------------------------------------------------
# Load trait data, tree, and list of climate summary statistic files
#-------------------------------------------------------------------------------
traits <- read.csv(file.path(curated_data_dir, "bee_traits_clean.csv"))

tree <- read.tree(file.path(curated_data_dir, "bee_tree_pruned.tre"))

all_climatic_vars <- list.files(
  file.path(curated_data_dir, "climate_summaries"),
  pattern = "summstats.csv",
  full.names = TRUE
)

head(all_climatic_vars)

#-------------------------------------------------------------------------------
# Run phylANOVA between sociality_binary and each climate variable
#-------------------------------------------------------------------------------
sink(file.path(results_dir, "phylANOVA_sociality_results.txt"))

for (climate_index in 1:length(all_climatic_vars)) {
  # Load and clean climate data
  climate <- read.csv(all_climatic_vars[climate_index])

  # Identify mean climate column
  mean_col <- grep("^mean_", colnames(climate), value = TRUE)

  # Remove species with NA climate means
  climate <- subset(climate, !is.na(climate[, mean_col]))

  # Identify species with complete data in all datasets
  sampled_species <- intersect(
    intersect(traits$tips, climate$species),
    tree$tip.label
  )

  # Subset datasets to include only sampled species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(
    tree,
    tree$tip.label[tree$tip.label %in% sampled_species]
  )

  # Merge traits and climate data by species
  merged_table <- merge(
    subset_traits,
    subset_climate,
    by.x = "tips",
    by.y = "species"
  )

  # Prepare named vectors for phylANOVA
  sociality <- merged_table$sociality_binary
  names(sociality) <- merged_table$tips

  one_clim_var <- merged_table[, mean_col]
  names(one_clim_var) <- merged_table$tips

  # Print description of analysis
  label <- gsub(
    "_climate_summstats.csv",
    "",
    basename(all_climatic_vars[climate_index])
  )
  print(paste0("sociality ~ ", label))

  # Run phylANOVA
  results <- phylANOVA(
    subset_tree,
    sociality,
    one_clim_var,
    p.adj = "bonferroni"
  )
  print(results)
}

sink()

#-------------------------------------------------------------------------------
# Run phylANOVA between nesting_binary and each climate variable
#-------------------------------------------------------------------------------
sink(file.path(results_dir, "phylANOVA_nesting_results.txt"))

for (climate_index in 1:length(all_climatic_vars)) {
  climate <- read.csv(all_climatic_vars[climate_index])

  mean_col <- grep("^mean_", colnames(climate), value = TRUE)

  climate <- subset(climate, !is.na(climate[, mean_col]))

  sampled_species <- intersect(
    intersect(traits$tips, climate$species),
    tree$tip.label
  )

  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(
    tree,
    tree$tip.label[tree$tip.label %in% sampled_species]
  )

  merged_table <- merge(
    subset_traits,
    subset_climate,
    by.x = "tips",
    by.y = "species"
  )

  boxplot(
    merged_table[, mean_col] ~ merged_table$nest_binary,
    xlab = "nest",
    ylab = "env var"
  )

  title(gsub(
    "_climate_summstats.csv",
    "",
    basename(all_climatic_vars[climate_index])
  ))

  nests <- merged_table$nest_binary
  names(nests) <- merged_table$tips

  one_clim_var <- merged_table[, mean_col]
  names(one_clim_var) <- merged_table$tips

  label <- gsub(
    "_climate_summstats.csv",
    "",
    basename(all_climatic_vars[climate_index])
  )

  print(paste0("nesting type ~ ", label))

  results <- phylANOVA(subset_tree, nests, one_clim_var, p.adj = "bonferroni")
  print(results)
}

sink()
dev.off()

#-------------------------------------------------------------------------------
# Run phylANOVA using 4-level combination of sociality and nesting traits
#-------------------------------------------------------------------------------
sink(file.path(results_dir, "phylANOVA_combined_results.txt"))

traits$comb_nest_soc <- paste(
  traits$sociality_binary,
  traits$nest_binary,
  sep = "_"
)
table(traits$comb_nest_soc)

for (climate_index in 1:length(all_climatic_vars)) {
  climate <- read.csv(all_climatic_vars[climate_index])

  mean_col <- grep("^mean_", colnames(climate), value = TRUE)

  climate <- subset(climate, !is.na(climate[, mean_col]))

  sampled_species <- intersect(
    intersect(traits$tips, climate$species),
    tree$tip.label
  )

  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(
    tree,
    tree$tip.label[tree$tip.label %in% sampled_species]
  )

  merged_table <- merge(
    subset_traits,
    subset_climate,
    by.x = "tips",
    by.y = "species"
  )

  nests <- merged_table$comb_nest_soc
  names(nests) <- merged_table$tips

  one_clim_var <- merged_table[, mean_col]
  names(one_clim_var) <- merged_table$tips

  label <- gsub(
    "_climate_summstats.csv",
    "",
    basename(all_climatic_vars[climate_index])
  )

  print(paste0("nesting type ~ ", label))

  results <- phylANOVA(subset_tree, nests, one_clim_var, p.adj = "bonferroni")
  print(results)
}

sink()

#-------------------------------------------------------------------------------
# Done!
#-------------------------------------------------------------------------------
