# ==============================================================================
# corHMM SIMMAPs
# ==============================================================================
# This script loads the best-fit corHMM model for binary sociality and
# nesting traits in bees, generates stochastic character maps (SIMMAPs),
# summarizes evolutionary transition frequencies across simulations,
# and exports transition summary statistics.
# ==============================================================================
# Setup
#-------------------------------------------------------------------------------
rm(list = ls())

data_dir <- "/home/lenarh/data/bee_nesting_sociality"

curated_data_dir <- file.path(data_dir, "curated_data")
results_dir <- file.path(data_dir, "results", "corHMM")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(corHMM)
# devtools::install_github("thej022214/corHMM")

#-------------------------------------------------------------------------------
# Organizing dataset
#-------------------------------------------------------------------------------
# Reloading traits, tree and climatic data
traits <- read.csv(
  file.path(curated_data_dir, "bee_traits_clean.csv")
)

phy <- read.tree(
  file.path(curated_data_dir, "bee_tree_pruned.tre")
)

# Drop tips in the tree that don’t match any in the trait data
phy <- keep.tip(phy, which(phy$tip.label %in% traits$tips))

# Find species shared by both datasets
dat <- traits
shared_species <- intersect(dat$tips, phy$tip.label)

all(shared_species %in% dat$tips)
all(shared_species %in% phy$tip.label)

# Reorder the data and tree to align species
dat <- dat[match(shared_species, dat$tips), ]
phy <- keep.tip(phy, shared_species)
dat <- dat[match(phy$tip.label, dat$tips), ]

# Keep relevant trait columns
dat <- dat[, c("tips", "sociality_binary", "nest_binary")]

#-------------------------------------------------------------------------------
# Load corHMM dredge results
#-------------------------------------------------------------------------------
# load(file.path(results_dir, "corHMM_dredge_binary_Aug4.Rsave"))
load(file.path(results_dir, "corhmm_dredge_binary_june8.Rsave"))

# model <- model_fits[[35]]$solution # best supported model
# root.p <- model_fits[[35]]$root.p
model <- dredge_sociality[[28]]$solution # best supported model
root.p <- dredge_sociality[[28]]$root.p

corHMM_results <- list(data, model, root.p)

# save(
#   corHMM_results,
#   file = file.path(results_dir, "corHMM_dredge_results_for_simmap_Aug4.Rsave")
# )
save(
  corHMM_results,
  file = file.path(results_dir, "corHMM_dredge_results_for_simmap_Aug5.Rsave")
)


# load(file.path(results_dir, "corhmm_dredge_results_for_simmap_june8.Rsave"))

#-------------------------------------------------------------------------------
# Generate SIMMAPs
#-------------------------------------------------------------------------------
simmaps <- makeSimmap(
  tree = phy,
  data = dat,
  model = model,
  rate.cat = 2,
  nSim = 100,
  nCores = 5,
  root.p = root.p
)

# save(simmaps, file = file.path(results_dir, "simmaps_Aug4.Rsave"))
save(simmaps, file = file.path(results_dir, "simmaps_Aug5.Rsave"))

#-------------------------------------------------------------------------------
# Load previously generated SIMMAPs
#-------------------------------------------------------------------------------
# load(
#   file.path(
#     results_dir,
#     "simmaps_Aug4.Rsave"
#   )
# )
load(
  file.path(
    results_dir,
    "simmaps_Aug5.Rsave"
  )
)

#-------------------------------------------------------------------------------
# Summarize transitions
#-------------------------------------------------------------------------------
simmap_summaries <- lapply(simmaps, summarize_single_simmap)
summary_df <- summarize_transition_stats(simmap_summaries)

print(summary_df)

# write.csv(
#   summary_df,
#   file = file.path(results_dir, "corHMM_transitions_summary_Aug4.csv"),
#   row.names = FALSE
# )
write.csv(
  summary_df,
  file = file.path(results_dir, "corHMM_transitions_summary_Aug5.csv"),
  row.names = FALSE
)

# Export as PDF
# pdf(
#   file.path(results_dir, "corHMM_transition_summary_Aug4.pdf"),
#   width = 15,
#   height = 8
# )
pdf(
  file.path(results_dir, "corHMM_transition_summary_Aug5.pdf"),
  width = 15,
  height = 8
)

plot_transition_summary(simmap_summaries)

dev.off()
