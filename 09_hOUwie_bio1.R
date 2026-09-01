# ==============================================================================
# hOUwie
# ==============================================================================
# Fits hOUwie models using mean annual temperature (bio_1) as the continuous trait
# ==============================================================================
# Setup
#-------------------------------------------------------------------------------
rm(list = ls())

# Define focal climate variable (BIO1)
focal_var <- "mean_bio_1"
focal_bio <- "bio_1"

data_dir <- "/home/lenarh/data/bee_nesting_sociality"
curated_data_dir <- file.path(data_dir, "curated_data")
corhmm_results_dir <- file.path(
  data_dir,
  "results",
  "corHMM",
  "Aug 31 Run (corHMM fixed)"
)
results_dir <- file.path(
  data_dir,
  "results",
  "hOUwie",
  "Aug 31 Run (corHMM fixed)"
)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(corHMM)
library(OUwie)
library(parallel)

source("00_utility_functions.R")

# Reloading traits, tree and climatic data
traits <- read.csv(
  file.path(curated_data_dir, "bee_traits_clean.csv")
)

phy <- read.tree(
  file.path(curated_data_dir, "bee_tree_pruned.tre")
)

all_climatic_vars <- list.files(
  path = file.path(curated_data_dir, "climate_summaries"),
  pattern = "summstats\\.csv$",
  full.names = TRUE
)

# Keep only target climate variables: temp, precip (1, 12), temp and precip seasonality (4, 15)
all_climatic_vars <- all_climatic_vars[grep(
  paste(c("bio01_", "bio12_", "bio04_", "bio15_"), collapse = "|"),
  all_climatic_vars
)]

climatic_list <- lapply(all_climatic_vars, read.csv)

# Merge all climate variables by species name
merged_climatic_vars <- climatic_list[[1]] # initializes merged climate dataframe using first climate variable in list (mean_bio_1)

for (i in 2:length(climatic_list)) {
  # adds each remaining climate variable one by one (series of merges by species name)
  one_climatic_var <- climatic_list[[i]]
  merged_climatic_vars <- merge(
    merged_climatic_vars,
    one_climatic_var,
    by = "species"
  )
}

# Keep only the mean columns
merged_climatic_vars <- merged_climatic_vars[, c(
  1,
  grep("mean", colnames(merged_climatic_vars))
)]

# Merge with trait data
merged_traits <- merge(
  traits,
  merged_climatic_vars,
  by.x = "tips",
  by.y = "species"
)

colnames(merged_traits)[8:11] <- c(
  "mean_bio_1",
  "mean_bio_4",
  "mean_bio_12",
  "mean_bio_15"
)

#-------------------------------------------------------------------------------
# Log-transform continuous variables
#-------------------------------------------------------------------------------
merged_traits$mean_bio_1 <- log((merged_traits$mean_bio_1) + 273) # convert °C to Kelvin for temp
merged_traits$mean_bio_12 <- log(merged_traits$mean_bio_12)
merged_traits$mean_bio_15 <- log(merged_traits$mean_bio_15)
merged_traits$mean_bio_4 <- log(merged_traits$mean_bio_4)

# Remove rows with missing/non-finite values for focal climate variable
merged_traits <- subset(
  merged_traits,
  is.finite(merged_traits[[focal_var]])
)

# Prune phylogeny to match data
phy <- keep.tip(phy, which(phy$tip.label %in% merged_traits$tips))

#-------------------------------------------------------------------------------
# Align data and tree tips
#-------------------------------------------------------------------------------
dat <- merged_traits
shared_species <- intersect(dat$tips, phy$tip.label)

all(shared_species %in% dat$tips)
all(shared_species %in% phy$tip.label)

dat <- dat[match(shared_species, dat$tips), ]
phy <- keep.tip(phy, shared_species)
dat <- dat[match(phy$tip.label, dat$tips), ]

# Keep only discrete traits and the focal continuous trait
# Change the focal continuous trait depending on which you want to analyze
dat <- dat[, c("tips", "sociality_binary", "nest_binary", focal_var)]

phy <- keep.tip(phy, dat$tips)

#-------------------------------------------------------------------------------
# Load best-fitting corHMM model
#-------------------------------------------------------------------------------
load(
  file.path(
    corhmm_results_dir,
    "corHMM_dredge_binary_best_fit_index_matrix_Aug31.Rsave"
  )
)

corhmm_tbl <- read.csv(
  file.path(
    corhmm_results_dir,
    "corHMM_tbl_dredge_Aug31.csv"
  )
)

cid_disc_model <- index_matrix

#-------------------------------------------------------------------------------
# Define OU parameter structures for continuous trait models
#-------------------------------------------------------------------------------
# CID: character-independent model (same optima across discrete states)
cid_oum_model <- nest_oum_model <- soc_oum_model <- full_oum_model <-
  getOUParamStructure("OUM", 4, 2, null.model = TRUE) # (i.e. all observed states have the same optima)

# Custom optima for each model
# Sociality model: shared optima by social state
#     Estimates one optima for all social species and one for all solitary species,
#     assuming nesting strategy does not affect the optimum.
soc_oum_model[3, c(1, 3, 5, 7)] <- 3 # Assigns optimum "3" to all social states (1, 3, 5, 7)
soc_oum_model[3, c(2, 4, 6, 8)] <- 4 # Assigns optimum "4" to all solitary states (2, 4, 6, 8)

# Nesting model: shared optima by nest type
nest_oum_model[3, c(1, 2, 5, 6)] <- 3
nest_oum_model[3, c(3, 4, 7, 8)] <- 4

# Full model: separate optima for all trait combinations
full_oum_model[3, c(1:8)] <- c(3:6)

#-------------------------------------------------------------------------------
# Bundle models to run (BM1, OU1, and various OUM configurations)
#-------------------------------------------------------------------------------
model_list <- list(
  list(2, cid_disc_model, "BM1"),
  list(2, cid_disc_model, "OU1"),
  list(2, cid_disc_model, soc_oum_model),
  list(2, cid_disc_model, nest_oum_model),
  list(2, cid_disc_model, full_oum_model),
  list(2, cid_disc_model, cid_oum_model)
) # the 2 is for two rate classes

names(model_list) <- paste0(
  c(
    "bm1_8states_run_",
    "ou1_8states_run_",
    "oum_soc_8states_run_",
    "oum_nest_8states_run_",
    "oum_full_8states_run_",
    "oum_cid_8states_run_"
  ),
  focal_bio
)
#-------------------------------------------------------------------------------
# Wrapper function to run and save each model
#-------------------------------------------------------------------------------
quickFunc <- function(model_list, model_name) {
  res <- hOUwie(
    phy,
    dat,
    2,
    model_list[[2]],
    model_list[[3]],
    nSim = 50,
    diagn_msg = TRUE,
    # adaptive_sampling = FALSE,
    n_starts = 10,
    ncores = 10
  )
  file.name <- file.path(results_dir, paste0(model_name, ".Rsave"))
  save(res, file = file.name)
}

#-------------------------------------------------------------------------------
# Run all models in parallel for focal climatic variable
#-------------------------------------------------------------------------------
mclapply(
  1:6,
  function(x) quickFunc(model_list[[x]], names(model_list)[x]),
  mc.cores = 100
)
