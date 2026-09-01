# ==============================================================================
# Ancestral state reconstruction of sociality and nesting strategy
# ==============================================================================
# Uses corHMM to model the correlated evolution of sociality and nesting strategy
# in bees under hidden rate models, selects the best-fitting model via AIC, and
# reconstructs ancestral states across the phylogeny.
# ==============================================================================
# Setup
#-------------------------------------------------------------------------------

# Note to self: Aug 31 run is AFTER re-downloading fixed corHMM; previous runs used a buggy version of corHMM
# Aug 31 run is the one to use in publication.

rm(list = ls())

repo_dir <- "/home/lenarh/repos/bee_nesting_sociality"
data_dir <- "/home/lenarh/data/bee_nesting_sociality"

curated_data_dir <- file.path(data_dir, "curated_data")
results_dir <- file.path(data_dir, "results", "corHMM")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(corHMM)
library(phytools)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(scales)
# devtools::install_github("thej022214/corHMM")

# Load utility functions
source(file.path(repo_dir, "00_utility_functions.R"))

# Reloading traits and tree
traits <- read.csv(
  file.path(curated_data_dir, "bee_traits_clean.csv")
)

phy <- read.tree(
  file.path(curated_data_dir, "bee_tree_pruned.tre")
)

#-------------------------------------------------------------------------------
# Preparing dataset
#-------------------------------------------------------------------------------
# Reorder data to match tree
dat <- traits[match(phy$tip.label, traits$tips), ]

# Keep relevant trait columns
dat <- dat[, c("tips", "sociality_binary", "nest_binary")]

head(dat)
nrow(dat) # 4214

setequal(phy$tip.label, traits$tips) # TRUE (meaning tree/traits contain same species)
all(dat$tips == phy$tip.label) # TRUE (meaning tips in dat are in same order as tips in tree)

#-------------------------------------------------------------------------------
# Model selection (run or load corHMMdredge) with root fixed as solitary/ground
#-------------------------------------------------------------------------------
sort(unique(paste(dat$sociality_binary, dat$nest_binary, sep = "_"))) # check which state solitary/ground is for root
# [1] "social_aboveground"   "social_ground"
# [3] "solitary_aboveground" "solitary_ground"

# Run model selection with corHMMdredge (if not already done)
model_fits <- corHMM::corHMMDredge(
  phy,
  dat,
  max.rate.cat = 2,
  root.p = c(0, 0, 0, 1),
  n.cores = 8,
  nstarts = 5,
  use_RTMB = TRUE,
  seed = 1234 # set seed for reproducibility
)

save(
  model_fits,
  file = file.path(results_dir, "corHMM_dredge_binary_Aug31.Rsave")
)

# load(file.path(results_dir, corHMM_dredge_binary_Aug31.Rsave"))

# Make model comparison table
corhmm_model_comparison <- corHMM:::getModelTable(model_fits)

# # Make model comparison table
# corhmm_model_comparison <- corHMM:::getModelTable(dredge_sociality)

write.csv(
  corhmm_model_comparison,
  file = file.path(results_dir, "corHMM_tbl_dredge_Aug31.csv")
)
# write.csv(
#   corhmm_model_comparison,
#   file = file.path(results_dir, "corHMM_tbl_dredge_Aug5.csv")
# )

corhmm_model_comparison

# Best supported model (7/27/26): Model 28
# np (number of parameters): 10
# nRateCat: 2 (two-rate-class hidden rate model)
# Fit:
# lnLik     AIC
# -503.4647 1026.929

# Best supported model (8/4/26): Model 35
# np (number of parameters): 9
# nRateCat: 2 (two-rate-class hidden rate model)
# Fit:
# lnLik     AIC
# -500.87   1019.735

# Extract the transition rate matrix from the best model
rates_mat <- model_fits[[which.min(corHMM:::getModelTable(model_fits)$AIC)]]
# rates_mat <- dredge_sociality[[which.min(
#   corHMM:::getModelTable(dredge_sociality)$AIC
# )]]

write.csv(
  rates_mat$solution,
  file = file.path(results_dir, "corHMM_transition_rates_Aug31.csv"),
  row.names = TRUE
)
# write.csv(
#   rates_mat$solution,
#   file = file.path(results_dir, "corHMM_transition_rates_Aug5.csv"),
#   row.names = TRUE
# )

rates_mat

#===============================================================================
# Ancestral State Reconstruction & Plotting
#===============================================================================
#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
WtQ <- function(Q, Weights) {
  # Weighted average of transition rates using AIC weights
  Weights <- Weights[!is.na(Q)] / sum(Weights[!is.na(Q)])
  AvgQ <- sum(Q[!is.na(Q)] * Weights)
  return(AvgQ)
}

getModelAvgRate <- function(file) {
  # Loads model set, computes AIC weights, and creates a model-averaged rate matrix
  load(file)
  rate.mat <- obj$res_ARD.ARD.2$index.mat
  AICc <- unlist(lapply(obj, function(x) x$AICc))
  AICwt <- exp(-0.5 * AICc - min(AICc)) / sum(exp(-0.5 * AICc - min(AICc)))
  # Solutions <- lapply(obj, function(x) x$solution)
  Solutions <- lapply(obj, function(x) c(na.omit(c(x$solution))))
  Solutions[[1]] <- c(
    Solutions[[1]][1],
    NA,
    Solutions[[1]][2],
    NA,
    NA,
    Solutions[[1]][1],
    NA,
    Solutions[[1]][2]
  )
  Solutions[[2]] <- c(
    Solutions[[2]][1],
    NA,
    Solutions[[2]][2],
    NA,
    NA,
    Solutions[[2]][1],
    NA,
    Solutions[[2]][2]
  )
  Rates <- do.call(rbind, Solutions)
  p.wt <- apply(Rates, 2, function(x) WtQ(x, AICwt))
  rate.mat[!is.na(rate.mat)] <- p.wt
  return(rate.mat)
}

getTipRecon <- function(file) {
  # Reconstructs states using model-averaged transition matrix, returns marginal likelihoods
  load(file)
  phy <- obj$res_ER$phy
  data <- obj$res_ER$data
  root.p <- obj$res_ARD$root.p
  index.mat <- obj$res_ARD.ARD.2$index.mat
  p <- getModelAvgRate(file)[sapply(
    1:max(index.mat, na.rm = TRUE),
    function(x) match(x, index.mat)
  )]
  res <- corHMM(
    phy = phy,
    data = data,
    rate.cat = 2,
    rate.mat = index.mat,
    node.states = "marginal",
    p = p,
    root.p = root.p,
    get.tip.states = TRUE
  )
  return(res)
}

# Plots ancestral state reconstructions with colored pie charts at nodes
plotRECON <- function(
  phy,
  likelihoods,
  piecolors = NULL,
  cex = 0.5,
  pie.cex = 0.25,
  file = NULL,
  height = 11,
  width = 8.5,
  show.tip.label = TRUE,
  title = NULL,
  ...
) {
  if (is.null(piecolors)) {
    piecolors <- rev(c(
      "#00204DFF",
      "#575C6DFF",
      "#A69D75FF",
      "#FFEA46FF",
      "#00204DFF",
      "#575C6DFF",
      "#A69D75FF",
      "#FFEA46FF"
    ))
  }
  if (!is.null(file)) {
    pdf(file, height = height, width = width, useDingbats = FALSE)
  }
  plot(phy, cex = cex, show.tip.label = show.tip.label, ...)

  if (!is.null(title)) {
    title(main = title)
  }
  nodelabels(pie = likelihoods, piecol = piecolors, cex = pie.cex)
  states <- colnames(likelihoods)
  legend(
    x = "topleft",
    states,
    cex = 0.8,
    pt.bg = piecolors,
    col = "black",
    pch = 21
  )

  if (!is.null(file)) {
    dev.off()
  }
}

#-------------------------------------------------------------------------------
# Plot ancestral state reconstruction
#-------------------------------------------------------------------------------
# Extract likelihoods for internal nodes from best model
mod_table <- corHMM:::getModelTable(model_fits)
best_fit <- model_fits[[which.min(mod_table$dAIC)]]
index_matrix <- best_fit$index.mat
save(
  index_matrix,
  file = file.path(
    results_dir,
    "corHMM_dredge_binary_best_fit_index_matrix_Aug31.Rsave"
  )
)
# mod_table <- corHMM:::getModelTable(dredge_sociality)
# best_fit <- dredge_sociality[[which.min(mod_table$dAIC)]]
# index_matrix <- best_fit$index.mat
# save(
#   index_matrix,
#   file = file.path(
#     results_dir,
#     "corHMM_dredge_binary_best_fit_index_matrix_Aug5.Rsave"
#   )
# )

p <- corHMM:::MatrixToPars(best_fit)
anc_recon <- ancRECON(
  phy = phy,
  data = best_fit$data,
  p = p,
  method = "marginal",
  rate.cat = best_fit$rate.cat,
  ntraits = NULL,
  rate.mat = best_fit$index.mat,
  root.p = best_fit$root.p
)

anc_recon$lik.anc.states[1, ] # checking root state
