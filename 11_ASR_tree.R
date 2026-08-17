# ============================================================================
# Circular Phylogeny: ASR Branch Coloring + Trait Ring + Family Labels
# ============================================================================
# Circular phylogeny with branches colored by ancestral state reconstruction (ASR),
# a trait ring at the tips, time axis, and labeled families with MRCA nodes denoted.
# ==============================================================================
# Setup
# ------------------------------------------------------------------------------
rm(list = ls())

data_dir <- "/home/lenarh/data/bee_nesting_sociality"
results_dir <- file.path(data_dir, "results", "Phylogeny")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(ape)
library(phytools)
library(viridis)
library(dplyr)

traits <- read.csv(
  file.path(data_dir, "curated_data", "bee_traits_clean.csv")
)
phy <- read.tree(
  file.path(data_dir, "curated_data", "bee_tree_pruned.tre")
)
Ntip(phy) # 4214
nrow(traits) # 4214

# Match trait rows to tree tips
traits <- traits[match(phy$tip.label, traits$tips), ]

# Create combined trait categories
traits$trait_combo <- paste(
  traits$sociality_binary,
  traits$nest_binary,
  sep = "/"
)

# Rename trait categories for readability
traits$trait_combo <- recode(
  traits$trait_combo,
  "solitary/ground" = "Solitary/Ground",
  "solitary/aboveground" = "Solitary/Above-ground",
  "social/ground" = "Social/Ground",
  "social/aboveground" = "Social/Above-ground"
)

# Create named vector for plotting
trait_combo_vector <- setNames(traits$trait_combo, traits$tips)

# Set colors for each trait combo
combo_colors <- viridis::viridis(n = 4, option = "cividis")
names(combo_colors) <- c(
  "Solitary/Ground",
  "Solitary/Above-ground",
  "Social/Ground",
  "Social/Above-ground"
)

# ------------------------------------------------------------------------------
# Load corHMM results and reconstruct ancestral states
# ------------------------------------------------------------------------------

load(
  file.path(
    data_dir,
    "results",
    "corHMM",
    "corhmm_dredge_binary_june8.Rsave"
  )
)

# Find the best-supported model
mod_table <- corHMM:::getModelTable(dredge_sociality)
best_model_number <- which.min(mod_table$AIC)
best_fit <- dredge_sociality[[best_model_number]]

message(
  "Using model ",
  best_model_number,
  " with AIC = ",
  mod_table$AIC[best_model_number]
)

# Run marginal ancestral-state reconstruction
p <- corHMM:::MatrixToPars(best_fit)

anc_recon_result <- corHMM::ancRECON(
  phy = phy,
  data = best_fit$data,
  p = p,
  method = "marginal",
  rate.cat = best_fit$rate.cat,
  ntraits = NULL,
  rate.mat = best_fit$index.mat,
  root.p = best_fit$root.p
)

# Matrix of marginal state probabilities for internal nodes
anc_recon <- anc_recon_result$lik.anc.states

# Check that the ASR matrix matches the tree
stopifnot(
  is.matrix(anc_recon),
  nrow(anc_recon) == phy$Nnode
)

# Check dimensions
print(dim(anc_recon))

# The ancRECON output lacks column names, so get state names
# from the fitted transition-rate matrix
state_names <- rownames(best_fit$solution)

print(state_names)
stopifnot(length(state_names) == ncol(anc_recon))

# Assign state names to the ASR probability matrix
colnames(anc_recon) <- state_names

print(colnames(anc_recon))

asr_state_key <- c(
  "R1 social|aboveground" = "Social/Above-ground",
  "R1 solitary|aboveground" = "Solitary/Above-ground",
  "R1 social|ground" = "Social/Ground",
  "R1 solitary|ground" = "Solitary/Ground",
  "R2 social|aboveground" = "Social/Above-ground",
  "R2 solitary|aboveground" = "Solitary/Above-ground",
  "R2 social|ground" = "Social/Ground",
  "R2 solitary|ground" = "Solitary/Ground"
)

# Confirm that all reconstructed-state names occur in the lookup key
stopifnot(all(colnames(anc_recon) %in% names(asr_state_key)))

n_tips <- Ntip(phy)
branch_colors <- rep("gray80", nrow(phy$edge))

# Return one character value for every internal node
node_states <- apply(anc_recon, 1, function(x) {
  if (max(x, na.rm = TRUE) > 0.5) {
    colnames(anc_recon)[which.max(x)]
  } else {
    NA_character_
  }
})

# Ensure node_states is a character vector, not a list
node_states <- as.character(node_states)

table(node_states, useNA = "ifany")

# # ------------------------------------------------------------------------------
# # Load ASR from corHMM
# # ------------------------------------------------------------------------------
# load(
#   file.path(
#     data_dir,
#     "results",
#     "corHMM",
#     # "corHMM_dredge_binary_best_fit_index_matrix_Aug4.Rsave"
#     "corHMM_dredge_binary_Aug4.Rsave"
#   )
# )

# anc_recon <- model_fits[[35]]$states

# asr_state_key <- c(
#   "R1 social|aboveground" = "Social/Above-ground",
#   "R1 solitary|aboveground" = "Solitary/Above-ground",
#   "R1 social|ground" = "Social/Ground",
#   "R1 solitary|ground" = "Solitary/Ground",
#   "R2 social|aboveground" = "Social/Above-ground",
#   "R2 solitary|aboveground" = "Solitary/Above-ground",
#   "R2 social|ground" = "Social/Ground",
#   "R2 solitary|ground" = "Solitary/Ground"
# )

# n_tips <- length(phy$tip.label)

# branch_colors <- rep("gray80", nrow(phy$edge))

# node_states <- apply(anc_recon, 1, function(x) {
#   if (max(x) > 0.5) names(which.max(x)) else NA
# })

for (i in 1:nrow(phy$edge)) {
  parent_node <- phy$edge[i, 1]
  if (parent_node > n_tips) {
    state <- node_states[parent_node - n_tips]
    if (!is.na(state)) {
      trait_label <- asr_state_key[state]
      branch_colors[i] <- combo_colors[trait_label]
    }
  }
}

# ------------------------------------------------------------------------------
# Functions: Coordinate conversion and axis drawing
# ------------------------------------------------------------------------------
toCart <- function(r, th, deg = FALSE) {
  if (deg) {
    th <- th * pi / 180
  }
  list(x = r * cos(th), y = r * sin(th))
}
toPolar <- function(x, y) {
  list(r = sqrt(x^2 + y^2), th = atan2(y, x))
}

add_time_axis_radial <- function(
  phy,
  ring_times,
  angle_deg = 0,
  tick_length = 0.015,
  label_cex = 0.4
) {
  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  n_tips <- length(phy$tip.label)
  root_radius <- max(sqrt(obj$xx^2 + obj$yy^2)[(n_tips + 1):length(obj$xx)])
  phy_height <- max(nodeHeights(phy))
  reversed_radii <- (phy_height - ring_times) / phy_height * root_radius
  theta <- angle_deg * pi / 180
  dx <- cos(theta)
  dy <- sin(theta)
  segments(0, 0, root_radius * dx, root_radius * dy, col = "black", lwd = 1)

  for (i in seq_along(reversed_radii)) {
    r <- reversed_radii[i]
    x_tick <- r * dx
    y_tick <- r * dy
    perp_dx <- -dy
    perp_dy <- dx

    # Draw tick marks
    tick_x0 <- x_tick - tick_length * root_radius * perp_dx
    tick_y0 <- y_tick - tick_length * root_radius * perp_dy
    tick_x1 <- x_tick + tick_length * root_radius * perp_dx
    tick_y1 <- y_tick + tick_length * root_radius * perp_dy
    segments(tick_x0, tick_y0, tick_x1, tick_y1, col = "black", lwd = 1)

    # Skip 100 label
    if (ring_times[i] != 100) {
      # Position label above tick
      label_offset <- 2 * tick_length * root_radius
      label_x <- x_tick + label_offset * perp_dx
      label_y <- y_tick + label_offset * perp_dy

      label_text <- if (ring_times[i] == 0) {
        "0 Ma"
      } else {
        as.character(ring_times[i])
      }
      text(
        label_x,
        label_y,
        labels = label_text,
        cex = label_cex,
        adj = c(0.5, 0)
      )
    }
  }
}

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------
pdf(
  file.path(results_dir, "ASR_phylogeny_labeled_intervals_Aug5.pdf"),
  width = 8,
  height = 8
)

par(mar = c(2, 2, 2, 2), xpd = TRUE)

# Plot phy
plot(
  phy,
  type = "fan",
  edge.color = branch_colors,
  edge.width = 0.5,
  show.tip.label = FALSE,
  no.margin = FALSE,
  open.angle = 15,
  rotate.phy = 7.5
)

# Add time axis
ring_times <- c(0, 20, 40, 60, 80, 100)
add_time_axis_radial(phy, ring_times, angle_deg = 0)

# ------------------------------------------------------------------------------
# Trait ring at tips
# ------------------------------------------------------------------------------
obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tips_xx <- obj$xx[1:n_tips]
tips_yy <- obj$yy[1:n_tips]
ri <- 0.2
len <- 15
space <- 4

for (i in 0:(len + 2 * space)) {
  p <- toPolar(tips_xx, tips_yy)
  p$r <- p$r + ri * i
  c <- toCart(p$r, p$th)
  if (i < space || i >= space + len) {
    points(c$x, c$y, col = "white", pch = 16, cex = 0.15)
  } else {
    trait_colors <- combo_colors[trait_combo_vector[phy$tip.label]]
    points(c$x, c$y, col = trait_colors, pch = 16, cex = 0.25)
  }
}

# # ------------------------------------------------------------------------------
# # Tip labels (ALL species)
# # ------------------------------------------------------------------------------
#
# obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
#
# n_tips <- length(phy$tip.label)
# tips_xx <- obj$xx[1:n_tips]
# tips_yy <- obj$yy[1:n_tips]
#
# tip_idx <- 1:n_tips  # ALL tips
#
# r_tip <- sqrt(tips_xx^2 + tips_yy^2)
# r_out <- max(r_tip) + ri * (space + len) + 0.02 * max(r_tip)
#
# theta <- atan2(tips_yy[tip_idx], tips_xx[tip_idx])
# x_lab <- r_out * cos(theta)
# y_lab <- r_out * sin(theta)
#
# rot <- theta * 180 / pi
# rot_flip <- ifelse(rot > 90 | rot < -90, rot + 180, rot)
#
# # Left/right side positioning (optional but helps readability a bit)
# pos_vec <- ifelse(cos(theta) >= 0, 4, 2)  # 4=right, 2=left
#
# for (k in seq_along(tip_idx)) {
#   text(
#     x = x_lab[k],
#     y = y_lab[k],
#     labels = phy$tip.label[tip_idx[k]],
#     cex = 0.25,
#     srt = rot_flip[k],   # scalar per call
#     pos = pos_vec[k],
#     offset = 0
#   )
# }

# ------------------------------------------------------------------------------
# Tip labels at intervals
# ------------------------------------------------------------------------------
obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)

n_tips <- length(phy$tip.label)
tips_xx <- obj$xx[1:n_tips]
tips_yy <- obj$yy[1:n_tips]

label_every <- 5 # try 15, 20, 30, 40
tip_idx <- seq(1, n_tips, by = label_every)

r_tip <- sqrt(tips_xx^2 + tips_yy^2)

# push labels outside the trait ring; increase if still crowded
r_out <- max(r_tip) + ri * (space + len) + 0.03 * max(r_tip)

theta <- atan2(tips_yy[tip_idx], tips_xx[tip_idx])
x_lab <- r_out * cos(theta)
y_lab <- r_out * sin(theta)

rot <- theta * 180 / pi
rot_flip <- ifelse(rot > 90 | rot < -90, rot + 180, rot)

pos_vec <- ifelse(cos(theta) >= 0, 4, 2) # right vs left

tip_cex <- 0.20 # increase/decrease

for (k in seq_along(tip_idx)) {
  text(
    x = x_lab[k],
    y = y_lab[k],
    labels = phy$tip.label[tip_idx[k]],
    cex = tip_cex,
    srt = rot_flip[k],
    pos = pos_vec[k],
    offset = 0
  )
}

# sanity check: how many labels are we drawing?
n_labeled <- length(tip_idx)
expected <- ceiling(n_tips / label_every)

message("Total tips: ", n_tips)
message(
  "Labeling every ",
  label_every,
  " tips -> labeled: ",
  n_labeled,
  " (expected: ",
  expected,
  ")"
)
stopifnot(n_labeled == expected)

# optional: check proportion ~ 1/label_every
prop <- n_labeled / n_tips
message(
  "Proportion labeled: ",
  round(prop, 3),
  " (~",
  round(1 / label_every, 3),
  ")"
)

# ------------------------------------------------------------------------------
# Family labels and MRCAs
# ------------------------------------------------------------------------------
family_list <- unique(traits$family)

family_nodes <- sapply(family_list, function(fam) {
  species <- traits %>% filter(family == fam) %>% pull(tips)
  tryCatch(findMRCA(phy, species, type = "node"), error = function(e) NA)
})

family_nodes <- na.omit(family_nodes)

r_ring <- max(sqrt(tips_xx^2 + tips_yy^2)) + ri * (space + len)

r_label <- r_ring + 0.02 * r_ring

for (i in seq_along(family_nodes)) {
  fam <- names(family_nodes)[i]
  fam_tips <- which(traits$family == fam)
  fam_tip_labels <- traits$tips[fam_tips]
  fam_tip_indices <- match(fam_tip_labels, phy$tip.label)
  x_vals <- obj$xx[fam_tip_indices]
  y_vals <- obj$yy[fam_tip_indices]
  thetas <- atan2(y_vals, x_vals)
  mean_theta <- atan2(mean(sin(thetas)), mean(cos(thetas)))
  x_out <- r_label * cos(mean_theta)
  y_out <- r_label * sin(mean_theta)
  hjust <- ifelse(cos(mean_theta) >= 0, 0, 1)
  text(x_out, y_out, labels = fam, cex = 0.6, srt = 0, adj = c(hjust, 0.5))
}

mrca_color <- "#9467bd"
for (node in family_nodes) {
  points(obj$xx[node], obj$yy[node], pch = 16, cex = 1, col = mrca_color)
}

# ------------------------------------------------------------------------------
# Legend
# ------------------------------------------------------------------------------
legend(
  "topleft",
  legend = names(combo_colors),
  col = combo_colors,
  pch = 16,
  pt.cex = 0.9,
  title = expression(bold("Trait Combination")),
  bty = "n",
  cex = 0.6,
  x.intersp = 1.3
)

dev.off()

#-------------------------------------------------------------------------------
# Done!
#-------------------------------------------------------------------------------
