#===============================================================================
# Climatic niche breadth
#===============================================================================
# Evaluates how combinations of sociality and nesting strategy influence
# bees’ climatic niche breadth, both univariately (temperature and precipitation)
# and multivariately (using PCA)
#===============================================================================
# Setup
#-------------------------------------------------------------------------------
rm(list = ls())

data_dir <- "/home/lenarh/data/bee_nesting_sociality"
curated_data_dir <- file.path(data_dir, "curated_data")
results_dir <- file.path(data_dir, "results", "PCA")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(dplyr)
library(viridis)
library(colorspace)
library(ggrepel)
library(gridExtra)
library(factoextra)
library(caret)
library(corrplot)
library(vegan)
library(alphashape3d)
library(grid)

# Load trait and climate summary data
traits <- read.csv(file.path(curated_data_dir, "bee_traits_clean.csv"))

all_climatic_vars <- list.files(
  file.path(curated_data_dir, "climate_summaries"),
  pattern = "summstats\\.csv$",
  full.names = TRUE
)

#===============================================================================
# 1) Univariate niche breadths (temp and precip)
#===============================================================================
# Found by taking the difference between the min/max temperatures of the coldest and warmest months,
# and the min/max precipitations of the wettest and driest months
#===============================================================================
# Filter for BIOCLIM variables corresponding to:
# - BIO5: Max temperature of warmest month
# - BIO6: Min temperature of coldest month
# - BIO13: Precipitation of wettest month
# - BIO14: Precipitation of driest month

climatic_list <- lapply(all_climatic_vars, read.csv)

# Merge climate means into one dataframe based on species
merged_climatic_vars <- climatic_list[[1]]

for (i in 2:length(climatic_list)) {
  one_climatic_var <- climatic_list[[i]]
  merged_climatic_vars <- merge(
    merged_climatic_vars,
    one_climatic_var,
    by = "species"
  )
}

# Keep only mean values (averaged across occurrence points)
merged_climatic_vars <- merged_climatic_vars[, c(
  1,
  grep("mean", colnames(merged_climatic_vars))
)]

# Rename CHELSA climate columns to simpler BIO names
colnames(merged_climatic_vars) <- gsub(
  "mean_CHELSA_bio([0-9]{2})_1981.2010_V.2.1",
  "mean_bio_\\1",
  colnames(merged_climatic_vars)
)

# Merge climate data with trait data using species names
merged_traits <- merge(
  traits,
  merged_climatic_vars,
  by.x = "tips",
  by.y = "species"
)

# Combine trait states, e.g., "social_aboveground"
merged_traits$combined_trait <- paste(
  merged_traits$sociality_binary,
  merged_traits$nest_binary,
  sep = "_"
)

# Removing NAs and NANs
merged_traits <- subset(merged_traits, !is.nan(merged_traits$mean_bio_05))
merged_traits <- subset(merged_traits, !is.na(merged_traits$mean_bio_05))

# Check data
head(merged_traits)

# Calculate univariate niche breadths for each species
# temp_breadth = BIO5 - BIO6 (temperature range)
# prec_breadth = BIO13 - BIO14 (precipitation range)
merged_traits$prec_breadth <- NA
merged_traits$temp_breadth <- NA

for (i in 1:nrow(merged_traits)) {
  merged_traits$prec_breadth[i] <- mean(merged_traits$mean_bio_13[i]) -
    mean(merged_traits$mean_bio_14[i])
  merged_traits$temp_breadth[i] <- mean(merged_traits$mean_bio_05[i]) -
    mean(merged_traits$mean_bio_06[i])
}

head(merged_traits$prec_breadth)
head(merged_traits$temp_breadth)
table(merged_traits$combined_trait)

head(merged_traits)

#-------------------------------------------------------------------------------
# PERMANOVA on univariate climatic niche space
#-------------------------------------------------------------------------------
# Temperature
#-------------------------------------------------------------------------------
adonis_data_temp <- data.frame(x = merged_traits$temp_breadth)
groups <- merged_traits$combined_trait

# Run PERMANOVA
permanova_result_temp <- adonis2(
  adonis_data_temp ~ groups,
  method = "euclidean",
  permutations = 999
)

# Print results
print(permanova_result_temp)

# Convert to df
permanova_df_temp <- as.data.frame(permanova_result_temp)
head(permanova_df_temp)

# Save PERMANOVA
write.csv(
  permanova_df_temp,
  file = file.path(results_dir, "permanova_temp.csv"),
  row.names = TRUE
)

# R^2 = 0.177
# F = 266.27
# P = 0.001 with 999 permutations
#     With 999 permutations, the smallest possible p-value you can report is 0.001,
#     This means none of the 999 random groupings had as much between-group difference as the real data.

#-------------------------------------------------------------------------------
# Precipitation
#-------------------------------------------------------------------------------
adonis_data_precip <- data.frame(x = merged_traits$prec_breadth)
groups <- merged_traits$combined_trait

# Run PERMANOVA
permanova_result_precip <- adonis2(
  adonis_data_precip ~ groups,
  method = "euclidean",
  permutations = 999
)

# Print results
print(permanova_result_precip)

# Convert to df
permanova_df_precip <- as.data.frame(permanova_result_precip)
head(permanova_df_precip)

# Save PERMANOVA
write.csv(
  permanova_df_precip,
  file = file.path(results_dir, "permanova_precip.csv"),
  row.names = TRUE
)

# R^2 = 0.111
# F = 154.29
# P = 0.001

#-------------------------------------------------------------------------------
# Pairwise PERMANOVA on univariate climatic niche space
#-------------------------------------------------------------------------------
# Temperature
#-------------------------------------------------------------------------------
groups <- unique(merged_traits$combined_trait)

# Store results
pairwise_results_temp <- data.frame(
  Group1 = character(),
  Group2 = character(),
  R2 = numeric(),
  F_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through all unique pairs
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]

    # Subset to the two groups and keep finite temp_breadth
    subset_df <- merged_traits %>%
      dplyr::filter(combined_trait %in% c(g1, g2), is.finite(temp_breadth))

    # UNIVARIATE predictor as a 1-col data frame
    data <- data.frame(x = subset_df$temp_breadth)
    group <- droplevels(factor(subset_df$combined_trait))

    # Run PERMANOVA
    result <- adonis2(data ~ group, method = "euclidean", permutations = 999)

    # Store in results table
    pairwise_results_temp <- rbind(
      pairwise_results_temp,
      data.frame(
        Group1 = g1,
        Group2 = g2,
        R2 = result$R2[1],
        F_value = result$F[1],
        p_value = result$`Pr(>F)`[1]
      )
    )
  }
}

# Print all pairwise results
print(pairwise_results_temp)

# Save pairwise PERMANOVA
write.csv(
  pairwise_results_temp,
  file = file.path(results_dir, "permanova_pairwise_temp.csv"),
  row.names = FALSE
)

#                 Group1               Group2          R2    F_value p_value
# 1      solitary_ground solitary_aboveground 0.063078835 189.993009   0.001
# 2      solitary_ground   social_aboveground 0.278931076 833.231772   0.001
# 3      solitary_ground        social_ground 0.001347291   3.388962   0.069
# 4 solitary_aboveground   social_aboveground 0.147329839 209.417161   0.001
# 5 solitary_aboveground        social_ground 0.068671355 115.763675   0.001
# 6   social_aboveground        social_ground 0.398280156 597.036485   0.001

#-------------------------------------------------------------------------------
# Precipitation
#-------------------------------------------------------------------------------
groups <- unique(merged_traits$combined_trait)

# Store results
pairwise_results_precip <- data.frame(
  Group1 = character(),
  Group2 = character(),
  R2 = numeric(),
  F_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through all unique pairs
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]

    # Subset to the two groups and keep finite prec_breadth
    subset_df <- merged_traits %>%
      dplyr::filter(combined_trait %in% c(g1, g2), is.finite(prec_breadth))

    # Univariate predictor as a 1-col data frame
    data <- data.frame(x = subset_df$prec_breadth)
    group <- droplevels(factor(subset_df$combined_trait))

    # Run PERMANOVA
    result <- adonis2(data ~ group, method = "euclidean", permutations = 999)

    # Store in results table
    pairwise_results_precip <- rbind(
      pairwise_results_precip,
      data.frame(
        Group1 = g1,
        Group2 = g2,
        R2 = result$R2[1],
        F_value = result$F[1],
        p_value = result$`Pr(>F)`[1]
      )
    )
  }
}

# Print all pairwise results
print(pairwise_results_precip)

# Save pairwise PERMANOVA
write.csv(
  pairwise_results_precip,
  file = file.path(results_dir, "permanova_pairwise_precip.csv"),
  row.names = FALSE
)

#                 Group1               Group2          R2   F_value p_value
# 1      solitary_ground solitary_aboveground 0.049775833 147.82554   0.001
# 2      solitary_ground   social_aboveground 0.193681928 517.40236   0.001
# 3      solitary_ground        social_ground 0.009262984  23.48617   0.001
# 4 solitary_aboveground   social_aboveground 0.083621648 110.59781   0.001
# 5 solitary_aboveground        social_ground 0.012685066  20.17143   0.001
# 6   social_aboveground        social_ground 0.174559997 190.75053   0.001

#===============================================================================
# 2) Multivariate climatic niche space and niche volume estimation
#===============================================================================
# Position of each species in the climatic multivariate space,
# plus their volume and correlation between volume and life history.
#===============================================================================
# Reload all climate summary stats for full PCA analysis
all_climatic_vars <- list.files(
  file.path(curated_data_dir, "climate_summaries"),
  pattern = "summstats\\.csv$",
  full.names = TRUE
)

climatic_list <- lapply(all_climatic_vars, read.csv)

# Merge all climatic variables across species
merged_climatic_vars <- climatic_list[[1]]

for (i in 2:length(climatic_list)) {
  one_climatic_var <- climatic_list[[i]]
  merged_climatic_vars <- merge(
    merged_climatic_vars,
    one_climatic_var,
    by = "species"
  )
}

# Merge with trait data again
merged_traits_pca <- merge(
  traits,
  merged_climatic_vars,
  by.x = "tips",
  by.y = "species"
)

# Rename CHELSA climate columns to simpler names
names(merged_traits_pca) <- gsub(
  "mean_CHELSA_bio([0-9]+)_1981\\.2010_V\\.2\\.1.*",
  "mean_bio_\\1",
  names(merged_traits_pca)
)

names(merged_traits_pca) <- gsub(
  "sd_CHELSA_bio([0-9]+)_1981\\.2010_V\\.2\\.1.*",
  "sd_bio_\\1",
  names(merged_traits_pca)
)

names(merged_traits_pca) <- gsub(
  "se_CHELSA_bio([0-9]+)_1981\\.2010_V\\.2\\.1.*",
  "se_bio_\\1",
  names(merged_traits_pca)
)

names(merged_traits_pca) <- gsub(
  "within_sp_var_CHELSA_bio([0-9]+)_1981\\.2010_V\\.2\\.1.*",
  "within_sp_var_bio_\\1",
  names(merged_traits_pca)
)

names(merged_traits_pca) <- gsub(
  "n_CHELSA_bio([0-9]+)_1981\\.2010_V\\.2\\.1.*",
  "n_bio_\\1",
  names(merged_traits_pca)
)

head(merged_traits_pca)

# Remove rows with missing values
merged_traits_pca <- subset(merged_traits_pca, !is.na(mean_bio_01))
merged_traits_pca <- subset(merged_traits_pca, !is.na(mean_bio_12))

# Extract only columns with climate mean values
cols <- grep("^mean_bio_", colnames(merged_traits_pca))
head(cols)

# Check for NA or non-finite values in climate columns
any(is.na(merged_traits_pca[, cols]))
any(!is.finite(as.matrix(merged_traits_pca[, cols])))

# Keep only complete cases and re-define combined trait grouping
merged_traits_pca <- merged_traits_pca[
  complete.cases(merged_traits_pca[, cols]),
]
merged_traits_pca$combined_trait <- paste(
  merged_traits_pca$sociality_binary,
  merged_traits_pca$nest_binary,
  sep = "_"
)

# Extract cleaned climate matrix + traits
clim_vars <- merged_traits_pca[, cols]
clim_vars$tips <- merged_traits_pca$tips
clim_vars$combined_trait <- merged_traits_pca$combined_trait
cols <- grep("^mean_bio_", colnames(clim_vars))

# Remove highly collinear variables (r > 0.9)
cor_matrix <- cor(clim_vars[, cols])

corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8)

highly_correlated <- findCorrelation(cor_matrix, cutoff = 0.9)

cleaned_clim_vars <- clim_vars[, -highly_correlated]

# Principal component analysis (PCA)
# Standardizes variables (mean = 0, SD = 1)
pca_cols <- grep("^mean_bio_", colnames(cleaned_clim_vars))
pca_result <- prcomp(cleaned_clim_vars[, pca_cols], scale. = TRUE)

# Create a data frame with PCA scores for each species
pca_df <- as.data.frame(pca_result$x)

# Add the grouping variable
pca_df$combined_trait <- cleaned_clim_vars$combined_trait

#-------------------------------------------------------------------------------
# Calculate the volume of the climatic space occupied by each trait combination
#-------------------------------------------------------------------------------
groups <- unique(pca_df$combined_trait)

volumes <- c()

for (i in 1:length(groups)) {
  ashape3d.occ <- ashape3d(
    as.matrix(pca_df[pca_df$combined_trait == groups[i], 1:3]),
    alpha = 2
  )
  volumes <- c(volumes, volume_ashape3d(ashape3d.occ)) # niche size
}

names(volumes) <- groups

# Multivariate climatic niche volumes per trait group
print(volumes)

# Convert to df
volumes_df <- data.frame(
  combined_trait = names(volumes),
  niche_volume = as.numeric(volumes)
)

head(volumes_df)

# Save niche volumes
write.csv(
  volumes_df,
  file = file.path(results_dir, "multivariate_niche_volumes.csv"),
  row.names = FALSE
)

#         combined_trait niche_volume     Trait_Combination
# 1      solitary_ground    201.80420       Solitary/Ground
# 2 solitary_aboveground    247.00729 Solitary/Above-ground
# 3   social_aboveground     56.46831   Social/Above-ground
# 4        social_ground    197.58209         Social/Ground

#-------------------------------------------------------------------------------
# PERMANOVA on multivariate climatic niche space
#-------------------------------------------------------------------------------
# Tests whether species grouped by trait combination occupy significantly
# different regions of climate space (based on PCA scores).
#-------------------------------------------------------------------------------
groups <- pca_df$combined_trait

# Use the first 3 PCs
adonis_data <- pca_df[, c("PC1", "PC2", "PC3")]

# Run PERMANOVA
permanova_results_multivar <- adonis2(
  adonis_data ~ groups,
  method = "euclidean",
  permutations = 999
)

# Print results
print(permanova_results_multivar)

# Convert to df and save
permanova_df_multivar <- as.data.frame(permanova_results_multivar)
head(permanova_df_multivar)

# Save PERMANOVA
write.csv(
  permanova_df_multivar,
  file = file.path(results_dir, "permanova_multivariate.csv"),
  row.names = TRUE
)

# F = 180.88
# R^2 = 0.127
# p = 0.001

#-------------------------------------------------------------------------------
# Pairwise PERMANOVA on multivariate climatic niche space
#-------------------------------------------------------------------------------
groups <- unique(pca_df$combined_trait)

# Store results
pairwise_results_multivariate <- data.frame(
  Group1 = character(),
  Group2 = character(),
  R2 = numeric(),
  F_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through all unique pairs
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]

    # Subset data for just these two groups
    subset_df <- pca_df %>%
      filter(combined_trait %in% c(g1, g2))

    # Get PCA scores and group variable
    data <- subset_df[, c("PC1", "PC2", "PC3")]
    group <- subset_df$combined_trait

    # Run PERMANOVA
    result <- adonis2(data ~ group, method = "euclidean", permutations = 999)

    # Store in results table
    pairwise_results_multivariate <- rbind(
      pairwise_results_multivariate,
      data.frame(
        Group1 = g1,
        Group2 = g2,
        R2 = result$R2[1],
        F_value = result$F[1],
        p_value = result$`Pr(>F)`[1]
      )
    )
  }
}

# Print all pairwise results
print(pairwise_results_multivariate)

# Save pairwise PERMANOVA
write.csv(
  pairwise_results_multivariate,
  file = file.path(results_dir, "permanova_pairwise_multivariate.csv"),
  row.names = FALSE
)

#                 Group1               Group2         R2   F_value p_value
# 1      solitary_ground solitary_aboveground 0.03255442  94.95993   0.001
# 2      solitary_ground   social_aboveground 0.16766950 433.91430   0.001
# 3      solitary_ground        social_ground 0.03984395 104.24140   0.001
# 4 solitary_aboveground   social_aboveground 0.10859695 147.65432   0.001
# 5 solitary_aboveground        social_ground 0.05253319  87.05012   0.001
# 6   social_aboveground        social_ground 0.30580592 397.34844   0.001

#===============================================================================
# 3) Plotting
#===============================================================================
# Setup required for ALL plotting
#-------------------------------------------------------------------------------

# Set combined_trait factor level order
merged_traits$combined_trait <- factor(
  merged_traits$combined_trait,
  levels = c(
    "solitary_ground",
    "solitary_aboveground",
    "social_ground",
    "social_aboveground"
  )
)

# Trait labels
trait_labels <- c(
  "solitary_ground" = "Solitary/Ground",
  "solitary_aboveground" = "Solitary/Above-ground",
  "social_ground" = "Social/Ground",
  "social_aboveground" = "Social/Above-ground"
)

# Trait colors
trait_colors <- viridis::viridis(n = 4, option = "cividis")
names(trait_colors) <- levels(merged_traits$combined_trait)

#-------------------------------------------------------------------------------
# Plot niche volumes
#-------------------------------------------------------------------------------
# Apply labels
volumes_df$Trait_Combination <- factor(
  volumes_df$combined_trait,
  levels = names(trait_labels),
  labels = trait_labels
)

# Barplot of niche volumes with same aesthetics as plot1
plot_volumes <- ggplot(
  volumes_df,
  aes(x = Trait_Combination, y = niche_volume, fill = Trait_Combination)
) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(niche_volume, 1)), vjust = -0.5, size = 5) +
  scale_fill_viridis_d(option = "cividis", name = "Trait Combination") +
  labs(x = "", y = "Niche Volume") +
  theme_minimal() +
  theme(
    legend.position = c(1.00, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(
      size = 16,
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(size = 16, color = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

print(plot_volumes)

# Save to file
ggsave(
  file.path(results_dir, "multivariate_niche_volumes_barplot.pdf"),
  plot_volumes,
  width = 8,
  height = 7,
  dpi = 300
)

#-------------------------------------------------------------------------------
# Plot 2D PCA
#-------------------------------------------------------------------------------
pca_plot <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = combined_trait, fill = combined_trait)
) +
  stat_ellipse(
    aes(group = combined_trait),
    type = "norm",
    geom = "polygon",
    alpha = 0.2,
    color = NA
  ) +
  stat_ellipse(
    aes(group = combined_trait),
    type = "norm",
    size = 0.5,
    linetype = "solid"
  ) +
  geom_point(size = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  scale_color_manual(values = trait_colors, labels = trait_labels) +
  scale_fill_manual(values = trait_colors, guide = "none") +
  theme_classic() +
  labs(
    title = "",
    color = "Trait Combination",
    x = "PC1: Cool, dry, seasonal (T) ↔ Wet, warm, aseasonal (T)",
    y = ""
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.title = element_text(size = 0),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20),
    legend.position = "right",
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

pca_plot

# Save PCA plot
ggsave(file.path(results_dir, "PCA.pdf"), pca_plot, width = 14, height = 8)

#-------------------------------------------------------------------------------
# Plot univariate niches as violin plots
#-------------------------------------------------------------------------------
# Convert niche breadths to interpretable units (rather than raw raster units) for plotting - added 7/26/26
plot_traits <- merged_traits %>%
  mutate(
    prec_breadth_plot = prec_breadth / 10,
    temp_breadth_plot = temp_breadth / 10
  )

# Define point colors (darker versions of trait colors)
point_colors <- sapply(trait_colors, function(x) darken(x, amount = 0.4))

# Precipitation breadth violin plot
# precip_breadth <- ggplot(
#   merged_traits,
#   aes(x = prec_breadth, y = combined_trait, fill = combined_trait)
# ) +
precip_breadth <- ggplot(
  # changed 7/26/26
  plot_traits,
  aes(x = prec_breadth_plot, y = combined_trait, fill = combined_trait)
) +
  geom_violin(width = 0.9, adjust = 0.8, alpha = 0.6, show.legend = FALSE) +
  geom_jitter(
    aes(color = combined_trait),
    height = 0.2,
    alpha = 0.7,
    size = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = trait_colors) +
  scale_color_manual(values = point_colors) +
  labs(
    x = "BIO13 - BIO14 (mm)",
    y = "",
    title = "Precipitation Niche Breadth"
  ) +
  scale_y_discrete(labels = trait_labels) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

# Temperature breadth violin plot
# temp_breadth <- ggplot(
#   merged_traits,
#   aes(x = temp_breadth, y = combined_trait, fill = combined_trait)
# ) +
temp_breadth <- ggplot(
  # changed 7/26/26
  plot_traits,
  aes(x = temp_breadth_plot, y = combined_trait, fill = combined_trait)
) +
  geom_violin(width = 0.9, adjust = 0.8, alpha = 0.6, show.legend = FALSE) +
  geom_jitter(
    aes(color = combined_trait),
    height = 0.2,
    alpha = 0.7,
    size = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = trait_colors) +
  scale_color_manual(values = point_colors) +
  labs(x = "BIO5 - BIO6 (°C)", y = "", title = "Temperature Niche Breadth") +
  scale_y_discrete(labels = trait_labels) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

precip_breadth
temp_breadth

# Combine into a two-panel figure (vertical layout)
combined_univariate_niches <- grid.arrange(
  precip_breadth,
  temp_breadth,
  ncol = 1
)

# Save niche breadth violin plots
ggsave(
  file.path(results_dir, "niche_breadth_violins.pdf"),
  combined_univariate_niches,
  width = 10,
  height = 12
)

#-------------------------------------------------------------------------------
# Plot PCA loading arrows
#-------------------------------------------------------------------------------
# Extract PCA variable loadings and contributions
pca_vars <- get_pca_var(pca_result)
contributions <- pca_vars$contrib

loadings_df <- as.data.frame(pca_vars$coord)
loadings_df$variable <- rownames(loadings_df)

# Replace variable codes with descriptive labels
bio_labels <- c(
  mean_bio_01 = "Annual Temp",
  mean_bio_02 = "Diurnal Temp Range",
  mean_bio_03 = "Isothermality",
  mean_bio_04 = "Temp Seasonality",
  mean_bio_05 = "Max Temp (Warmest Month)",
  mean_bio_06 = "Min Temp (Coldest Month)",
  mean_bio_07 = "Annual Temp Range",
  mean_bio_08 = "Mean Temp (Wettest Qtr)",
  mean_bio_09 = "Mean Temp (Driest Qtr)",
  mean_bio_10 = "Mean Temp (Warmest Qtr)",
  mean_bio_11 = "Mean Temp (Coldest Qtr)",
  mean_bio_12 = "Annual Precip",
  mean_bio_13 = "Precip (Wettest Month)",
  mean_bio_14 = "Precip (Driest Month)",
  mean_bio_15 = "Precip Seasonality",
  mean_bio_16 = "Precip (Wettest Qtr)",
  mean_bio_17 = "Precip (Driest Qtr)",
  mean_bio_18 = "Precip (Warmest Qtr)",
  mean_bio_19 = "Precip (Coldest Qtr)"
)

loadings_df$label <- bio_labels[loadings_df$variable]

# Filter to top contributors (≥ 9% to PC1 or PC2)
top_vars <- contributions[, 1:2] >= 9
top_var_names <- rownames(contributions)[apply(top_vars, 1, any)]
loadings_df <- loadings_df[loadings_df$variable %in% top_var_names, ]

# Add contribution value for coloring
loadings_df$contrib <- apply(contributions[loadings_df$variable, 1:2], 1, max)

# Plot
max_radius <- sqrt(max(loadings_df$Dim.1^2 + loadings_df$Dim.2^2, na.rm = TRUE)) # Calculate radius based on longest arrow

arrow_plot <- ggplot(loadings_df, aes(x = 0, y = 0)) +
  # Arrows
  geom_segment(
    aes(xend = Dim.1, yend = Dim.2, color = contrib),
    arrow = arrow(length = unit(0.25, "cm")),
    size = 1
  ) +
  # Labels
  geom_text_repel(
    aes(x = Dim.1, y = Dim.2, label = label),
    size = 4,
    color = "black"
  ) +
  # Axes
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  # Reference circle (scaled to longest arrow)
  geom_path(
    data = data.frame(
      x = max_radius * cos(seq(0, 2 * pi, length.out = 300)),
      y = max_radius * sin(seq(0, 2 * pi, length.out = 300))
    ),
    aes(x = x, y = y),
    linetype = "dotted",
    color = "gray40",
    linewidth = 0.6
  ) +
  # Color scale with whole-number breaks
  scale_color_viridis_c(
    name = "Max Contribution\n(PC1 or PC2)",
    breaks = pretty(range(loadings_df$contrib), n = 5) |> floor() |> unique()
  ) +
  coord_equal() +
  theme_classic() +
  labs(x = "PC1", y = "PC2") +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "top",
    legend.justification = "right",
    legend.box = "horizontal",
    legend.background = element_blank(),
    plot.margin = margin(20, 20, 10, 10)
  )

# Save
ggsave(
  file.path(results_dir, "pca_loadings.pdf"),
  arrow_plot,
  width = 9,
  height = 8.5
)

# Arrow direction = the direction of the variable's influence on PC1 and PC2 (into which quadrant)

# Arrow length = the magnitude of the loading (i.e., how strongly the variable is associated with those axes)

# Arrow color = the contribution, i.e., how much that variable contributes to explaining variance along PC1 or PC2

#-------------------------------------------------------------------------------
# View PCA loadings and contributions
#-------------------------------------------------------------------------------
# Cumulative variance explained by first 2 PCs
eig <- (pca_result$sdev)^2
variance <- eig * 100 / sum(eig)
sum(variance[1:2]) # PC1 + PC2 explain 67.17% of variation
variance # PC1 explains 40.05%, PC2 explains 27.12%, PC3 explains 11.85%

# Get PCA variable information
pca_vars <- get_pca_var(pca_result)

# Loadings: directions of variables in PC space
print(round(pca_vars$coord[, 1:2], 2)) # Loadings for PC1–PC2

# Contributions: how much each variable contributes (%) to each PC
print(round(pca_vars$contrib[, 1:2], 2)) # Contributions for PC1–PC2

# Filter to top contributors (≥ 9% to PC1 or PC2)
top_vars <- contributions[, 1:2] >= 9
top_var_names <- rownames(contributions)[apply(top_vars, 1, any)]
top_var_names
#  [1] "mean_bio_01" "mean_bio_02" "mean_bio_03" "mean_bio_04" "mean_bio_08" "mean_bio_09" "mean_bio_13"
#  [8] "mean_bio_14" "mean_bio_15" "mean_bio_18" "mean_bio_19"

# Extract loadings and contributions for top variables
top_loadings <- as.data.frame(pca_vars$coord[top_var_names, 1:2])
top_contribs <- as.data.frame(pca_vars$contrib[top_var_names, 1:2])

# Add variable names as a column
top_loadings$Variable <- rownames(top_loadings)
top_contribs$Variable <- rownames(top_contribs)

# Merge both into one table
top_pca_summary <- merge(top_loadings, top_contribs, by = "Variable")

# Rename columns
colnames(top_pca_summary) <- c(
  "Variable",
  "PC1_Loading",
  "PC2_Loading",
  "PC1_Contribution",
  "PC2_Contribution"
)

# Round values to 2 decimal places
top_pca_summary[, 2:5] <- round(top_pca_summary[, 2:5], 2)

head(top_pca_summary)

#      Variable PC1_Loading PC2_Loading PC1_Contribution PC2_Contribution
# 1 mean_bio_01        0.90       -0.35            16.88             3.73
# 2 mean_bio_02       -0.28       -0.71             1.63            15.70
# 3 mean_bio_03        0.81       -0.18            13.79             1.03
# 4 mean_bio_04       -0.86        0.02            15.32             0.01
# 5 mean_bio_08        0.72       -0.10            10.78             0.29
# 6 mean_bio_09        0.71       -0.49            10.52             7.27

# Save PCA loadings/contributions
write.csv(
  top_pca_summary,
  file.path(results_dir, "pca_top_contributors_PC1_PC2.csv"),
  row.names = FALSE
)
