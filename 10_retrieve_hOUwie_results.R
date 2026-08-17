# ==============================================================================
# hOUwie Results Processing
# ==============================================================================
# Loop over the four climatic variables of interest (bio1, bio12, bio15, bio4)
#
# For each variable,
#   Load the relevant .Rsave files
#   Extract model comparison tables
#   Compute model-averaged parameters
#   Save CSVs of the above
#   Plot expected means and variances
# ==============================================================================
# Setup
# ------------------------------------------------------------------------------
rm(list = ls())

data_dir <- "/home/lenarh/data/bee_nesting_sociality"
repo_dir <- "/home/lenarh/repos/bee_nesting_sociality"
results_dir <- file.path(data_dir, "results", "hOUwie")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(OUwie)
library(corHMM)
library(parallel)

source(file.path(repo_dir, "00_utility_functions.R"))

# ------------------------------------------------------------------------------
# Load hOUwie results
# ------------------------------------------------------------------------------
# Define variables of interest
bio_vars <- c("bio_1", "bio_4", "bio_12", "bio_15")

# Loop through each climate variable
for (bio in bio_vars) {
  message("Processing: ", bio)
  # Define expected output file paths
  csv_table_path <- file.path(
    results_dir,
    paste0("hOUwie_model_table_", bio, ".csv")
  )
  csv_pars_path <- file.path(
    results_dir,
    paste0("hOUwie_average_pars_", bio, ".csv")
  )
  plot_mean_path <- file.path(
    results_dir,
    paste0("hOUwie_boxplot_", bio, "_mean.pdf")
  )
  plot_var_path <- file.path(
    results_dir,
    paste0("hOUwie_boxplot_", bio, "_variance.pdf")
  )

  # Skip if outputs already exist (comment this out if you want to re-run all)
  # if (all(file.exists(csv_table_path, csv_pars_path, plot_mean_path, plot_var_path))) {
  #   message("Skipping ", bio, " (results already exist)")
  #   next
  # }

  # Subset relevant .Rsave files
  all_model_results <- list.files(results_dir, full.names = TRUE)
  bio_model_files <- all_model_results[grepl(
    paste0("_", bio, "\\.Rsave$"),
    all_model_results
  )]
  model_names <- gsub(".Rsave", "", basename(bio_model_files))

  # Load each result into a list
  all_results <- list()
  for (i in seq_along(bio_model_files)) {
    tryCatch(
      {
        load(bio_model_files[i]) # Loads 'res'
        if (exists("res") && !is.null(res)) {
          all_results[[i]] <- res
          names(all_results)[i] <- model_names[i]
        }
        rm(res)
      },
      error = function(e) {
        warning(
          "Failed to load or process: ",
          bio_model_files[i],
          "\n",
          e$message
        )
      }
    )
  }

  # Skip if no models loaded
  if (length(all_results) == 0) {
    warning("No valid models found for ", bio)
    next
  }

  # Generate and save model comparison table
  model_table <- getModelTable(all_results, type = "AICc")
  print(model_table)
  write.csv(model_table, file = csv_table_path, row.names = TRUE)

  # Get model-averaged parameters
  average_pars <- getModelAvgParams(all_results, type = "AICc")

  # Convert expected means to original units
  if (bio == "bio_1") {
    average_pars$expected_mean <- exp(average_pars$expected_mean) - 273.15
    ylab <- "Expected mean (°C)"
  } else if (bio == "bio_4") {
    average_pars$expected_mean <- exp(average_pars$expected_mean)
    ylab <- "Expected mean (temperature seasonality)"
  } else if (bio == "bio_12") {
    average_pars$expected_mean <- exp(average_pars$expected_mean)
    ylab <- "Expected mean (precipitation)"
  } else if (bio == "bio_15") {
    average_pars$expected_mean <- exp(average_pars$expected_mean)
    ylab <- "Expected mean (precipitation seasonality)"
  } else {
    ylab <- "Expected mean (log-transformed units)"
  }

  # Round for readability
  average_pars$expected_mean <- round(average_pars$expected_mean, 2)
  print(range(average_pars$expected_mean))

  # Save average parameters (after converting back to original units)
  write.csv(average_pars, file = csv_pars_path, row.names = TRUE)

  # Save boxplots
  pdf(file = plot_mean_path)
  boxplot(
    average_pars$expected_mean ~ average_pars$tip_state,
    xlab = "Tip state",
    ylab = ylab,
    main = paste("Expected Means -", bio)
  )
  dev.off()

  pdf(file = plot_var_path)
  boxplot(
    average_pars$expected_var ~ average_pars$tip_state,
    xlab = "Tip state",
    ylab = "Expected variance",
    main = paste("Expected Variances -", bio)
  )
  dev.off()
}

#-------------------------------------------------------------------------------
# Explore data
#-------------------------------------------------------------------------------
library(dplyr)

bio1 <- read.csv(file.path(results_dir, "hOUwie_average_pars_bio_1.csv"))
bio4 <- read.csv(file.path(results_dir, "hOUwie_average_pars_bio_4.csv"))
bio12 <- read.csv(file.path(results_dir, "hOUwie_average_pars_bio_12.csv"))
bio15 <- read.csv(file.path(results_dir, "hOUwie_average_pars_bio_15.csv"))

head(bio1)

bio1_summary <- bio1 %>%
  group_by(tip_state) %>%
  summarize(
    mean_temp = round(mean(expected_mean)),
    .groups = "drop"
  )

print(bio1_summary)

# tip_state            mean_temp (C)
# social_aboveground        22.1
# social_ground             13.8
# solitary_aboveground      16.3
# solitary_ground           16.9

# 5/15/26:
# tip_state            mean_temp (C)
# social_aboveground          16
# social_ground               16
# solitary_aboveground        16
# solitary_ground             16

bio4_summary <- bio4 %>%
  group_by(tip_state) %>%
  summarize(
    mean_seasonality = round(mean(expected_mean)),
    .groups = "drop"
  )

print(bio4_summary)

# tip_state            mean_seasonality (SD x 100)
# social_aboveground               186
# social_ground                    618
# solitary_aboveground             351
# solitary_ground                  430

# 5/15/26:
# tip_state            mean_seasonality (SD x 100)
# social_aboveground                402
# social_ground                     268
# solitary_aboveground              400
# solitary_ground                   265

bio12_summary <- bio12 %>%
  group_by(tip_state) %>%
  summarize(
    mean_precip = round(mean(expected_mean)),
    .groups = "drop"
  )

print(bio12_summary)

# tip_state            mean_precip (mm)
# social_aboveground          665
# social_ground               727
# solitary_aboveground        729
# solitary_ground             500

# 5/15/26:
# tip_state            mean_precip (mm)
# social_aboveground          1230
# social_ground               1225
# solitary_aboveground         860
# solitary_ground              673

bio15_summary <- bio15 %>%
  group_by(tip_state) %>%
  summarize(
    mean_precip_seasonality = round(mean(expected_mean), 2),
    .groups = "drop"
  )

print(bio15_summary)


#-------------------------------------------------------------------------------
# Plotting
#-------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)
library(colorspace)
library(patchwork)

# Load results
bio1 <- read.csv(file.path(results_dir, "hOUwie_average_pars_bio_1.csv"))
bio4 <- read.csv(file.path(results_dir, "hOUwie_average_pars_bio_4.csv"))
bio12 <- read.csv(file.path(results_dir, "hOUwie_average_pars_bio_12.csv"))
bio15 <- read.csv(file.path(results_dir, "hOUwie_average_pars_bio_15.csv"))

# Add labels
bio1$variable <- "BIO1: Mean Annual Temperature (°C)"
bio4$variable <- "BIO4: Temperature Seasonality (SD × 100)"
bio12$variable <- "BIO12: Annual Precipitation (mm)"
bio15$variable <- "BIO15: Precipitation Seasonality (CV)"

# Combine datasets
all_data <- bind_rows(bio1, bio4, bio12, bio15)

# Define factor levels and color mapping
trait_order <- c(
  "solitary_ground",
  "solitary_aboveground",
  "social_ground",
  "social_aboveground"
)

trait_labels <- c(
  "solitary_ground" = "Solitary/Ground",
  "solitary_aboveground" = "Solitary/Above-ground",
  "social_ground" = "Social/Ground",
  "social_aboveground" = "Social/Above-ground"
)

all_data$tip_state <- factor(all_data$tip_state, levels = trait_order)

box_colors <- viridis::viridis(length(trait_order), option = "cividis")
names(box_colors) <- trait_order

# Generate combined boxplots for expected means
variables <- unique(all_data$variable)
mean_plots <- list()

for (v in variables) {
  subset_data <- filter(all_data, variable == v)

  p <- ggplot(
    subset_data,
    aes(x = tip_state, y = expected_mean, fill = tip_state)
  ) +
    geom_boxplot(
      width = 0.6,
      alpha = 0.6,
      outlier.size = 0.8,
      outlier.alpha = 0.4
    ) +
    scale_fill_manual(values = box_colors) +
    scale_x_discrete(labels = trait_labels) +
    labs(x = "", y = "Expected Mean", title = v) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.text.x = element_text(
        size = 14,
        color = "black",
        angle = 15,
        hjust = 1
      ),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  mean_plots[[v]] <- p
}

# Combine plots
combined_mean_plot <- wrap_plots(mean_plots, ncol = 2) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(size = 18, face = "bold")
    )
  )

# Save
ggsave(
  filename = file.path(results_dir, "hOUwie_combined_means.pdf"),
  plot = combined_mean_plot,
  width = 12,
  height = 12
)

#-------------------------------------------------------------------------------
# Save mean expected climatic optima
#-------------------------------------------------------------------------------
library(dplyr)
library(tidyr)

# (Re)load the model-averaged parameter CSVs
bio1 <- read.csv(file.path(results_dir, "hOUwie_average_pars_bio_1.csv"))
bio4 <- read.csv(file.path(results_dir, "hOUwie_average_pars_bio_4.csv"))
bio12 <- read.csv(file.path(results_dir, "hOUwie_average_pars_bio_12.csv"))
bio15 <- read.csv(file.path(results_dir, "hOUwie_average_pars_bio_15.csv"))

# Summarize means by state for each variable
bio1_summary <- bio1 %>%
  group_by(tip_state) %>%
  summarize(
    mean_bio1 = round(mean(expected_mean, na.rm = TRUE), 1),
    .groups = "drop"
  )

bio4_summary <- bio4 %>%
  group_by(tip_state) %>%
  summarize(
    mean_bio4 = round(mean(expected_mean, na.rm = TRUE), 0),
    .groups = "drop"
  )

bio12_summary <- bio12 %>%
  group_by(tip_state) %>%
  summarize(
    mean_bio12 = round(mean(expected_mean, na.rm = TRUE), 0),
    .groups = "drop"
  )

bio15_summary <- bio15 %>%
  group_by(tip_state) %>%
  summarize(
    mean_bio15 = round(mean(expected_mean, na.rm = TRUE), 2),
    .groups = "drop"
  )

# Join into one table
clim_optima_means <- bio1_summary %>%
  full_join(bio4_summary, by = "tip_state") %>%
  full_join(bio12_summary, by = "tip_state") %>%
  full_join(bio15_summary, by = "tip_state") %>%
  mutate(
    tip_state = factor(
      tip_state,
      levels = c(
        "solitary_ground",
        "solitary_aboveground",
        "social_ground",
        "social_aboveground"
      )
    )
  ) %>%
  arrange(tip_state)

# Save to results_dir
out_path <- file.path(results_dir, "hOUwie_mean_climatic_optima_by_state.csv")
write.csv(clim_optima_means, out_path, row.names = FALSE)
print(clim_optima_means)

# ORIGINAL SCRIPT BELOW
#
# #-------------------------------------------------------------------------------
# # Setup: Load libraries, set working directory, and source helper functions
# #-------------------------------------------------------------------------------
# #rm(list=ls())
# #install.packages("devtools")
# #devtools::install_github("thej022214/OUwie")
# library(OUwie)
# library(corHMM)
# library(parallel)
# setwd("/Users/lenarh/Desktop/bee_sociality_dist")
# #setwd("/Users/tvasc/Desktop/bee_sociality_dist")
# source("00_utility_functions.R") # provides getModelTable() and getModelAvgParams()?
#

# #-------------------------------------------------------------------------------
# # Load hOUwie model results (use comments to select BIO4 or BIO1)
# #-------------------------------------------------------------------------------
# all_model_results <- list.files("houwie_results") # List result files
#
# all_model_results <- subset(all_model_results, grepl("_run2", all_model_results)) # Keep BIO1 only
# #all_model_results <- subset(all_model_results, grepl("_bio4", all_model_results))  # Keep BIO4 only
#
# model_names <- gsub(".Rsave", "", all_model_results) # Strip extensions for naming
#
# #-------------------------------------------------------------------------------
# # Load each result into a list named by model
# #-------------------------------------------------------------------------------
# all_results <- list() # Initialize
#
# for (i in seq_along(all_model_results)) {
#   load(paste0("houwie_results/", all_model_results[i])) # Loads object named 'res'
#   if (exists("res")) {
#     all_results[[i]] <- res
#     names(all_results)[i] <- model_names[i]
#   }
#   remove(res)
# }
#
# #-------------------------------------------------------------------------------
# # Generate AIC model comparison table
# #-------------------------------------------------------------------------------
# model_table <- getModelTable(all_results, type = "AICc")
# print(model_table)
#
# write.csv(model_table, file="model_table_results_bio1.csv") # Keep for BIO1
# #write.csv(model_table, file = "model_table_results_bio4.csv") # Keep for BIO4
#
#
# # 7/22/25
# # Model comparison table summarizing the performance of your six hOUwie models for temperature seasonality (BIO4)
# #                       np     lnLik   DiscLik   ContLik     AICc     dAICc        AICcwt
# # bm1_8states_run_bio4  12 -3665.862 -649.4574 -3011.799 7355.808 749.88796 1.458472e-163
# # ou1_8states_bio4      13 -3314.317 -657.8155 -2651.897 6654.733  48.81264  2.514601e-11
# # oum_cid_8states_bio4  14 -3325.218 -663.4269 -2657.396 6678.550  72.62935  1.693313e-16
# # oum_full_8states_bio4 16 -3286.887 -660.7572 -2621.426 6605.920   0.00000  1.000000e+00
# # oum_nest_8states_bio4 14 -3317.567 -662.8891 -2649.841 6663.247  57.32636  3.562379e-13
# # oum_soc_8states_bio4  14 -3307.305 -666.6537 -2634.220 6642.723  36.80296  1.019387e-08
#
# # Best model: oum_full_8states_bio4
# #   - Each trait combination has its own optimum
# #   - AICc = 6605.920, ΔAICc = 0, AICc weight = 1.0
# # Supports hypothesis that trait combinations (sociality × nesting) are tied to distinct climatic niches
# # Simpler models (e.g., BM1, OU1, sociality-only, nesting-only) fit worse
#
# #-------------------------------------------------------------------------------
# # Get model-averaged parameter estimates for BIO4 (temperature seasonality)
# # These reflect expected means and variances across species under the best-fitting models
# #-------------------------------------------------------------------------------
# average_pars <- getModelAvgParams(all_results, type = "AICc")
# head(average_pars)
# range(average_pars$expected_mean)
#     # 5.166410 <--> 6.467966 - BIO4
#
# write.csv(average_pars, file = "average_pars_bio1.csv") # Keep for BIO1
# #write.csv(average_pars, file="average_pars_bio4.csv") # Keep for BIO4
#
# # These are the expected trait values at the tips (species), averaged across models
# # based on AICc weights. They account for the full hOUwie model:
# # - The history of discrete trait evolution (e.g., sociality × nesting)
# # - The process of continuous trait evolution (e.g., optima, variance, strength of selection)
# # The values are in log scale and need to be transformed back to their original units
#
# #-------------------------------------------------------------------------------
# # Visualize expected trait values by tip state
# #-------------------------------------------------------------------------------
# # BIO4: Convert log-transformed values back to original scale (standard deviation × 100)
# # average_pars$expected_mean <- exp(average_pars$expected_mean)
# # average_pars$expected_mean <- round(average_pars$expected_mean, 3)
# # head(average_pars$expected_mean)
# # range(average_pars$expected_mean)
# #     # 175.284 <--> 644.172 - BIO4
#
# # BIO1: Convert Kelvins back to original scale (°C)
# average_pars$expected_mean <- exp(average_pars$expected_mean) - 273
# average_pars$expected_mean <- round(average_pars$expected_mean, 2)
# head(average_pars$expected_mean)
# range(average_pars$expected_mean)
#
# # Boxplots of expected values and variances across tip states (BIO1)
# boxplot(average_pars$expected_mean ~ average_pars$tip_state,
#         xlab = "Tip state", ylab = "Expected mean (temperature)")
#
# boxplot(average_pars$expected_var ~ average_pars$tip_state,
#         xlab = "Tip state", ylab = "Expected variance (temperature)")
#
# # Boxplots of expected values and variances across tip states (BIO4)
# # boxplot(average_pars$expected_mean ~ average_pars$tip_state,
# #         xlab = "Tip state", ylab = "Expected mean (temperature seasonality)")
# #
# # boxplot(average_pars$expected_var ~ average_pars$tip_state,
# #         xlab = "Tip state", ylab = "Expected variance (temperature seasonality)")
