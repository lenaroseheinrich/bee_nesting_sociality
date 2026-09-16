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
results_dir <- file.path(
  data_dir,
  "results",
  "hOUwie",
  "Aug 31 Run (corHMM fixed)"
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(OUwie)
library(corHMM)
library(parallel)

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
  if (
    all(file.exists(
      csv_table_path,
      csv_pars_path,
      plot_mean_path,
      plot_var_path
    ))
  ) {
    message("Skipping ", bio, " (results already exist)")
    next
  }

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
  model_table <- OUwie::getModelTable(all_results, type = "AICc")
  print(model_table)
  write.csv(model_table, file = csv_table_path, row.names = TRUE)

  # Get model-averaged parameters
  average_pars <- OUwie::getModelAvgParams(all_results, type = "AICc")

  # Convert expected means
  if (bio == "bio_1") {
    average_pars$expected_mean <- exp(average_pars$expected_mean) - 273.15 # convert kelvin to °C
    ylab <- "Expected mean (°C)"
  } else if (bio == "bio_4") {
    average_pars$expected_mean <- exp(average_pars$expected_mean) / 10
    ylab <- "Expected mean (temperature seasonality)"
  } else if (bio == "bio_12") {
    average_pars$expected_mean <- exp(average_pars$expected_mean)
    ylab <- "Expected mean (precipitation)"
  } else if (bio == "bio_15") {
    average_pars$expected_mean <- exp(average_pars$expected_mean)
    ylab <- "Expected mean (precipitation seasonality)"
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

bio1_summary <- bio1 %>%
  group_by(tip_state) %>%
  summarize(
    mean_temp = round(mean(expected_mean)),
    .groups = "drop"
  )

print(bio1_summary)

bio4_summary <- bio4 %>%
  group_by(tip_state) %>%
  summarize(
    mean_seasonality = round(mean(expected_mean)),
    .groups = "drop"
  )

print(bio4_summary)

bio12_summary <- bio12 %>%
  group_by(tip_state) %>%
  summarize(
    mean_precip = round(mean(expected_mean)),
    .groups = "drop"
  )

print(bio12_summary)

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

# Done!

#   tip_state            mean_bio1 mean_bio4 mean_bio12 mean_bio15
#   <fct>                    <dbl>     <dbl>      <dbl>      <dbl>
# 1 solitary_ground           16.3       441        604       555
# 2 solitary_aboveground      17.2       343        843       555
# 3 social_ground             13.2       499        907       469
# 4 social_aboveground        21.6       180       1388       472
