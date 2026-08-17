#===============================================================================
# Kruskal-Wallis test
#===============================================================================
# Runs Kruskal-Wallis tests between bee traits and climatic variables,
# and generates violin plots.
#===============================================================================
# Setup
#-------------------------------------------------------------------------------
rm(list = ls())

data_dir <- "/home/lenarh/data/bee_nesting_sociality"
results_dir <- file.path(data_dir, "results", "kruskal-wallis")
climate_summaries_dir <- file.path(
  data_dir,
  "curated_data",
  "climate_summaries"
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(phytools)
library(ggplot2)
library(gridExtra)
library(FSA)
library(multcompView)
library(viridis)
library(dunn.test)
library(colorspace)

#-------------------------------------------------------------------------------
# Load trait data, phylogeny, and list of climatic variable summary files
#-------------------------------------------------------------------------------
traits <- read.csv(
  file.path(data_dir, "curated_data", "bee_traits_clean.csv")
)

tree <- read.tree(
  file.path(data_dir, "curated_data", "bee_tree_pruned.tre")
)

all_climatic_vars <- list.files(
  # Selects summstats.csv files
  file.path(climate_summaries_dir),
  pattern = "summstats.csv",
  full.names = TRUE
)

# ------------------------------------------------------------------------------
# Kruskal-Wallis test for combined sociality-nesting strategies vs. all 19 climate variables
# ------------------------------------------------------------------------------
# trait combinations: solitary_ground, solitary_aboveground,
# social_ground, social_aboveground.
# ------------------------------------------------------------------------------
sink(file.path(results_dir, "kruskal_results_all_vars.txt"))

plot_list <- list()

# Trait levels and labels
trait_levels <- c(
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

# CHELSA variable names
bio_names <- c(
  bio01 = "Annual Mean Temperature",
  bio02 = "Mean Diurnal Range",
  bio03 = "Isothermality",
  bio04 = "Temperature Seasonality",
  bio05 = "Max Temperature of Warmest Month",
  bio06 = "Min Temperature of Coldest Month",
  bio07 = "Temperature Annual Range",
  bio08 = "Mean Temperature of Wettest Quarter",
  bio09 = "Mean Temperature of Driest Quarter",
  bio10 = "Mean Temperature of Warmest Quarter",
  bio11 = "Mean Temperature of Coldest Quarter",
  bio12 = "Annual Precipitation",
  bio13 = "Precipitation of Wettest Month",
  bio14 = "Precipitation of Driest Month",
  bio15 = "Precipitation Seasonality",
  bio16 = "Precipitation of Wettest Quarter",
  bio17 = "Precipitation of Driest Quarter",
  bio18 = "Precipitation of Warmest Quarter",
  bio19 = "Precipitation of Coldest Quarter"
)

# Colors
trait_colors <- viridis::viridis(
  n = 4,
  option = "cividis"
)

names(trait_colors) <- trait_levels

point_colors <- sapply(
  trait_colors,
  function(x) colorspace::darken(x, amount = 0.4)
)

for (climate_index in seq_along(all_climatic_vars)) {
  climate_file <- all_climatic_vars[climate_index]

  # Read climate data
  climate <- read.csv(
    file.path(climate_summaries_dir, basename(climate_file))
  )

  climate <- subset(
    climate,
    !is.na(climate[, 3])
  )

  # Match species
  sampled_species <- intersect(
    intersect(traits$tips, climate$species),
    tree$tip.label
  )

  subset_traits <- subset(
    traits,
    traits$tips %in% sampled_species
  )

  subset_climate <- subset(
    climate,
    climate$species %in% sampled_species
  )

  subset_tree <- keep.tip(
    tree,
    tree$tip.label[
      tree$tip.label %in% sampled_species
    ]
  )

  merged_table <- merge(
    subset_traits,
    subset_climate,
    by.x = "tips",
    by.y = "species"
  )

  # Combined trait group
  merged_table$comb_nest_soc <- interaction(
    merged_table$sociality_binary,
    merged_table$nest_binary,
    sep = "_"
  )

  merged_table$value <- merged_table[, 9]

  merged_table <- merged_table[
    !is.na(merged_table$value),
  ]

  merged_table$comb_nest_soc <- factor(
    merged_table$comb_nest_soc,
    levels = trait_levels
  )

  # Extract BIO variable name
  bio_id <- sub(
    ".*(bio[0-9]{2}).*",
    "\\1",
    basename(climate_file)
  )

  label <- bio_names[bio_id]

  # Kruskal-Wallis test
  kruskal_result <- kruskal.test(
    value ~ comb_nest_soc,
    data = merged_table
  )

  p_value <- kruskal_result$p.value

  significance <- ifelse(
    p_value < 0.0001,
    "****",
    ifelse(
      p_value < 0.001,
      "***",
      ifelse(
        p_value < 0.01,
        "**",
        ifelse(
          p_value < 0.05,
          "*",
          "ns"
        )
      )
    )
  )

  print(paste0("comb_nest_soc ~ ", label))
  print(kruskal_result)

  # Position significance stars
  y_max <- max(merged_table$value, na.rm = TRUE)
  y_min <- min(merged_table$value, na.rm = TRUE)

  y_range <- y_max - y_min

  star_y <- y_max + (0.08 * y_range)

  # Plot
  plot <- ggplot(
    merged_table,
    aes(x = comb_nest_soc, y = value, fill = comb_nest_soc)
  ) +
    geom_violin(width = 0.9, adjust = 0.8, alpha = 0.6, show.legend = FALSE) +
    geom_jitter(
      aes(color = comb_nest_soc),
      width = 0.15,
      alpha = 0.7,
      size = 0.5,
      show.legend = FALSE
    ) +
    annotate(
      "text",
      x = 2.5,
      y = star_y,
      label = significance,
      size = 5,
      fontface = "bold"
    ) +
    scale_fill_manual(values = trait_colors) +
    scale_color_manual(values = point_colors) +
    scale_x_discrete(labels = trait_labels) +
    scale_y_continuous(
      breaks = pretty(merged_table$value, n = 6),
      expand = expansion(mult = c(0.05, 0.15))
    ) +
    labs(x = "", y = "Raster units", title = label) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.ticks.y = element_line(color = "black", linewidth = 0.5)
    )

  print(plot)

  plot_list[[climate_index]] <- plot
}

# Drop "all" climate plot if present
drop_plot <- which(
  tolower(
    sub(
      "_climate_summstats\\.csv$",
      "",
      basename(all_climatic_vars)
    )
  ) ==
    "all"
)

if (length(drop_plot)) {
  plot_list <- plot_list[-drop_plot]
}

# Arrange plots
n_plots <- length(plot_list)

ncols <- 3
nrows <- ceiling(n_plots / ncols)

combined_plot <- grid.arrange(
  grobs = plot_list,
  ncol = ncols,
  nrow = nrows
)

grid.arrange(
  grobs = plot_list,
  ncol = ncols,
  nrow = nrows
)

# Save plots
ggsave(
  file.path(results_dir, "violins_all_vars.pdf"),
  plot = combined_plot,
  width = 15,
  height = 40
)

sink()

# ------------------------------------------------------------------------------
# Kruskal-Wallis test for combined sociality-nesting strategies vs. focal climate variables
# ------------------------------------------------------------------------------
# trait combinations: solitary_ground, solitary_aboveground,
# social_ground, social_aboveground.
# ------------------------------------------------------------------------------
# bio_1, bio_4, bio_12, bio_15
# ------------------------------------------------------------------------------
sink(file.path(results_dir, "kruskal_results_focal_vars.txt"))

# Select variables
selected_vars <- c(
  "CHELSA_bio01_1981-2010_V.2.1_climate_summstats.csv",
  "CHELSA_bio04_1981-2010_V.2.1_climate_summstats.csv",
  "CHELSA_bio12_1981-2010_V.2.1_climate_summstats.csv",
  "CHELSA_bio15_1981-2010_V.2.1_climate_summstats.csv"
)

# Variable labels
variable_labels <- c(
  "CHELSA_bio01_1981-2010_V.2.1_climate_summstats.csv" = "Mean annual temperature",

  "CHELSA_bio04_1981-2010_V.2.1_climate_summstats.csv" = "Temperature seasonality",

  "CHELSA_bio12_1981-2010_V.2.1_climate_summstats.csv" = "Annual precipitation",

  "CHELSA_bio15_1981-2010_V.2.1_climate_summstats.csv" = "Precipitation seasonality"
)

variable_units <- c(
  "CHELSA_bio01_1981-2010_V.2.1_climate_summstats.csv" = "Mean annual temperature (°C)",

  "CHELSA_bio04_1981-2010_V.2.1_climate_summstats.csv" = "Temperature seasonality (SD x 100)",

  "CHELSA_bio12_1981-2010_V.2.1_climate_summstats.csv" = "Annual precipitation (mm)",

  "CHELSA_bio15_1981-2010_V.2.1_climate_summstats.csv" = "Precipitation seasonality (% CV)"
)

# Convert CHELSA raster values to interpretable units
convert_climate_units <- function(values, climate_file) {
  if (grepl("bio01", climate_file)) {
    # Tenths of Kelvin to degrees Celsius
    values <- values / 10 - 273.15
  } else if (grepl("bio04", climate_file)) {
    # Tenths of SD × 100 to SD × 100
    values <- values / 10
  } else if (grepl("bio12", climate_file)) {
    # Already expressed in millimetres
    values <- values
  } else if (grepl("bio15", climate_file)) {
    # Tenths of percent CV to percent CV
    values <- values / 10
  }

  return(values)
}

# Trait order and labels
trait_levels <- c(
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

# Colors
violin_colors <- viridis::viridis(
  n = 4,
  option = "cividis"
)

names(violin_colors) <- trait_levels

point_colors <- sapply(
  violin_colors,
  function(x) darken(x, amount = 0.4)
)

# Initialize plot storage
plot_list <- list()

# Compute global min/max for each climate variable
climate_ranges <- list()

for (file in selected_vars) {
  climate <- read.csv(
    file.path(climate_summaries_dir, file)
  )

  climate_values <- convert_climate_units(
    climate[, 3],
    file
  )

  climate_ranges[[file]] <- c(
    min(climate_values, na.rm = TRUE),
    max(climate_values, na.rm = TRUE)
  )
}

for (climate_index in seq_along(selected_vars)) {
  climate_file <- selected_vars[climate_index]

  climate <- read.csv(
    file.path(climate_summaries_dir, climate_file)
  )

  climate <- subset(
    climate,
    !is.na(climate[, 3])
  )

  sampled_species <- intersect(
    intersect(traits$tips, climate$species),
    tree$tip.label
  )

  subset_traits <- subset(
    traits,
    traits$tips %in% sampled_species
  )

  subset_climate <- subset(
    climate,
    climate$species %in% sampled_species
  )

  subset_tree <- keep.tip(
    tree,
    tree$tip.label[
      tree$tip.label %in% sampled_species
    ]
  )

  merged_table <- merge(
    subset_traits,
    subset_climate,
    by.x = "tips",
    by.y = "species"
  )

  # Combined trait group
  merged_table$comb_nest_soc <- interaction(
    merged_table$sociality_binary,
    merged_table$nest_binary,
    sep = "_"
  )

  merged_table$value <- convert_climate_units(
    merged_table[, 9],
    climate_file
  )

  merged_table <- merged_table[
    !is.na(merged_table$value),
  ]

  merged_table$comb_nest_soc <- factor(
    merged_table$comb_nest_soc,
    levels = trait_levels
  )

  # Kruskal-Wallis tests
  kruskal_result <- kruskal.test(
    value ~ comb_nest_soc,
    data = merged_table
  )

  p_value <- kruskal_result$p.value

  label <- variable_labels[climate_file]

  print(paste0("comb_nest_soc ~ ", label))
  print(kruskal_result)

  if (p_value < 0.05) {
    dunn_result <- dunnTest(
      value ~ comb_nest_soc,
      data = merged_table,
      method = "bonferroni"
    )

    print("Dunn's post-hoc test (Bonferroni corrected):")
    print(dunn_result$res)
  }

  # Axis scaling
  x_limits <- climate_ranges[[climate_file]]

  x_min <- x_limits[1]
  x_max <- x_limits[2]

  pretty_breaks <- pretty(
    x_limits,
    n = 6
  )

  tick_step <- pretty_breaks[2] - pretty_breaks[1]

  tick_min <- floor(x_min / tick_step) * tick_step
  tick_max <- ceiling(x_max / tick_step) * tick_step

  tick_breaks <- seq(
    tick_min,
    tick_max,
    by = tick_step
  )

  # Plot
  plot <- ggplot(
    merged_table,
    aes(
      x = value,
      y = comb_nest_soc,
      fill = comb_nest_soc
    )
  ) +

    geom_violin(
      width = 0.9,
      adjust = 0.8,
      alpha = 0.6,
      show.legend = FALSE
    ) +

    geom_jitter(
      aes(color = comb_nest_soc),
      height = 0.2,
      alpha = 0.7,
      size = 0.5,
      show.legend = FALSE
    ) +

    scale_fill_manual(
      values = violin_colors
    ) +

    scale_color_manual(
      values = point_colors
    ) +

    scale_y_discrete(
      labels = trait_labels
    ) +

    scale_x_continuous(
      limits = c(tick_min, tick_max),
      breaks = tick_breaks
    ) +

    labs(y = "", x = variable_units[climate_file], title = label) +

    theme_classic() +

    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.y = element_text(size = 24, color = "black"),
      axis.text.x = element_text(size = 20, color = "black"),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
      plot.margin = margin(t = 10, r = 40, b = 10, l = 40)
    )

  print(plot)

  plot_list[[climate_index]] <- plot
}

# Arrange plots
ncols <- 1
nrows <- ceiling(length(plot_list) / ncols)

combined_plot <- grid.arrange(
  grobs = plot_list,
  ncol = ncols,
  nrow = nrows
)

# Save plots
ggsave(
  file.path(results_dir, "violins_focal_vars.pdf"),
  plot = combined_plot,
  width = 12,
  height = 24
)

sink()

#-------------------------------------------------------------------------------
# Done!
#-------------------------------------------------------------------------------
