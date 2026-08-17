# ==============================================================================
# Heatmaps
# ==============================================================================
# Plots heatmaps of the Americas showing proportion of social and above-ground species
# And makes scatterplots of absolute latitude vs. proportion social/above-ground nesting
# ==============================================================================
# Setup
# ------------------------------------------------------------------------------
rm(list = ls())

data_dir <- "/home/lenarh/data/bee_nesting_sociality"
curated_data_dir <- file.path(data_dir, "curated_data")
results_dir <- file.path(data_dir, "results", "heatmaps")
polygons_dir <- file.path(
  data_dir,
  "original_data",
  "TWDG",
  "wgsrpd-master",
  "level3",
  "level3.shp"
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(sp)
library(sf)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(data.table)
library(patchwork)
library(lwgeom)
library(nlme)

# Load polygon shapefiles
twdg_data <- st_read(polygons_dir) %>% as("Spatial")
proj4string(twdg_data) <- CRS("")

# ------------------------------------------------------------------------------
# Load trait and occurrence data
# ------------------------------------------------------------------------------
# Load trait data
traits <- read.csv(file.path(curated_data_dir, "bee_traits_clean.csv"))
colnames(traits)[3] <- "species"
traits$species <- gsub("_", " ", traits$species)

# Load thinned GBIF occurrence points
load(file.path(curated_data_dir, "thinned_points_res1.Rsave"))
colnames(thinned_points) <- c("species", "lat", "lon")
thinned_points$species <- gsub("_", " ", thinned_points$species)
head(thinned_points)

# ------------------------------------------------------------------------------
# Calculate species richness per region (as defined by TWDG polyons)
# ------------------------------------------------------------------------------
# Load function
organize.bubble.plot2 <- function(points, twdg_data) {
  focal_areas <- as.character(twdg_data$LEVEL3_COD) # Names of areas in global map shapefile
  species <- points[, 1] # Vector of species names
  points <- points[, c(3, 2)] # Vector of lat and long
  sp::coordinates(points) <- ~ lon + lat # Transforming into sp format coordinates
  results <- matrix(nrow = 0, ncol = 4) # Create empty results matrix
  list_result1 <- list()
  for (i in 1:length(focal_areas)) {
    one_area <- focal_areas[i]
    area_plus_buffer <- twdg_data[
      which(as.character(twdg_data$LEVEL3_COD) %in% one_area),
    ]
    if (nrow(area_plus_buffer) > 0) {
      res <- sp::over(points, area_plus_buffer)
      if (any(!is.na(res$LEVEL1_COD))) {
        # Should we use LEVEL3?
        sp_rich <- unique(species[which(!is.na(res$LEVEL1_COD))])
        list_result1[[i]] <- sp_rich
        names(list_result1)[i] <- one_area
        n_points <- length(sp_rich)
        coords <- sp::coordinates(area_plus_buffer)
        centroid_x <- mean(coords[, 1])
        centroid_y <- mean(coords[, 2])
        centroids <- sp::SpatialPoints(matrix(
          c(centroid_x, centroid_y),
          ncol = 2
        )) # Making centroid
        lon <- raster::extent(centroids)[1]
        lat <- raster::extent(centroids)[3]
        results <- rbind(results, cbind(n_points, one_area, lon, lat)) # Adds row for each iteration that is species richness in one_area for that iteration
      }
      cat(i, "\r")
    }
  }
  results <- as.data.frame(results)
  results$n_points <- as.numeric(results$n_points)
  results$lon <- as.numeric(results$lon)
  results$lat <- as.numeric(results$lat)
  save(list_result1, file = "list_result1.Rsave")
  return(results)
}

# Generate species richness per region using organize.bubble.plot() and thinned GBIF points
richness_per_area <- organize.bubble.plot2(points = thinned_points, twdg_data) # 369
colnames(richness_per_area)[1] <- "spp_rich"
head(richness_per_area)

# ------------------------------------------------------------------------------
# Compute proportion social and above-ground nesting per polygon
# ------------------------------------------------------------------------------
# Social
social_subset <- thinned_points[
  thinned_points$species %in%
    traits$species[traits$sociality_binary == "social"],
]
richness_social_per_area <- organize.bubble.plot2(social_subset, twdg_data)
colnames(richness_social_per_area)[1] <- "social_rich"

# Nesting
abvgrnd_subset <- thinned_points[
  thinned_points$species %in%
    traits$species[traits$nest_binary == "aboveground"],
]
richness_abvgrnd_per_area <- organize.bubble.plot2(abvgrnd_subset, twdg_data)
colnames(richness_abvgrnd_per_area)[1] <- "aboveground_rich"

# Merge with total richness
all_rich_social <- richness_per_area %>%
  left_join(
    richness_social_per_area,
    by = "one_area"
  ) %>%
  mutate(
    social_rich = ifelse(is.na(social_rich), 0, social_rich),
    prop_social = social_rich / spp_rich
  ) %>%
  select(
    one_area,
    prop_social
  )

all_rich_nest <- richness_per_area %>%
  left_join(
    richness_abvgrnd_per_area,
    by = "one_area"
  ) %>%
  mutate(
    aboveground_rich = ifelse(is.na(aboveground_rich), 0, aboveground_rich),
    prop_aboveground = aboveground_rich / spp_rich
  ) %>%
  select(
    one_area,
    prop_aboveground
  )

# ------------------------------------------------------------------------------
# Join trait proportion data to spatial polygons
# ------------------------------------------------------------------------------
twdg_data_sf <- st_as_sf(twdg_data)
st_crs(twdg_data_sf) <- 4326

# Join polygons with species richness
bee_twdg <- merge(
  twdg_data_sf,
  richness_per_area,
  by.x = "LEVEL3_COD",
  by.y = "one_area"
) %>%
  filter(LEVEL1_COD %in% c(7, 8))

# Join polygons with sociality proportion
bee_twdg_sociality <- merge(
  twdg_data_sf,
  all_rich_social,
  by.x = "LEVEL3_COD",
  by.y = "one_area"
) %>%
  filter(LEVEL1_COD %in% c(7, 8)) %>%
  st_as_sf()

# Join polygons with nesting proportion
bee_twdg_nest <- merge(
  twdg_data_sf,
  all_rich_nest,
  by.x = "LEVEL3_COD",
  by.y = "one_area"
) %>%
  filter(LEVEL1_COD %in% c(7, 8)) %>%
  st_as_sf()

# See which regions are in bee_twdg but not in bee_twdg_sociality
anti_join(
  bee_twdg %>% st_drop_geometry(),
  bee_twdg_sociality %>% st_drop_geometry(),
  by = "LEVEL3_COD"
) %>%
  select(LEVEL3_COD, LEVEL3_NAM)

# ------------------------------------------------------------------------------
# Plot heatmaps
# ------------------------------------------------------------------------------
plot_heatmap <- function(data, var, title, viridis_option = "viridis") {
  ggplot(data) +
    geom_sf(aes(fill = !!sym(var))) +
    labs(x = "Longitude", y = "Latitude", fill = title) +
    scale_fill_viridis_c(
      option = viridis_option,
      limits = c(0, 1)
    ) +
    coord_sf(
      xlim = c(-175, 0),
      ylim = c(-60, 90),
      expand = FALSE
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.95, 0.80),
      legend.justification = "left",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key.height = unit(0.6, "cm"),
      legend.key.width = unit(0.6, "cm"),
      legend.title = element_text(size = 12, margin = margin(b = 6)),
      legend.text = element_text(size = 10),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.background = element_blank(),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line = element_line(color = "black")
    )
}

prop_social_heatmap <- plot_heatmap(
  bee_twdg_sociality,
  "prop_social",
  "Proportion Social",
  viridis_option = "cividis"
)
prop_aboveground_heatmap <- plot_heatmap(
  bee_twdg_nest,
  "prop_aboveground",
  "Proportion\nAbove-Ground\nNesting",
  viridis_option = "cividis"
)

heatmaps <- prop_aboveground_heatmap +
  prop_social_heatmap +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 16, face = "bold"))

print(heatmaps)

ggsave(
  file.path(results_dir, "heatmaps.pdf"),
  heatmaps,
  width = 12,
  height = 12,
  dpi = 1500
)

ggsave(
  file.path(results_dir, "heatmaps.png"),
  heatmaps,
  width = 12,
  height = 12,
  dpi = 1500
)

# ------------------------------------------------------------------------------
# Scatterplots: abs(latitude) vs trait proportions with spatial autocorrelation correction
# ------------------------------------------------------------------------------
# Drop geometry from both heatmap sf objects to get raw trait values
social_df <- bee_twdg_sociality %>%
  st_drop_geometry() %>%
  select(one_area = LEVEL3_COD, prop_social)
social_df

nesting_df <- bee_twdg_nest %>%
  st_drop_geometry() %>%
  select(one_area = LEVEL3_COD, prop_aboveground)
nesting_df

# Repair invalid geometries first
bee_twdg_valid <- st_make_valid(bee_twdg)

# Then compute centroids
centroid_coords <- bee_twdg_valid %>%
  st_centroid() %>%
  st_coordinates() %>%
  as.data.frame()

# Combine with region ID and species richness
richness_clean <- bee_twdg_valid %>%
  st_drop_geometry() %>%
  select(one_area = LEVEL3_COD, spp_rich) %>%
  bind_cols(centroid_coords) %>%
  rename(lon = X, lat = Y)

head(richness_clean)

# Merge into dataframe for plotting
scatterplot_data <- social_df %>%
  inner_join(nesting_df, by = "one_area") %>%
  inner_join(richness_clean, by = "one_area")

head(scatterplot_data)
colnames(scatterplot_data)
nrow(scatterplot_data) # 113

# Calculate absolute latitude
scatterplot_data$abs_lat <- abs(scatterplot_data$lat)

# Remove NAs and inspect
scatterplot_data_clean <- scatterplot_data %>%
  filter(
    !is.na(abs_lat),
    !is.na(prop_social),
    !is.na(prop_aboveground),
    !is.na(lon),
    !is.na(lat)
  )

nrow(scatterplot_data_clean) # 113

scatterplot_data_clean %>% # Check which polygons have prop_social == 1
  filter(prop_social == 1) # ALU, ARU, BER, GNL

bee_twdg %>% # Get country names for polygons which prop_social == 1
  filter(LEVEL3_COD %in% c("ALU", "ARU", "BER", "GNL")) %>%
  select(LEVEL3_COD, LEVEL3_NAM)

scatterplot_data_clean %>% # Check which polygons have prop_aboveground == 1
  filter(prop_aboveground == 1) # ARU, BER, NLA

bee_twdg %>% # Get country names for polygons which prop_aboveground == 1
  filter(LEVEL3_COD %in% c("ARU", "BER", "NLA")) %>%
  select(LEVEL3_COD, LEVEL3_NAM)

# Fit GLS models using cleaned dataset (quadratic, since data has multiple modes)
gls_social <- gls(
  prop_social ~ abs_lat + I(abs_lat^2),
  data = scatterplot_data_clean,
  correlation = corSpher(form = ~ lon + lat, nugget = TRUE),
  method = "REML"
)

gls_aboveground <- gls(
  prop_aboveground ~ abs_lat + I(abs_lat^2),
  data = scatterplot_data_clean,
  correlation = corSpher(form = ~ lon + lat, nugget = TRUE),
  method = "REML"
)

summary(gls_social)
summary(gls_aboveground)

# Extract predicted values from GLS models
scatterplot_data_clean$gls_fit_social <- predict(gls_social)
scatterplot_data_clean$gls_fit_aboveground <- predict(gls_aboveground)

# ------------------------------------------------------------------------------
# Plot scatterplots
# ------------------------------------------------------------------------------
plot_trait_scatter <- function(
  data,
  trait_var,
  color_var,
  y_label,
  viridis_option = "viridis",
  show_size_legend = TRUE,
  gls_fit_var = NULL
) {
  p <- ggplot(data, aes(x = abs_lat, y = !!sym(trait_var))) +
    geom_point(aes(size = spp_rich, color = !!sym(color_var)), alpha = 0.8) +
    # Loess line
    geom_smooth(
      mapping = aes(),
      se = FALSE,
      formula = y ~ x,
      method = "loess",
      span = 1.5,
      color = scales::alpha("darkgray", 0.8),
      linewidth = 1
    )
  # Include line from GLS models
  if (!is.null(gls_fit_var)) {
    p <- p +
      geom_line(
        aes(y = !!sym(gls_fit_var)),
        color = scales::alpha("gray", 0.8),
        linetype = "dashed",
        linewidth = 1
      )
  }

  p +
    # Color scale (no legend shown by default)
    scale_color_viridis_c(
      option = viridis_option,
      limits = c(0, 1),
      guide = "none"
    ) +
    # Size scale
    scale_size_continuous(
      name = if (show_size_legend) "Species\nRichness" else NULL,
      breaks = if (show_size_legend) c(5, 50, 150, 300, 500, 700) else waiver(),
      limits = range(data$spp_rich, na.rm = TRUE),
      guide = if (show_size_legend) "legend" else "none"
    ) +
    # Labels
    labs(x = "Absolute Latitude", y = y_label) +
    # Theme
    theme_classic() +
    theme(
      aspect.ratio = 1,
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.key.height = unit(0.3, "in"),
      legend.key.width = unit(0.3, "in")
    )
}

social_scatter <- plot_trait_scatter(
  data = scatterplot_data_clean,
  trait_var = "prop_social",
  color_var = "prop_social",
  y_label = "Proportion Social",
  viridis_option = "cividis",
  show_size_legend = TRUE,
  gls_fit_var = "gls_fit_social"
) +
  theme(
    legend.position = c(1.3, 0.7),
    legend.justification = "right",
    legend.background = element_rect(fill = "white", color = NA)
  )

nesting_scatter <- plot_trait_scatter(
  data = scatterplot_data_clean,
  trait_var = "prop_aboveground",
  color_var = "prop_aboveground",
  y_label = "Proportion Above-Ground Nesting",
  viridis_option = "cividis",
  show_size_legend = FALSE,
  gls_fit_var = "gls_fit_aboveground"
)

combined_scatterplots <- social_scatter +
  nesting_scatter +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

print(combined_scatterplots)

ggsave(
  filename = file.path(results_dir, "combined_scatterplots.pdf"),
  plot = combined_scatterplots,
  width = 10,
  height = 5,
  dpi = 1500
)

ggsave(
  filename = file.path(results_dir, "combined_scatterplots.png"),
  plot = combined_scatterplots,
  width = 10,
  height = 5,
  dpi = 1500
)

# ------------------------------------------------------------------------------
# Assemble four-panel figure of heatmaps + scatterplots
# ------------------------------------------------------------------------------
# Add tags and formatting manually to each plot
prop_social_heatmap <- prop_social_heatmap +
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

social_scatter <- social_scatter +
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

prop_aboveground_heatmap <- prop_aboveground_heatmap +
  labs(tag = "C") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

nesting_scatter <- nesting_scatter +
  labs(tag = "D") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

# Combine plots into a 2x2 grid
combined_mapscatter <- ((prop_social_heatmap +
  social_scatter +
  plot_layout(widths = c(1, 1))) /
  (prop_aboveground_heatmap + nesting_scatter + plot_layout(widths = c(1, 1))))

print(combined_mapscatter)

ggsave(
  filename = file.path(results_dir, "combined_mapscatter.png"),
  plot = combined_mapscatter,
  width = 12,
  height = 8,
  dpi = 1500
)

ggsave(
  filename = file.path(results_dir, "combined_mapscatter.pdf"),
  plot = combined_mapscatter,
  width = 12,
  height = 8,
  dpi = 1500
)
