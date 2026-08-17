#===============================================================================
# Plotting points to inspect species distributions
#===============================================================================
# Setup
#-------------------------------------------------------------------------------
rm(list = ls())

data_dir <- "/home/lenarh/data/bee_nesting_sociality"
results_dir <- file.path(data_dir, "results", "supplemental")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(sf)
library(rnaturalearth)

# Load world polygons
world <- ne_countries(scale = "medium", returnclass = "sf")

#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------
load(file.path(data_dir, "curated_data", "thinned_points_res1.Rsave"))

species <- unique(thinned_points$species)
species <- species[species != ""]

#-------------------------------------------------------------------------------
# Plot species distributions
#-------------------------------------------------------------------------------
pdf(file.path(results_dir, "per_species_maps.pdf"))

for (spp_index in seq_along(species)) {
  tmp_subset <- thinned_points[
    thinned_points$species == species[spp_index],
  ]

  plot(st_geometry(world), col = "lightgray", border = "gray")

  points(
    tmp_subset$lon,
    tmp_subset$lat,
    col = "red",
    pch = 16,
    cex = 0.5
  )

  title(species[spp_index])

  cat(spp_index, "\r")
}

dev.off()

#-------------------------------------------------------------------------------
# Done!
#-------------------------------------------------------------------------------
