#===============================================================================
# Occurrence data curation, thinning, and climate summaries
#===============================================================================
# Workflow:
# PART 1: Create GeoPackages for species introduced outside their native ranges.
# PART 2: Re-import the manually reviewed GeoPackages, remove flagged occurrence
#         points and all Apis mellifera records, and create a supplemental file
#         containing the manually removed points.
# PART 3: Restrict occurrences to species in the phylogeny, save the fully
#         curated pre-thinning occurrence dataset, spatially thin occurrence
#         records and save for use in downstream analyses, extract climate values,
#         and calculate per-species climate summary statistics.
#===============================================================================
# Setup
#-------------------------------------------------------------------------------
rm(list = ls())

data_dir <- "/home/lenarh/data/bee_nesting_sociality"
repo_dir <- "/home/lenarh/repos/bee_nesting_sociality"

raw_data_dir <- file.path(data_dir, "original_data")
curated_data_dir <- file.path(data_dir, "curated_data")
results_dir <- file.path(data_dir, "results")

gpkg_dir <- file.path(raw_data_dir, "invasives_gpkg")
climate_dir <- file.path(raw_data_dir, "climate_layers")
climate_summary_dir <- file.path(curated_data_dir, "climate_summaries")
supplemental_dir <- file.path(results_dir, "supplemental")

dir.create(curated_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gpkg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplemental_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(climate_summary_dir, recursive = TRUE, showWarnings = FALSE)

library(data.table)
library(sf)
library(ape)
library(phytools)
library(raster)
library(sp)
library(rworldmap)

source(file.path(repo_dir, "00_utility_functions.R"))

# Custom functions used in Part 3:
# - Thinning(): spatially thin occurrence records per species
# - DataFromPoints(): extract raster values at occurrence coordinates
# - GetClimateSummStats_custom(): calculate per-species climate summaries

#===============================================================================
# PART 1: Make QGIS files for species with introductions outside native ranges
#===============================================================================

#-------------------------------------------------------------------------------
# Load occurrence data and invasive-species list
#-------------------------------------------------------------------------------
dorey_dat <- fread(
  file.path(raw_data_dir, "bee_occurrences_Dorey.csv")
)

invasives <- fread(
  file.path(raw_data_dir, "invasive_species.csv")
)

setDT(dorey_dat)
setDT(invasives)

# Assign a unique ID to every row in the original Dorey dataset.
# These IDs will later identify occurrence points selected for removal in QGIS.
dorey_dat[, point_id := .I]

#-------------------------------------------------------------------------------
# Prepare invasive-species list
#-------------------------------------------------------------------------------
inv <- invasives[
  Introduction_type != "N" &
    !is.na(Non_native_species)
]

# Remove leading and trailing whitespace from species names
inv[, species := trimws(Non_native_species)]

# Keep species present in both the invasive-species list and Dorey dataset
species_keep <- intersect(
  inv$species,
  dorey_dat$species
)

species_list <- sort(unique(species_keep))

#-------------------------------------------------------------------------------
# Build dataset containing all occurrence points for introduced species
#-------------------------------------------------------------------------------
all_points <- dorey_dat[
  species %in% species_keep,
  .(
    species,
    decimalLongitude,
    decimalLatitude,
    point_id
  )
]

# Add introduction metadata
all_points <- merge(
  all_points,
  inv[, .(species, From, Found_in)],
  by = "species",
  all.x = TRUE
)

# Remove missing or invalid coordinates
all_points <- all_points[
  !is.na(decimalLongitude) &
    !is.na(decimalLatitude) &
    decimalLongitude >= -180 &
    decimalLongitude <= 180 &
    decimalLatitude >= -90 &
    decimalLatitude <= 90
]

# Add the column that will be edited manually in QGIS
all_points[, remove_flag := FALSE]

#-------------------------------------------------------------------------------
# Export one GeoPackage per species for inspection in QGIS
#
# Uncomment this section when GeoPackages need to be generated or regenerated.
#-------------------------------------------------------------------------------
# species_list <- sort(unique(all_points$species))
#
# for (sp in species_list) {
#
#   message("Exporting: ", sp)
#
#   sp_dat <- all_points[species == sp]
#
#   if (nrow(sp_dat) == 0) {
#     next
#   }
#
#   sp_sf <- st_as_sf(
#     sp_dat,
#     coords = c("decimalLongitude", "decimalLatitude"),
#     crs = 4326,
#     remove = FALSE
#   )
#
#   file_name <- gsub("[^A-Za-z0-9]", "_", sp)
#
#   out_file <- file.path(
#     gpkg_dir,
#     paste0(file_name, ".gpkg")
#   )
#
#   st_write(
#     sp_sf,
#     out_file,
#     delete_dsn = TRUE,
#     quiet = TRUE
#   )
# }

#-------------------------------------------------------------------------------
# Optional sanity check for a manually edited GeoPackage
#-------------------------------------------------------------------------------
# species_to_check <- "Megachile_rotundata"
#
# gpkg_path <- file.path(
#   gpkg_dir,
#   paste0(species_to_check, ".gpkg")
# )
#
# occ_sf <- st_read(
#   dsn = gpkg_path,
#   quiet = TRUE
# )
#
# occ_df <- st_drop_geometry(occ_sf)
#
# n_remove_true <- sum(
#   occ_df$remove_flag == TRUE,
#   na.rm = TRUE
# )
#
# n_remove_true

#===============================================================================
# STOP HERE
# Inspect and edit the GeoPackages in QGIS before running the remaining sections.
#
# When continuing later, rerun the script from the beginning. Because point_id is
# assigned from the original row order each time, the same original records will
# receive the same point_id values.
#===============================================================================

#===============================================================================
# PART 2: Remove points outside native ranges and all Apis mellifera records
#===============================================================================

#-------------------------------------------------------------------------------
# Re-import GeoPackages after manual flagging in QGIS
#-------------------------------------------------------------------------------
gpkg_files <- list.files(
  gpkg_dir,
  pattern = "\\.gpkg$",
  full.names = TRUE
)

gpkg_list <- lapply(
  gpkg_files,
  function(f) {
    message("Reading: ", basename(f))
    st_read(f, quiet = TRUE)
  }
)

# Combine occurrence data from all reviewed GeoPackages
reviewed_points <- rbindlist(
  lapply(
    gpkg_list,
    function(x) as.data.table(x)
  ),
  fill = TRUE
)

#-------------------------------------------------------------------------------
# Inspect flagged points
#-------------------------------------------------------------------------------
flagged_points <- reviewed_points[
  !is.na(remove_flag) &
    remove_flag %in% c(TRUE, "TRUE", 1)
]

flagged_df <- as.data.table(
  st_drop_geometry(flagged_points)
)

head(flagged_df)

nrow(flagged_df)
# Expected: 51,101 points flagged for removal

length(unique(flagged_df$species))

#-------------------------------------------------------------------------------
# Identify species with no flagged points
#-------------------------------------------------------------------------------
expected_species <- sort(unique(species_keep))

flagged_species <- sort(
  unique(flagged_df$species)
)

no_flag_species <- setdiff(
  expected_species,
  flagged_species
)

length(expected_species) # 75
length(flagged_species) # 66
length(no_flag_species) # 9

no_flag_species
# These species lacked introduced points in the Dorey dataset or did not meet
# the criteria for point removal.

#-------------------------------------------------------------------------------
# Prepare point IDs for removal
#-------------------------------------------------------------------------------
dorey_dat[, point_id := as.integer(point_id)]
flagged_df[, point_id := as.integer(point_id)]

# Confirm that none of the flagged records has a missing point ID
sum(is.na(flagged_df$point_id))

# Extract unique IDs to remove
points_to_remove <- unique(flagged_df$point_id)

length(points_to_remove)
# Expected: 51,101

# Confirm that every selected ID exists in the original Dorey dataset
missing_ids <- setdiff(
  points_to_remove,
  dorey_dat$point_id
)

length(missing_ids)
# Should be zero

#-------------------------------------------------------------------------------
# Remove manually flagged occurrence points
#-------------------------------------------------------------------------------
cleaned_dorey <- dorey_dat[
  !point_id %in% points_to_remove
]

# Remove all Apis mellifera occurrence points
cleaned_dorey <- cleaned_dorey[
  species != "Apis mellifera"
]

#-------------------------------------------------------------------------------
# Sanity checks
#-------------------------------------------------------------------------------
sum(cleaned_dorey$species == "Apis mellifera")
# Should be zero

n_removed_apis <- sum(
  dorey_dat$species == "Apis mellifera"
)

n_removed_apis
# Expected: 630,049 Apis mellifera points removed

n_before <- nrow(dorey_dat)

n_removed <- length(points_to_remove) + n_removed_apis

n_after <- nrow(cleaned_dorey)

n_before
# Expected: 6,785,927

n_removed
# Expected: 681,150 total removed

n_after
# Expected: 6,104,777

#-------------------------------------------------------------------------------
# Make supplemental file containing manually removed introduced points
#-------------------------------------------------------------------------------
removed_points_full <- dorey_dat[
  point_id %in% points_to_remove
]

fwrite(
  removed_points_full,
  file.path(
    supplemental_dir,
    "removed_points_with_metadata.csv"
  )
)

head(removed_points_full)

nrow(removed_points_full)
# Expected: 51,101

length(unique(removed_points_full$species))


#===============================================================================
# PART 3: Final occurrence filtering, thinning, and climate summaries
#===============================================================================

#-------------------------------------------------------------------------------
# Keep only occurrence points for species present in the phylogeny
#-------------------------------------------------------------------------------
bee_tree <- read.tree(
  file.path(
    curated_data_dir,
    "bee_tree_pruned.tre"
  )
)

# Convert occurrence species names to the underscore format used in the tree
cleaned_dorey[, species := gsub(" ", "_", species)]

# Restrict occurrences to species represented in the phylogeny
all_cleaned_points <- cleaned_dorey[
  species %in% bee_tree$tip.label
]

# point_id was only needed for identifying records selected in QGIS
all_cleaned_points[, point_id := NULL]

#-------------------------------------------------------------------------------
# Save fully curated occurrence dataset before spatial thinning
#-------------------------------------------------------------------------------
fwrite(
  all_cleaned_points,
  file.path(
    curated_data_dir,
    "bee_occurrences_clean.csv"
  )
)

#-------------------------------------------------------------------------------
# Spatially thin occurrence records
# Keep one occurrence point per grid cell per species
#-------------------------------------------------------------------------------
thinned_points <- Thinning(
  all_cleaned_points,
  species = "species",
  lat = "decimalLatitude",
  lon = "decimalLongitude",
  n = 1
)

colnames(thinned_points) <- c(
  "species",
  "lat",
  "lon"
)

# Save filtered and spatially thinned occurrence points
save(
  thinned_points,
  file = file.path(
    curated_data_dir,
    "thinned_points_res1.Rsave"
  )
)

nrow(thinned_points)
# Expected from the previous run: 185,921 occurrence points

# This thinned_points dataset will be used in downstream analyses.

#-------------------------------------------------------------------------------
# Load climate rasters
#-------------------------------------------------------------------------------
climate_files <- list.files(
  climate_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

# Use raster filenames without the .tif extension as layer labels
climate_labels <- tools::file_path_sans_ext(
  basename(climate_files)
)

# Read each climate file as a RasterLayer
all_layers <- lapply(
  climate_files,
  raster
)

names(all_layers) <- climate_labels

#-------------------------------------------------------------------------------
# Extract climate values and calculate per-species summaries
#-------------------------------------------------------------------------------
# For each climate layer:
#   1. Extract values at all thinned occurrence coordinates.
#   2. Calculate per-species summary statistics.
#   3. Save one summary-statistics CSV.

for (i in seq_along(all_layers)) {
  one_layer <- all_layers[[i]]
  one_label <- names(all_layers)[i]

  message("Processing: ", one_label)

  # Extract climate values for every species-coordinate row
  allpoints <- DataFromPoints(
    thinned_points,
    one_layer
  )

  # Calculate per-species climate summary statistics
  summstats <- GetClimateSummStats_custom(
    allpoints,
    type = "raw"
  )

  # Save one summary file per climate layer
  fwrite(
    summstats,
    file.path(
      climate_summary_dir,
      paste0(
        one_label,
        "_climate_summstats.csv"
      )
    )
  )

  message(one_label, " done.")
}

#===============================================================================
# Done!
#===============================================================================
