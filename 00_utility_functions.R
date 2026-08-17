# Some utility functions

#' #' Taxize GBIF
#' #' @param name A character vector with species names.
#' resolveGBIF <- function(name) {
#'   gnr_resolve_x <- function(x) {
#'     sources <- taxize::gnr_datasources()
#'     tmp.name <- suppressWarnings(taxize::gnr_resolve(names=x, data_source_ids=sources$id[sources$title == "GBIF Backbone Taxonomy"], best_match_only=TRUE)$matched_name)
#'     if(is.null(tmp.name)) {
#'       tmp.name <- paste0("UNMATCHED_",x)
#'     }
#'     return(tmp.name)
#'   }
#'   new.names <- pbapply::pbsapply(name, gnr_resolve_x)
#'   return(as.character(new.names))
#' }
#'
#' #' Removes points in the sea
#' #' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' #' @param lat The name of the column in the data.frame with latitudes
#' #' @param lon The name of the column in the data.frame with longitudes
#' #' @param buffer A number of degrees around continental areas where points are still kept after filtering
#' RemoveSeaPoints <- function(points, lon="decimalLongitude", lat="decimalLatitude", buffer=0) {
#'   npoints_start <- nrow(points)
#'   tmp_points = points
#'   colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
#'   colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
#'   wrld_map <- rworldmap::getMap(resolution="low") # leaving both maps in the code for now, should probably drop one of them later
#'   coords <- tmp_points[,c("x","y")]
#'   sp::coordinates(coords) <- ~ x + y
#'   sp::proj4string(coords) <- sp::proj4string(wrld_map)
#'   country_plus_buffer <- raster::buffer(wrld_map, buffer) # adding buffer around polygons
#'   answer <- which(is.na(sp::over(coords, country_plus_buffer)))
#'   if(length(answer) > 0) {
#'     points <- points[-answer,]
#'     npoints_end <- nrow(points)
#'     print(paste0(npoints_start - npoints_end, " points removed."))
#'     return(points)
#'   } else {
#'     print("no points removed")
#'     return(points) }
#' }
#'
#' #' Removes points that have 0 for both latitude and longitude
#' #' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' #' @param lat The name of the column in the data.frame with latitudes
#' #' @param lon The name of the column in the data.frame with longitudes
#' RemoveZeros <- function(points, lon="decimalLongitude", lat="decimalLatitude") {
#'   if(!inherits(points, "data.frame")) {
#'     stop("Argument points is not a data.frame.")
#'   }
#'   npoints_start <- nrow(points)
#'   tmp_points = points
#'   colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
#'   colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
#'   if(any(tmp_points$x==0 & tmp_points$y==0)) {
#'     points <- points[-which(tmp_points$x==0 & tmp_points$y==0),]
#'   }
#'   npoints_end <- nrow(points)
#'   print(paste0(npoints_start - npoints_end, " points removed."))
#'   return(points)
#' }
#'
#' #' Removes outliers
#' #' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' #' @param species The name of the column in the data.frame with the names of species
#' #' @param lat The name of the column in the data.frame with latitudes
#' #' @param lon The name of the column in the data.frame with longitudes
#' RemoveOutliers <- function(points, species="species", lon="decimalLongitude", lat="decimalLatitude") {
#'   npoints_start <- nrow(points)
#'   tmp_points = points
#'   colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
#'   colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
#'   colnames(tmp_points)[colnames(tmp_points)==species] <- "species"
#'   spp <- unique(tmp_points$species)
#'   all_points <- list()
#'   for(species_index in 1:length(spp)){
#'     sp0 <- tmp_points[tmp_points$species==spp[species_index],]
#'     out_lat <- grDevices::boxplot.stats(sp0$y)$out
#'     out_lon <- grDevices::boxplot.stats(sp0$x)$out
#'     sp0 <- sp0[!sp0$y %in% out_lat, ]
#'     sp0 <- sp0[!sp0$x %in% out_lon, ]
#'     all_points[[species_index]] <- sp0
#'   }
#'   points <- do.call(rbind, all_points)
#'   colnames(points)[colnames(points)=="x"] <- lon
#'   colnames(points)[colnames(points)=="y"] <- lat
#'   npoints_end <- nrow(points)
#'   print(paste0(npoints_start - npoints_end, " points removed."))
#'   return(points)
#' }
#'
#' #' Removes points that are located in the wrong country according to their GBIF labels
#' #' That will remove points that are not located in the countries where their labels say they were collected
#' #' @param points A data.frame of distribution points with at least five columns where one column represents species names and other two decimal coordinates.
#' #' @param lat The name of the column in the data.frame with latitudes
#' #' @param lon The name of the column in the data.frame with longitudes
#' #' @param buffer A number of degrees around each country where points are still considered part of that country
#' #' @details The input data.frame must have a column named countryCode and one named gbifID, as the .csv files downloaded directly from GBIF.
#' RemoveWrongCountries <- function(points, lon="decimalLongitude", lat="decimalLatitude", buffer=5, wrld_simpl="") {
#'   npoints_start <- nrow(points)
#'   tmp_points = points
#'   colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
#'   colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
#'   countries <- as.character(wrld_simpl[2]$ISO2)
#'   dubiousGBIF_ids <- c()
#'   for(country_index in 1:length(countries)) {
#'     tmp_country <- countries[country_index]
#'     if(length(which(tmp_points$countryCode %in% tmp_country)) > 0) {
#'       tmp_subset <- tmp_points[tmp_points$countryCode==tmp_country,]
#'       coords <- stats::na.omit(tmp_subset[,c("x","y")])
#'       sp::coordinates(coords) <- ~ x + y
#'       sp::proj4string(coords) <- sp::proj4string(wrld_simpl) <- "+proj=longlat +ellps=WGS84 +no_defs"
#'       country_plus_buffer <- raster::buffer(wrld_simpl[country_index,], buffer) # adding buffer around country
#'       sp::proj4string(country_plus_buffer) <- "+proj=longlat +ellps=WGS84 +no_defs"
#'       answer <- which(is.na(sp::over(coords, country_plus_buffer)))
#'       dubiousGBIF_ids <- c(dubiousGBIF_ids, tmp_subset$gbifID[answer])
#'     }
#'   }
#'   if(!is.null(dubiousGBIF_ids)) {
#'     points_cleaned <- points[-which(points$gbifID %in% dubiousGBIF_ids),]
#'   }
#'   npoints_end <- nrow(points_cleaned)
#'   print(paste0(npoints_start - npoints_end, " points removed."))
#'   return(points_cleaned)
#' }
#'
#' #' Removes points that are located in country centroids
#' #' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' #' @param lat The name of the column in the data.frame with latitudes
#' #' @param lon The name of the column in the data.frame with longitudes
#' #' @param buffer A number in meters around each country centroid for points to be removed
#' RemoveCentroids <- function(points, lon="decimalLongitude", lat="decimalLatitude", buffer=75000) {
#'   npoints_start <- nrow(points)
#'   tmp_points = points
#'   colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
#'   colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
#'   wrld_map <- rworldmap::getMap(resolution="low")
#'   # here the buffer is in meters
#'   # Note: probably should do something about small countries where the buffer of 75km may be too broad
#'   coords <- tmp_points[,c("x","y")]
#'   sp::coordinates(coords) <- ~ x + y
#'   sp::proj4string(coords) <- sp::proj4string(wrld_map)
#'   centroids <- rgeos::gCentroid(wrld_map, byid=TRUE)
#'   centroids_plus_buffer <- raster::buffer(centroids, buffer) # adding buffer around polygons
#'   answer <- which(!is.na(sp::over(coords, centroids_plus_buffer)))
#'   if(length(answer) > 0) {
#'     points <- points[-answer,]
#'     npoints_end <- nrow(points)
#'     print(paste0(npoints_start - npoints_end, " points removed."))
#'     return(points)
#'   } else {
#'     print("no points removed")
#'     return(points) }
#' }
#'
#'
#' #' Removes duplicated latitudes and longitudes for the same species
#' #' @param points A data.frame of distribution points
#' #' @param species The name of the column in the data.frame with the names of species
#' #' @param lat The name of the column in the data.frame with latitudes
#' #' @param lon The name of the column in the data.frame with longitudes
#' RemoveDuplicates <- function(points, species="species", lon="decimalLongitude", lat="decimalLatitude") {
#'   npoints_start <- nrow(points)
#'   tmp_points = points
#'   colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
#'   colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
#'   colnames(tmp_points)[colnames(tmp_points)==species] <- "species"
#'   spp <- unique(tmp_points$species)
#'   all_points <- list()
#'   for(species_index in 1:length(spp)){
#'     tmp_subset <- as.data.frame(tmp_points[tmp_points$species==spp[species_index],])
#'     all_points[[species_index]] <- tmp_subset[-which(duplicated(tmp_subset[,c("x","y")])),]
#'   }
#'   points <- do.call(rbind, all_points)
#'   colnames(points)[colnames(points)=="x"] <- lon
#'   colnames(points)[colnames(points)=="y"] <- lat
#'   npoints_end <- nrow(points)
#'   print(paste0(npoints_start - npoints_end, " points removed."))
#'   return(points)
#' }
#'
#'
#' #' Removes points with coordinates without decimal cases (probably inaccurate)
#' #' @param points A data.frame of distribution points
#' #' @param lat The name of the column in the data.frame with latitudes
#' #' @param lon The name of the column in the data.frame with longitudes
#' RemoveNoDecimal <- function(points, lon="decimalLongitude", lat="decimalLatitude") {
#'   npoints_start <- nrow(points)
#'   tmp_points = points
#'   colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
#'   colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
#'   length_decimal_lat <- nchar(sub("^[^.]*", "", tmp_points$y))
#'   length_decimal_lon <- nchar(sub("^[^.]*", "", tmp_points$x))
#'   points <- points[which(length_decimal_lat>=1 & length_decimal_lon>=1),]
#'   npoints_end <- nrow(points)
#'   print(paste0(npoints_start - npoints_end, " points removed."))
#'   return(points)
#' }
#'
#'
#' #' Removes points around large herbaria.
#' #' @param points A data.frame of distribution points with at least three columns where one column represents species names and other two decimal coordinates.
#' #' @param lat The name of the column in the data.frame with latitudes
#' #' @param lon The name of the column in the data.frame with longitudes
#' RemoveHerbariaLocalities <- function(points, lon="decimalLongitude", lat="decimalLatitude") {
#'   npoints_start <- nrow(points)
#'   tmp_points = points
#'   colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
#'   colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
#'   herbaria_localities <- readRDS("R/allHerbaria_ADM1_badCoords.Rdata") # Jeremy provided this list
#'   # herbaria_localities <- read.csv("~/Desktop/MiSSEgradient/MiSSEGradient/Sampling/data/allHerbaria_ADM1_badCoords.txt", header=FALSE)
#'   for(herb.index in 1:length(herbaria_localities)){
#'     points <- tmp_points[which(!round(tmp_points$y,2)==herbaria_localities[herb.index,1] | !round(tmp_points$x,2)==herbaria_localities[herb.index,2]),]
#'   }
#'   colnames(points)[colnames(points)=="x"] <- lon
#'   colnames(points)[colnames(points)=="y"] <- lat
#'   npoints_end <- nrow(points)
#'   print(paste0(npoints_start - npoints_end, " points removed."))
#'   return(points)
#' }
#'
#'
#' #' Flag species that have "weird" distributions in the dataset
#' #' @param points A data.frame of distribution points with at least three columns where one column represents species names and other two decimal coordinates.
#' #' @param lat The name of the column in the data.frame with latitudes
#' #' @param lon The name of the column in the data.frame with longitudes
#' UnusualDistributions <- function (points, species="species", lat="decimalLatitude", lon="decimalLongitude", buffer.polygon=c(5,10)) {
#'   tmp_points = as.data.frame(points)
#'   tmp_points = tmp_points[,c(which(colnames(tmp_points)==species),which(colnames(tmp_points)==lon),which(colnames(tmp_points)==lat))]
#'   colnames(tmp_points) <- c("species","lon","lat")
#'   spp <- unique(tmp_points$species)
#'
#'   # WWFload is taken from speciesgeocodeR; all credit goes to the original authors
#'   WWFload <- function(x = NULL) {
#'     if (missing(x)) {
#'       x <- getwd()
#'     }
#'     download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
#'                   destfile = file.path(x, "wwf_ecoregions.zip"), quiet=TRUE)
#'     unzip(file.path(x, "wwf_ecoregions.zip"), exdir = file.path(x, "WWF_ecoregions"))
#'     file.remove(file.path(x, "wwf_ecoregions.zip"))
#'     wwf <- maptools::readShapeSpatial(file.path(x, "WWF_ecoregions", "official",
#'                                                 "wwf_terr_ecos.shp"))
#'     return(wwf)
#'   }
#'   wwf <- WWFload(tempdir())
#'   result <- matrix(nrow=length(spp), ncol=2)
#'   for(species_index in 1:length(spp)) {
#'     one_species <- tmp_points[tmp_points$species==spp[species_index],]
#'     locations.spatial <- sp::SpatialPointsDataFrame(coords=one_species[,c("lon", "lat")], data=one_species)
#'     mappedregions <- sp::over(locations.spatial, wwf)
#'
#'     print(species_index)
#'     if(length(which(table(mappedregions$REALM)!=0)) > 2) {
#'       result[species_index,1] <- spp[species_index]
#'       result[species_index,2] <- "flagged"
#'     } else {
#'       result[species_index,1] <- spp[species_index]
#'       result[species_index,2] <- "ok"
#'
#'     }
#'   }
#'   flagged <- result[,1][which(result[,2]=="flagged")]
#'   answer <- c()
#'   for(flagged_index in 1:length(flagged)) {
#'     #print(flagged_index)
#'     subset <- points[points$species==flagged[flagged_index],]
#'     subset1 <- subset[,c("decimalLongitude","decimalLatitude")]
#'     distances <- as.matrix(dist(subset1))
#'     mean_distances <- apply( distances, 1, function(x) mean( x[order(x)][2:4] ) )
#'     outliers <- subset[which(mean_distances > mean(distances)),]
#'     if(nrow(outliers)>0){
#'       subset <- subset[-which(mean_distances > mean(distances)),]
#'     }
#'     max.lat <- ceiling(max(subset[,"decimalLatitude"])) + buffer.polygon[1]
#'     min.lat <- floor(min(subset[,"decimalLatitude"])) - buffer.polygon[1]
#'     max.lon <- ceiling(max(subset[,"decimalLongitude"])) + buffer.polygon[2]
#'     min.lon <- floor(min(subset[,"decimalLongitude"])) - buffer.polygon[2]
#'     if(any(c(max.lat > min.lat + 75, max.lon > min.lon + 125))){
#'       answer <- c(answer,flagged[flagged_index])
#'     }
#'   }
#'   return(answer)
#' }

#' Thinning distribution data to smooth sampling bias
#' @param points A data.frame of three columns containing species and coordinates
#' @param species A character string indicating name of column with species names
#' @param lat A character string indicating name of column with latitudes
#' @param lon character string indicating name of column with longitudes
#' @param n A number indicating how many points to keep in each cell after thinning
Thinning <- function(
  points,
  species = "scientificName",
  lat = "decimalLatitude",
  lon = "decimalLongitude",
  n = 1
) {
  tmp_points <- points
  colnames(tmp_points)[colnames(tmp_points) == lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points) == lat] <- "y"
  colnames(tmp_points)[colnames(tmp_points) == species] <- "tmp_sp"
  spp <- unique(tmp_points[, tmp_sp])
  results <- list()
  for (species_index in 1:length(spp)) {
    coords <- tmp_points[
      tmp_points[, tmp_sp] == spp[species_index],
      c("y", "x")
    ]
    coords <- coords[!duplicated(coords[, "x"]) & !duplicated(coords[, "y"]), ]
    if (nrow(coords) > 1) {
      sp::coordinates(coords) <- ~ y + x
      raster::crs(
        coords
      ) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      r0 <- raster::raster(coords)
      raster::res(r0) <- 1 # cell resolution
      r0 <- raster::extend(r0, raster::extent(r0) + 5)
      res <- cbind(
        spp[species_index],
        as.data.frame(dismo::gridSample(coords, r0, n))
      ) # n = maximum number of points per cell
      colnames(res) <- c("tmp_sp", "lat", "lon")
      results[[species_index]] <- res
    } else {
      res <- cbind(spp[species_index], coords)
      colnames(res) <- c("tmp_sp", "lat", "lon")
      results[[species_index]] <- res
    }
  }
  results <- do.call(rbind, results)
  colnames(results) <- c(species, lat, lon)
  return(results)
}


#' Get values of selected climatic variables from filtered points
#' @param points A data.frame of three columns containing species and coordinates
#' @param species A character string indicating name of column with species names
#' @param lat A character string indicating name of column with latitudes
#' @param lon character string indicating name of column with longitudes
#' @param layerdir name of directory with climate layers
ClimateFromPoint_custom <- function(
  points,
  species = "species",
  lon = "lon",
  lat = "lat",
  layerdir = ""
) {
  tmp_points <- points
  colnames(tmp_points)[which(colnames(tmp_points) == lon)] <- "lon"
  colnames(tmp_points)[which(colnames(tmp_points) == lat)] <- "lat"
  colnames(tmp_points)[which(colnames(tmp_points) == species)] <- "species"
  tmp_points <- tmp_points[, c("species", "lat", "lon")]
  # Load climatic layers
  temp <- raster::raster(paste0(layerdir, "/current_30sec/bio_1.tif"))
  prec <- raster::raster(paste0(layerdir, "/current_30sec/bio_12.tif"))
  pet <- raster::raster(paste0(layerdir, "/et0_yr/et0_yr.tif"))
  aridity <- raster::raster(paste0(layerdir, "/ai_et0/ai_et0.tif"))
  bio <- list(temp, prec, pet, aridity)
  vars <- c(temp, prec, pet, aridity)
  names(vars) <- c("temp", "prec", "pet", "aridity")
  final_matrix <- matrix(nrow = nrow(tmp_points), ncol = length(vars))
  cat("Extracting climatic information of", nrow(tmp_points), "points", "\n")
  sp::coordinates(tmp_points) <- ~ lon + lat
  for (var_index in 1:length(vars)) {
    layer <- bio[[var_index]]
    cat("\r", names(vars)[var_index])
    cat("", "\n")
    values <- raster::extract(layer, tmp_points)
    final_matrix[, var_index] <- values
  }
  colnames(final_matrix) <- names(vars)
  result <- cbind(tmp_points, final_matrix)
  return(as.data.frame(result))
}

GetClimateSummStats_custom <- function(points, type = c("raw", "transformed")) {
  tmp_points <- points[, -which(colnames(points) %in% c("lon", "lat"))]
  spp <- unique(tmp_points$species)
  vars <- colnames(tmp_points)[2]
  n_i <- c()
  sigma2_wi <- c()
  summ_stats <- matrix(nrow = length(spp), ncol = 5)
  for (species_index in 1:length(spp)) {
    sp1 <- tmp_points[tmp_points$species == spp[species_index], ]
    cat("Now doing species", species_index, "\r")
    values <- sp1[, 2]
    values <- values[!is.na(values)]
    if (type == "raw") {
      if (vars %in% c("bio_1", "bio_4", "bio_5", "bio_6")) {
        values <- (values / 10)
      }
      n_i[species_index] <- length(values) # sample size
      sigma2_wi[species_index] <- ifelse(is.na(var(values)), 0, var(values)) # sample variance
    }
    if (type == "transformed") {
      if (vars %in% c("bio_1", "bio_4", "bio_5", "bio_6")) {
        values <- (values / 10) + 273.15 # transforms to Kelvin
      }
      values <- log(values) # log
      n_i[species_index] <- length(values) # sample size
      sigma2_wi[species_index] <- ifelse(is.na(var(values)), 0, var(values)) # sample variance

      if (any(values == -Inf)) {
        values <- values[-which(values == -Inf)]
      }
    }
    n0 <- length(values)
    mean0 <- round(mean(values), 6)
    sd0 <- round(stats::sd(values), 6)
    se0 <- round(sd0 / sqrt(n0), 6)
    tmp_summ_stats <- c(n0, mean0, sd0, se0)
    summ_stats[species_index, ] <- c(spp[species_index], tmp_summ_stats)
  }

  sigma2_w <- sum(sigma2_wi * (n_i - 1)) / sum(n_i - 1)
  within_sp_var <- round(sigma2_w / n_i, 6)
  summ_stats <- cbind(summ_stats, within_sp_var)
  colnames(summ_stats)[6] <- paste0("within_sp_var_", vars)
  colnames(summ_stats) <- c(
    "species",
    paste0("n_", vars),
    paste0("mean_", vars),
    paste0("sd_", vars),
    paste0("se_", vars),
    paste0("within_sp_var_", vars)
  )
  return(as.data.frame(summ_stats))
}

#'
#' #' Estimates range and range size based on species distribution modeling.
#' #' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' #' @param species The name of the column in the data.frame with the names of species
#' #' @param lat The name of the column in the data.frame with latitudes
#' #' @param lon The name of the column in the data.frame with longitudes
#' #' @param threshold A number between 0 and 1 for the threshold to make models binary
#' #' @param buffer The radius in kilometers around points to estimate range size when there are three or less valid points
#' #' @param res The resolution (2.5, 5 or 10) of climatic variables used for modeling
#' #' @return A list containing information about the species distribution modeling performance and a raster with the possible range of each species
#' #' @export
#' #'
#' GetRanges <- function(points, species="species", lat="decimalLatitude", lon="decimalLongitude", threshold=0.75, buffer=25, res=10) {
#'   tmp_points = as.data.frame(points)
#'   tmp_points = tmp_points[,c(which(colnames(tmp_points)==species),which(colnames(tmp_points)==lon),which(colnames(tmp_points)==lat))]
#'   colnames(tmp_points) <- c("species","lon","lat")
#'   spp <- unique(tmp_points[,1])
#'   #list_of_ranges <- list()
#'   predictors <- LoadWcLayers(res.layers=res)
#'   for(species_index in 1:length(spp)) {
#'     points_for_range <- tmp_points[tmp_points$species==spp[species_index],]
#'     tmp_list <- try(GetOneRange(points_for_range, threshold, buffer, res, predictors))
#'     if(exists("tmp_list")) {
#'       saveRDS(tmp_list, file=paste0("shapefile_bees/",spp[species_index],".Rdata"))
#'       #list_of_ranges[[species_index]] <- tmp_list
#'       #names(list_of_ranges)[species_index] <- spp[species_index]
#'       rm(tmp_list)
#'     } else {
#'       next }
#'     #cat("","\n")
#'     #cat(species[species_index], "done.")
#'     #cat("","\n")
#'     cat("\r", species_index, " out of ",length(spp))
#'   }
#'   #try(unlink("wc2-5", recursive = TRUE))
#'   #try(unlink("wc5", recursive = TRUE))
#'   #try(unlink("wc10", recursive = TRUE))
#'   #return(list_of_ranges)
#' }
#'
#' GetRanges_angio <- function(points, species="species", lat="decimalLatitude", lon="decimalLongitude", threshold=0.75, buffer=25, res=10) {
#'   tmp_points = as.data.frame(points)
#'   tmp_points = tmp_points[,c(which(colnames(tmp_points)==species),which(colnames(tmp_points)==lon),which(colnames(tmp_points)==lat))]
#'   colnames(tmp_points) <- c("species","lon","lat")
#'   spp <- unique(tmp_points[,1])
#'   #list_of_ranges <- list()
#'   predictors <- LoadWcLayers(res.layers=res)
#'   for(species_index in 1:length(spp)) {
#'     points_for_range <- tmp_points[tmp_points$species==spp[species_index],]
#'     tmp_list <- try(GetOneRange(points_for_range, threshold, buffer, res, predictors))
#'     if(exists("tmp_list")) {
#'       saveRDS(tmp_list, file=paste0("shapefile_angio_scndtry/",spp[species_index],".Rdata"))
#'       #list_of_ranges[[species_index]] <- tmp_list
#'       #names(list_of_ranges)[species_index] <- spp[species_index]
#'       rm(tmp_list)
#'     } else {
#'       next }
#'     #cat("","\n")
#'     #cat(species[species_index], "done.")
#'     #cat("","\n")
#'     cat("\r", species_index, " out of ",length(spp))
#'   }
#'   #try(unlink("wc2-5", recursive = TRUE))
#'   #try(unlink("wc5", recursive = TRUE))
#'   #try(unlink("wc10", recursive = TRUE))
#'   #return(list_of_ranges)
#' }
#'
#' Bg_cHull <- function (thinned_points, width=2.5) {
#'   geo <- thinned_points
#'   sp::coordinates(geo) <- ~ lon + lat
#'   x <- rgeos::gConvexHull(geo)
#'   buffer <- rgeos::gBuffer(x, width, byid=F)
#'   return(buffer)
#' }
#'
#' #make.new.dir <- function(dir, namedir, overwrite.dir) {
#' #  new.dir <- paste0(getwd(),"/", namedir)
#' #  if (overwrite.dir) {
#' #    unlink(new.dir, recursive=TRUE) # overwrites output directory
#' #    dir.create(file.path(new.dir))
#' #  } else { suppressWarnings(dir.create(new.dir)) }
#' #  return(new.dir)
#' #}
#'
#' #' @param points_for_range points_for_range
#' #' @param threshold threshold
#' #' @param buffer buffer
#' #' @param res res
#' GetOneRange <- function(points_for_range, threshold, buffer, res, predictors) {
#'   cat("Thinning points...")
#'   if(nrow(points_for_range) > 2) {
#'     thinned_points <- Thinning_internal(points_for_range)
#'   } else { thinned_points = points_for_range}
#'   #thinned_points = points_for_range
#'   #cat("Loading environmental predictors...")
#'   cat("Creating background polygon...")
#'   bg <- BgPolygon(thinned_points)
#'   bg_chull <- Bg_cHull(thinned_points, width=2.5)
#'   if(nrow(thinned_points) < 3) { # If two or fewer valid points, the range will be retrieved from a circle around these points
#'     list_of_model_results <- RangeFromFewPoints(thinned_points, predictors, buffer)
#'   } else {
#'     cat("Removing predictors with colinearity problems ...")
#'     predictors_final <- ColinearityTest(bg, predictors)
#'     cat("Performing species distribution modeling...")
#'     list_of_model_results <- RangeFromSDM_dismo(thinned_points, predictors_final, bg)
#'   }
#'   cat("Making models binary based on threshold...")
#'   fullresults <- RangeSize(list_of_model_results, threshold, bg_chull)
#'   cat("Adding alerts...")
#'   fullresults_w_alerts <- AddAlerts(fullresults, bg)
#'   return(fullresults_w_alerts)
#' }
#'
#'
#' #' Performs sdm using dismo's maxent
#' #' For great maxent and sdm tutorials see:
#' #' jcoliver.github.io/learn-r/011-species-distribution-models.html
#' #' cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf
#' #' etc
#' #' @importFrom dismo randomPoints kfold maxent predict evaluate
#' #' @importFrom raster raster
#' #' @param thinned_points thinned_points
#' #' @param predictors_final predictors_final
#' #' @param bg bg
#' RangeFromSDM_dismo <- function (thinned_points, predictors_final, bg) {
#'   points_for_sdm <- thinned_points[,c("lon","lat")]
#'   species_name <- as.character(thinned_points[1,1])
#'   mask <- raster::raster(predictors_final[[1]])
#'   background <- dismo::randomPoints(mask = mask, n = abs((bg[1] - bg[2]) * (bg[3] - bg[4])), # Number of random background points
#'                                     ext = bg, extf = 1.25)
#'   background <- as.data.frame(background)
#'   testing_group <- 1 # arbitrarily assign group 1 as the testing data group
#'   group_presence <- dismo::kfold(points_for_sdm, k = 2) # randomly divide dataset into two groups
#'   presence_train <- points_for_sdm[group_presence != testing_group, ]
#'   presence_test <- points_for_sdm[group_presence == testing_group, ]
#'
#'   group_background <- dismo::kfold(background, k = 2) # background points
#'   background_train <- background[group_background != testing_group, ]
#'   background_test <- background[group_background == testing_group, ]
#'
#'   sdm_raw <- dismo::maxent(x = predictors_final, p = presence_train, nbg=1000)
#'   predict <- dismo::predict(object = sdm_raw, x = predictors_final, ext = bg)
#'
#'   eval <- dismo::evaluate(p = presence_test, a = background_test, model = sdm_raw, x = predictors_final)
#'   auc0 <- methods::slot(eval, "auc")
#'
#'   results <- list()
#'   results[[1]] <- paste0(species_name)
#'   results[[2]] <- predict
#'   results[[3]] <- paste("maxent (dismo)")
#'   results[[4]] <- names(predictors_final)
#'   results[[5]] <- round(auc0, 2)
#'   names(results) <- c("species_name","original_model", "sdm_method", "predictors","auc")
#'
#'   return(results)
#' }
#'
#' #' Detecting collinearity in predictors
#' #' @importFrom raster crop stack
#' #' @importFrom usdm exclude vifcor
#' #' @param bg bg
#' #' @param predictors predictors
#' ColinearityTest <- function(bg, predictors) {
#'   predictors_final <- raster::crop(predictors, raster::extent(bg))
#'   predictors_final <- raster::stack(predictors_final)
#'   layers1 <- as.data.frame(predictors_final)
#'   v0 <- suppressWarnings(usdm::vifcor(layers1, th=0.8)) # established threshold
#'   layers1 <- usdm::exclude(layers1, v0) # excludes variables that have collinearity problems
#'   predictors_final <- predictors_final[[which(names(predictors_final)%in%colnames(layers1))]]
#'   return(predictors_final)
#' }
#'
#' #' Internal function to add alerts to the results
#' #' @param fullresults fullresults
#' #' @param bg bg
#' AddAlerts <- function(fullresults, bg) {
#'   if(length(fullresults)==5){
#'     # is the number of good points three or less? (data deficiency alert)
#'     fullresults[[6]] <- "Data deficiency alert: species with three or less valid points."
#'     names(fullresults)[6] <- "alerts"
#'     return(fullresults)
#'   } else {
#'     alerts <- c()
#'     #cross check with list of crops and invasive species
#'     #crops <- readRDS("R/crops_taxized_gbif.Rdata")
#'     #invasive <- readRDS("R/invasive_taxized_gbif.Rdata")
#'     #if(fullresults$species_name %in% crops) {
#'     #crop_alert <- "Crop alert: this species is listed as a crop at either fao.org, hort.purdue.edu/newcrop/ or Meyer et al. (2012)."
#'     # Meyer, R. S., DuVal, A. E. & Jensen, H. R. Patterns and processes in crop domestication: an historical review and quantitative analysis of 203 global food crops. New Phytol. 196, 29–48 (2012).
#'     # alerts <- c(alerts, crop_alert)
#'     #}
#'     #if(fullresults$species_name %in% invasive) {
#'     #  invasive_alert <- "Invasive alert: this species is listed as invasive at GISD. Check www.iucngisd.org for more information."
#'     #  alerts <- c(alerts, invasive_alert)
#'     #}
#'     # is AUC very low?
#'     if(fullresults$auc < 0.75) {
#'       low_auc_alert <- "Alert of low AUC: model with AUC lower than 0.75."
#'       alerts <- c(alerts, low_auc_alert)
#'     }
#'     # is range too wide? (proxy for non-natural distribution)
#'     lims <- bg[]
#'     if(lims[2] > lims[1] + 100 & lims[4] > lims[3] + 50) {
#'       wide_range_alert <- "Alert of wide range: range seems too wide for a natural distribution. Check quality of distribution data."
#'       alerts <- c(alerts, wide_range_alert)
#'     } else if(fullresults$range_size > 10000000) {
#'       wide_range_alert <- "Alert of wide range: range seems too wide for a natural distribution. Check quality of distribution data."
#'       alerts <- c(alerts, wide_range_alert)
#'     }
#'     if(is.null(alerts)) {
#'       fullresults[[9]] <- "No alerts returned."
#'       names(fullresults)[9] <- "alerts"
#'     } else {
#'       fullresults[[9]] <- alerts
#'       names(fullresults)[9] <- "alerts"
#'     }
#'     return(fullresults)
#'   }
#' }
#'
#' #' Performs spatial thinning before modeling by selecting one occurence per grid cell
#' #' @importFrom dismo gridSample
#' #' @importFrom raster crs raster res extend extent
#' #' @importFrom sp coordinates
#' #' @param points_for_range points_for_range
#' Thinning_internal <- function(points_for_range) {
#'   coords <- points_for_range[,c(3:2)]
#'   sp::coordinates(coords) <- ~ lat + lon
#'   raster::crs(coords) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
#'   coords@bbox[1:4] <- c(coords@bbox[1] - 1, coords@bbox[2] - 1, coords@bbox[3] + 1, coords@bbox[4] + 1) # that's to avoid errors when min and max lon and lat are the same
#'   r0 <- raster::raster(coords)
#'   raster::res(r0) <- 1 # cell resolution
#'   r0 <- raster::extend(r0, raster::extent(r0) + 5) # expand the extent of the RasterLayer a little
#'   thinned_points <- as.data.frame(dismo::gridSample(coords, r0, n = 1)) # n = maximum number of points per cell
#'   species <- points_for_range$species[sequence(nrow(thinned_points))]
#'   thinned_points <- cbind(species, thinned_points)
#'   return(thinned_points)
#' }
#'
#' #' Getting a proxy for range size when you have three points or less.
#' #' When there are three or less valid points, the distribution will be estimated from a circle around each point of radius equal to the argument buffer.
#' #' @importFrom sp coordinates polygons
#' #' @importFrom raster raster crop mask
#' #' @importFrom dismo circles
#' #' @param thinned_points thinned_points
#' #' @param predictors predictors
#' #' @param buffer buffer
#' RangeFromFewPoints <- function(thinned_points, predictors, buffer) {
#'   point <- thinned_points
#'   species_name <- as.character(point[1,1])
#'   point[,1] <- 1
#'   sp::coordinates(point) <- ~ lon + lat
#'   area_mask <- raster::raster(res=raster::res(predictors))
#'   area_mask[] <- 1
#'   buffer_mask = buffer # Radius of the circle in kilometers
#'   circle_around_point <- dismo::circles(point, d=buffer_mask*1000, lonlat=TRUE)
#'   circle_around_point <- sp::polygons(circle_around_point)
#'   area_mask <- raster::crop(area_mask, circle_around_point)
#'   area_mask <- raster::mask(area_mask, circle_around_point)
#'   results <- list()
#'   results[[1]] <- paste0(species_name)
#'   results[[2]] <- area_mask
#'   results[[3]] <- "no AUC"
#'   names(results) <- c("species_name","range","auc")
#'   return(results)
#' }
#'
#' #' Estimates range size from binary model
#' #' @importFrom raster area
#' #' @param list_of_model_results list_of_model_results
#' #' @param threshold threshold
#' RangeSize <- function (list_of_model_results, threshold, bg_chull) {
#'   if(length(list_of_model_results) == 3) {
#'     range_size <- round(sum(raster::area(list_of_model_results$range, na.rm=T)[], na.rm=T), 2)
#'     list_of_model_results[[4]] <- range_size
#'     list_of_model_results[[5]] <- paste("Range estimated from three points or less.")
#'     names(list_of_model_results)[c(4,5)] <- c("range_size", "note")
#'     return(list_of_model_results)
#'   } else {
#'     model <- list_of_model_results$original_model
#'     model[model[] < threshold] <- NA
#'     model[!is.na(model)] <- 1
#'     model <- raster::mask(model, bg_chull)
#'     range_size <- round(sum(raster::area(model, na.rm=T)[], na.rm=T), 2)
#'     list_of_model_results[[6]] <- model
#'     list_of_model_results[[7]] <- paste(threshold)
#'     list_of_model_results[[8]] <- range_size
#'     names(list_of_model_results)[c(6,7,8)] <- c("range", "threshold","range_size")
#'     return(list_of_model_results)
#'   }
#' }
#'
#' #' Internal function -- gets polygon around distribution to crop predictors before modeling
#' #' @importFrom raster extent
#' #' @param thinned_points thinned_points
#' #' @param buffer.polygon buffer.polygon
#' BgPolygon <- function (thinned_points, buffer.polygon=c(5, 10)) {
#'   max.lat <- ceiling(max(thinned_points[,"lat"])) + buffer.polygon[1]
#'   min.lat <- floor(min( thinned_points[,"lat"])) - buffer.polygon[1]
#'   max.lon <- ceiling(max(thinned_points[,"lon"])) + buffer.polygon[2]
#'   min.lon <- floor(min(thinned_points[,"lon"])) - buffer.polygon[2]
#'   #if(max.lat > min.lat + 50 & max.lon > min.lon + 100){
#'   #  warning("Very wide background detected - this may not correspond to a natural distribution.")
#'   #}
#'   bg.area <- raster::extent(x = c(min.lon, max.lon, min.lat, max.lat))
#'   if(bg.area[1] > 180 | bg.area[2] > 180 |  bg.area[2] < -180 | bg.area[1] < -180){
#'     bg.area[which(bg.area[] < -180)] <- -179.9
#'     bg.area[which(bg.area[] > 180)] <- 179.9
#'   }
#'   #if(plot.polygons) {
#'   #  polygon.dir <- paste0(getwd(),"/polygon")
#'   #  if (overwrite.dir) {
#'   #    unlink(polygon.dir, recursive=TRUE) # overwrites output directory for models
#'   #    dir.create(file.path(polygon.dir))
#'   #  } else { suppressWarnings(dir.create(polygon.dir)) }
#'   #  pdf(file=paste0(polygon.dir,"/", thinned_points[1,1],"__background_polygon.pdf"), height=6, width=10)
#'   #  plot(wrld_simpl)
#'   #  plot(bg.area, col="orange", add=T)
#'   #  graphics::title(main=paste0(thinned_points[1,1]))
#'   #  dev.off()
#'   #}
#'   return(bg.area)
#' }

#' Internal function for now -- loads climatic layers from Worldclim
#' @importFrom raster getData
#' @param res.layers res.layers
LoadWcLayers <- function(res.layers) {
  bio <- raster::getData(
    "worldclim",
    var = "bio",
    res = res.layers,
    path = getwd()
  ) # climatic layers
  return(bio)
}

#' A function to overlay distributions and create a map of species richness
#' @param list_of_ranges A list in the format of the output of GetRanges()
#' @return A raster of species richness
#' @importFrom raster resample calc stack mask
#' @export

GetSpRichness <- function(list_of_ranges) {
  template.map <- readRDS("original_data/template.map.Rdata")
  res(template.map)[] <- res(list_of_ranges[[1]])
  if (length(list_of_ranges) > 1) {
    r0 <- list_of_ranges[[1]]
    r0 <- raster::resample(r0, template.map)
    cat(1, "\r")
    sum_raster <- r0
    for (j in 2:length(list_of_ranges)) {
      r1 <- list_of_ranges[[j]]
      r1 <- raster::resample(r1, template.map)
      sum_raster <- calc(stack(list(sum_raster, r1)), sum, na.rm = T)
      cat(j, "\r")
    }
  } else {
    r0 <- list_of_ranges[[1]]
    r0 <- raster::resample(r0, template.map)
    sum_raster <- calc(stack(list(r0, template.map)), sum, na.rm = T)
    cat("Only one valid range in list. Returning original list.")
  }
  return(sum_raster)
}

#' A function to map distribution of traits (still with no PW)
#' @param list_of_ranges A list in the format of the output of GetRanges()
#' @param trait_data A data.frame where the first column contain species names and the second contains trait data
#' @param type "binary" or "continuous", depending on the kind of trait data
#' @return A raster with the traits mapped in space
#' @importFrom raster resample calc stack mask
#' @export
GetTraitDistribution <- function(
  list_of_ranges,
  trait_data,
  type = c("binary", "continuous")
) {
  if (length(list_of_ranges) != nrow(trait_data)) {
    stop("Length of list_of_ranges and trait_data do not match.")
  }
  if (any(!names(list_of_ranges) %in% trait_data[, 1])) {
    stop("Some species in the list_of_ranges are not in the trait_data")
  }
  ranges <- lapply(list_of_ranges, function(x) x$range)
  xmin <- min(unlist(lapply(lapply(ranges, "extent"), "xmin")))
  xmax <- max(unlist(lapply(lapply(ranges, "extent"), "xmax")))
  ymin <- min(unlist(lapply(lapply(ranges, "extent"), "ymin")))
  ymax <- max(unlist(lapply(lapply(ranges, "extent"), "ymax")))
  template.map <- readRDS("R/template.map.Rdata")
  #template.map <- raster::getData("worldclim", var="bio", download=TRUE, res=10)[[1]]
  #template.map[!is.na(template.map)] <- 0
  if (type == "continuous") {
    tmp.raster.list_traits <- list()
    tmp.raster.list_sprich <- list()
    for (range_index in 1:length(ranges)) {
      cat("\r", range_index, " out of ", length(ranges))
      species <- names(ranges)[[range_index]]
      r1 <- ranges[[range_index]]
      r1 <- raster::resample(r1, template.map)
      r1[is.na(r1)] <- 0
      tmp.raster.list_sprich[[range_index]] <- raster::mask(r1, template.map)
      r1[r1 == 1] <- trait_data[trait_data$species == species, 2]
      tmp.raster.list_traits[[range_index]] <- raster::mask(r1, template.map)
    }
    traits_sum <- raster::calc(raster::stack(tmp.raster.list_traits), sum)
    sp_sum <- raster::calc(raster::stack(tmp.raster.list_sprich), sum)
    result <- r1
    result[!is.na(result)] <- traits_sum / sp_sum
    result <- mask(result, template.map)
    result <- crop(result, extent(xmin, xmax, ymin, ymax) + +10)
    return(result)
  } else if (type == "binary") {
    trait_states_list <- list()
    states <- unique(trait_data[, 2])
    if (length(states) != 2) {
      stop("Trait is not binary")
    }
    #trait_data[trait_data[,2]==states[1],2] <- 0
    #trait_data[trait_data[,2]==states[2],2] <- 1
    for (state_index in 1:length(states)) {
      t0 <- trait_data[which(trait_data[, 2] == states[state_index]), ]
      r0 <- ranges[which(names(ranges) %in% t0[, 1])]
      tmp.raster.list <- list()
      for (range_index in 1:length(r0)) {
        r1 <- r0[[range_index]]
        r1 <- raster::resample(r1, template.map)
        r1[is.na(r1)] <- 0
        tmp.raster.list[[range_index]] <- raster::mask(r1, template.map)
      }
      trait_states_list[[state_index]] <- calc(stack(tmp.raster.list), sum)
      names(trait_states_list)[state_index] <- states[state_index]
    }
    total_sprich <- raster::calc(raster::stack(trait_states_list), sum)
    test <- total_sprich
    proportion_trait1 <- ((trait_states_list[[1]] * 100) / total_sprich) / 100
    proportion_trait2 <- ((trait_states_list[[2]] * 100) / total_sprich) / 100
    final <- c(proportion_trait1, proportion_trait2)
    names(final) <- names(trait_states_list)
    final <- crop(final, extent(xmin, xmax, ymin, ymax) + +10)
    return(final)
  }
}

#remove_data_raster <- function(x) {
#  x[x[]==0] <- NA
#  x[!is.na(x)] <- 0
#  x
#}

#range01 <- function(x){(x-min(x))/(max(x)-min(x))}

#' A function to map distribution of traits
#' @param raster1 A raster of mapped trait distribution or species-richness
#' @param raster2 A raster of mapped trait distribution or species-richness
#' @return A raster with mapped residuals of a linear regression between raster1 and raster2
#' @importFrom raster getValues crop extent
#' @importFrom stats lm na.exclude residuals.lm
#' @export
GetResDistribution <- function(raster1, raster2) {
  #if(is.null(template.map)) {
  template.map <- readRDS("data/template.map.Rdata")
  #template.map <- raster::getData("worldclim", var="bio", download=TRUE, res=10)[[1]]
  #template.map[!is.na(template.map)] <- 0
  #} else { template.map=template.map }
  template <- crop(template.map, raster::extent(raster1))
  # set pallete
  raster1[raster1[] == 0] <- NA
  raster2[raster2[] == 0] <- NA
  l.model <- stats::lm(
    raster::getValues(raster1) ~ raster::getValues(raster2),
    na.action = na.exclude
  )
  res.raster <- template
  res.raster[] <- as.numeric(stats::residuals.lm(l.model))
  return(res.raster)
}

FilterRanges <- function(rangers_output, label_fail = "species_fail") {
  list_of_ranges_tmp <- list()
  for (i in 1:length(rangers_output)) {
    if (length(rangers_output[[i]]) > 1) {
      list_of_ranges_tmp[[i]] <- rangers_output[[i]]$range
      names(list_of_ranges_tmp)[i] <- names(rangers_output)[i]
    }
  }
  species_fail <- names(rangers_output)[which(unlist(lapply(
    list_of_ranges_tmp,
    is.null
  )))]
  write.csv(species_fail, file = paste0(label_fail, ".csv"), row.names = F)
  list_of_ranges_tmp[which(unlist(lapply(list_of_ranges_tmp, is.null)))] <- NULL
  rangers_output <- list_of_ranges_tmp
  return(rangers_output)
}

DataFromPoints <- function(points, layer) {
  if (any(colnames(points) != c("species", "lat", "lon"))) {
    stop(
      "Columns have to be in the order of taxon, latitude and longitude and named as 'species', 'lat', and 'lon"
    )
  }
  if (ncol(points) != 3) {
    stop(
      "Dataset should be of class data.frame and organized in three columns named as 'species', 'lat', and 'lon'"
    )
  }
  if (!is.data.frame(points)) {
    stop(
      "Dataset should be of class data.frame and organized in three columns named as 'species', 'lat', and 'lon'"
    )
  }
  cat("Extracting climatic information of", nrow(points), "points", "\n")
  colnames(points) <- c("species", "lat", "lon")
  sp::coordinates(points) <- ~ lon + lat
  values <- raster::extract(layer, points)
  result <- cbind(points, values)
  out <- as.data.frame(result)
  if (class(layer@data@names) == "character") {
    colnames(out)[2] <- layer@data@names
  }
  return(out)
}

GetClimateSummStats <- function(climate_by_point, convert = TRUE) {
  original_species_order <- unique(climate_by_point[, 1])
  if (convert) {
    summ_stats <- aggregate(
      climate_by_point[, 2],
      by = list(climate_by_point[, 1]),
      function(x) BasicSummStats.k(x)
    )
  } else {
    summ_stats <- aggregate(
      climate_by_point[, 2],
      by = list(climate_by_point[, 1]),
      function(x) BasicSummStats(x)
    )
  }
  summ_stats <- cbind(summ_stats[, 1], as.data.frame(summ_stats[, -1]))
  colnames(summ_stats) <- c("species", "n", "mean", "sd", "se")
  summ_stats <- summ_stats[match(original_species_order, summ_stats[, 1]), ]
  rownames(summ_stats) <- NULL
  return(summ_stats)
}


CalcSpRich <- function(list_of_shapes, template.map) {
  if (length(list_of_shapes) > 1) {
    r0 <- list_of_shapes[[1]]
    r0 <- raster::resample(r0, template.map)
    cat(1, "\r")
    sum_raster <- r0
    for (j in 2:length(list_of_shapes)) {
      r1 <- list_of_shapes[[j]]
      r1 <- raster::resample(r1, template.map)
      sum_raster <- calc(stack(list(sum_raster, r1)), sum, na.rm = T)
      cat(j, "\r")
    }
    mean_raster <- sum_raster
  } else {
    mean_raster <- list_of_shapes
  }
  return(mean_raster)
}

# setting up corhmm model
corhmm.model.setup <- function(dataset) {
  RateCat1 <- getStateMat4Dat(dataset)$rate.mat # R1
  RateCat2 <- getStateMat4Dat(dataset)$rate.mat # R1
  RateClassMat <- getRateCatMat(2) #
  RateClassMat <- equateStateMatPars(RateClassMat, c(1, 2))
  StateMats <- list(RateCat1, RateCat2)
  FullMat <- getFullMat(StateMats, RateClassMat)
  return(FullMat)
}

# # setting up houwie
# houwie.model.setup <- function(corhmm_set, model_names="") {
#   # getting model structure for continuous trait
#   # setting up both character independent and OUMV models
#   cid_oumv_model <- full_oum_model <- full_ouv_model <- full_oumv_model <- getOUParamStructure("OUMV", 3, 2, null.model = TRUE) # character independent model (i.e. all observed states have the same optima)
#
#   # manually adjusting full OUMV model matrix (i.e. variable theta and sigma2 between observed characters)
#   full_oumv_model[2,c(1:6)] <- c(2:4)
#   full_oumv_model[3,c(1:6)] <- c(5:7)
#   # manually adjusting OUM model matrix
#   full_oum_model[3,c(1:6)] <- c(4:6)
#   # manually adjusting OUV model matrix
#   full_ouv_model[2,c(1:6)] <- c(2:4)
#
#   # All models that were run
#   model_list <- list(list(2, corhmm_set, "BM1"),
#                      list(2, corhmm_set, "OU1"),
#                      list(2, corhmm_set, cid_oumv_model),
#                      list(2, corhmm_set, full_ouv_model),
#                      list(2, corhmm_set, full_oum_model),
#                      list(2, corhmm_set, full_oumv_model)) # the two is for two rate classes
#
#   names(model_list) <- c(paste0("bm1_",model_names), paste0("ou1_",model_names),
#                          paste0("cid_oumv_",model_names), paste0("ouv_",model_names),
#                          paste0("oum_",model_names), paste0("oumv_",model_names))
#   return(model_list)
# }
#
#
# #  full houwie run
# one.full.houwie.run <- function(dat, phy, disc_model, model_names="", ncores=20){
#   # organize data
#   shared_species <- intersect(dat$tips, phy$tip.label)
#   dat <- dat[match(shared_species, dat$tips),]
#   phy <- keep.tip(phy, shared_species)
#   dat <- dat[match(phy$tip.label, dat$tips),]
#   # model structure
#   #corhmm_set <- corhmm.model.setup(dat[,c(1,2)])
#   corhmm_set <- disc_model
#   model_list <- houwie.model.setup(corhmm_set, model_names)
#   # setting up parallel
#   quickFunc <- function(model_list, model_name){
#     res <- hOUwie(phy, dat, model_list[[1]], model_list[[2]], model_list[[3]], nSim = 100, diagn_msg = TRUE, adaptive_sampling = FALSE, n_starts = 10, ncores = 10)
#     file.name <- paste0("houwie_results/",model_name, ".Rsave")
#     save(res, file=file.name)
#   }
#   mclapply(1:6, function(x) quickFunc(model_list[[x]], names(model_list)[x]), mc.cores = 10)
# }
