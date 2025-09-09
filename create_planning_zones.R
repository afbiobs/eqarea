# create_planning_zones.R
# Construct planning zones oriented by terrain aspect from a DEM.

library(terra)
library(sf)
library(dplyr)
library(lwgeom)

# 1. Load DEM in EPSG:2157
# Update the path as necessary
raster_path <- "NI_ZONE_PROJ2157.tif"
dem <- rast(raster_path)

# Fill small no-data holes using a 3x3 focal mean
fill_na <- function(r) {
  filled <- focal(r, w = matrix(1, 3, 3), fun = mean, na.policy = "only", na.rm = TRUE)
  r[is.na(r)] <- filled[is.na(r)]
  r
}

dem <- fill_na(dem)

# 2. Derive slope and aspect rasters
terrain_stack <- terrain(dem, v = c("slope", "aspect"), unit = "degrees")
slope <- terrain_stack$slope
aspect <- terrain_stack$aspect

# Smooth aspect to reduce noise
aspect <- focal(aspect, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)

# 3. Create a 10 km reference grid covering the DEM extent
ext_polygon <- st_as_sfc(ext(dem))
st_crs(ext_polygon) <- 2157
grid <- st_make_grid(ext_polygon, cellsize = c(10000, 10000))
grid <- st_sf(id = seq_along(grid), geometry = grid)

# 4. Compute mean elevation and circular mean aspect per grid square
mean_elev <- terra::extract(dem, vect(grid), fun = mean, na.rm = TRUE)

circular_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  r <- x * pi / 180
  m <- atan2(mean(sin(r)), mean(cos(r))) * 180 / pi
  (m + 360) %% 360
}

asp_vals <- terra::extract(aspect, vect(grid))
mean_aspect <- vapply(asp_vals, circular_mean, numeric(1))

grid$mean_elev <- mean_elev[,2]
grid$mean_aspect <- mean_aspect

# 5. Rectangle orientation and width
# Orientation runs parallel to contours (aspect + 90) modulo 180
# Width obeys min(mean_elev * 100, 5000)

grid$orientation <- (grid$mean_aspect + 90) %% 180
grid$width <- pmin(grid$mean_elev * 100, 5000)

grid <- grid[order(-grid$mean_elev), ]

# Helper: create an oriented rectangle of 10 km x width
create_rect <- function(center, orient_deg, width_m, crs) {
  base <- st_as_sfc(st_bbox(c(
    xmin = -5000, xmax = 5000,
    ymin = -width_m/2, ymax = width_m/2
  ), crs = crs))
  lwgeom::st_rotate(base, orient_deg * pi / 180, center)
}

# 6. Process rectangles sequentially, removing overlaps
accum <- st_sfc(crs = 2157)
rects <- list()
for (i in seq_len(nrow(grid))) {
  center <- st_centroid(grid[i,])
  rect <- create_rect(center, grid$orientation[i], grid$width[i], 2157)
  if (length(accum) > 0) {
    rect <- st_difference(rect, st_union(accum))
  }
  if (!st_is_empty(rect)) {
    rects[[length(rects) + 1]] <- rect
    accum <- c(accum, rect)
  }
}

zones <- st_sf(zone_id = seq_along(rects), geometry = st_sfc(rects, crs = 2157))

# Remove slivers and validate geometry
zones <- st_make_valid(zones)
zones <- st_buffer(zones, 0)

# 7. Export the zone layer to GeoPackage
st_write(zones, "planning_zones.gpkg", layer = "zones", driver = "GPKG", delete_dsn = TRUE)
