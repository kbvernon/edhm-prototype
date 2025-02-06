library(arcgislayers)
library(awdb)
library(dplyr)
library(httr2)
library(sf)
library(terra)
library(tidyr)

# path to geopackage (sqlite) database
gpkg <- file.path("data", "edhm.gpkg")

# CRS
epsg <- 3857

# time frame (ten years and one month)
# we include extra month at beginning as buffer for lag effects
time_frame <- list(start = "2010-12-01", end = "2020-12-31")

# project window ---------------------------------------------------------

# use US Census TIGERweb AGOL REST API to get state boundaries
service <- file.path(
  "https://tigerweb.geo.census.gov",
  "arcgis/rest/services",
  "TIGERweb",
  "State_County",
  "MapServer/0"
) |> arc_open()

targets <- sprintf(
  "'%s'",
  c(
    "Arizona",
    "California",
    "Colorado",
    "Idaho",
    "Montana",
    "Nevada",
    "New Mexico",
    "Oregon",
    "Utah",
    "Washington",
    "Wyoming"
  )
) |> paste(collapse = ",")

states <- arc_select(
  service,
  fields = "NAME",
  where = paste0("NAME IN (", targets, ")")
)

# merge state boundaries and remove islands off coast
# - need to validate geometries to ensure that states with islands
#   are represented as MULTIPOLYGON, not POLYGON
# - merge, cast back to POLYGON, and pull the largest one to exclude islands
western_us <- states |>
  st_make_valid() |>
  summarize() |>
  st_cast("POLYGON") |>
  slice(which.max(st_area(geometry))) |>
  st_geometry()

# plot(western_us)

write_sf(
  western_us,
  dsn = gpkg,
  layer = "window"
)

remove(service, states, targets)

# soil moisture station metadata -----------------------------------------

# primary key = id:state:network

# get stations from USDA NWCC AWDB REST API that are in the aoi and have soil
# moisture data
stations <- western_us |>
  st_geometry() |>
  get_station_metadata(elements = "SMS:*") |>
  rename(
    "id" = station_id,
    "state" = state_code,
    "network" = network_code
  ) |>
  select(station_triplet, id, state, network, name) |>
  st_transform(epsg)

# TODO: get USCRN data

# # get stations from US Climate Reference Network run by NOAA
# # - can find spatial metadata in the NOAA NCEI Historical Observing
# #   Metadata Repository, but it's a poorly structured text file
# # - so we grab the table put together by the National Soil Moisture Network and
# #   subset to the USCRN stations in the project area
# # - who knows how long this link will work...
# path <- file.path(
#   "http://nationalsoilmoisture.com",
#   "test/VWC_QAQC",
#   "Stationlist.csv"
# )

# # "https://gis.ncdc.noaa.gov/arcgis/rest/services/cdo/crn/MapServer/0"

# uscrn <- read.csv(path) |>
#   rename_with(tolower) |>
#   filter(network == "USCRN") |>
#   st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
#   st_transform(epsg) |>
#   select(id, state, network) |>
#   st_filter(western_us)

# stations <- rbind(stations, uscrn)

# plot(western_us)
# plot(
#   st_geometry(stations),
#   pch = 19,
#   cex = 0.7,
#   col = adjustcolor("gray35", alpha.f = 0.7),
#   add = TRUE
# )

write_sf(
  stations,
  dsn = gpkg,
  layer = "stations"
)

# remove(path, uscrn)

# soil moisture station data ---------------------------------------------

# SMS:* - soil moisture percent, all depths
station_data <- western_us |>
  get_station_data(
    elements = "SMS:*",
    begin_date = time_frame[["start"]],
    end_date = time_frame[["end"]],
    request_size = 35L
  ) |>
  select(station_triplet:duration_name, values) |>
  unnest(values)

write_sf(
  station_data,
  dsn = gpkg,
  layer = "soil"
)

remove(station_data)

# elevation and terrain --------------------------------------------------

# location of the virtual raster pointing to the
# USGS one arc second elevation data (~30 m resolution)
service <- file.path(
  "https://prd-tnm.s3.amazonaws.com/StagedProducts",
  "Elevation/1/TIFF",
  "USGS_Seamless_DEM_1.vrt"
)

# connect to vrt
rr <- rast(service, vsi = TRUE)

terrain_data <- vector(mode = "list", length = nrow(stations))

for (i in 1:nrow(stations)) {
  feature <- stations[i, ] |>
    st_buffer(1000) |>
    vect() |>
    project(crs(rr))

  focal_rr <- crop(rr, ext(feature), snap = "out")

  # compute terrain measures and stack
  rr_stack <- c(
    focal_rr,
    terrain(focal_rr, "TRIrmsd"),
    terrain(focal_rr, "TPI"),
    terrain(focal_rr, "slope", unit = "degrees") <= 5
  )

  names(rr_stack) <- c("elevation", "tri_rmsd", "tpi", "p_flat")

  terrain_data[[i]] <- terra::extract(
    rr_stack,
    feature,
    fun = mean,
    weights = TRUE,
    na.rm = TRUE,
    ID = FALSE
  )

  # simple progress bar
  p <- i / nrow(stations)

  sprintf(
    "%-50s| %3s%% \r",
    strrep("=", floor(50 * p)),
    floor(100 * p)
  ) |> cat()

  flush.console()

  remove(i, feature, focal_rr, rr_stack, p)
}

terrain_data <- do.call(rbind, terrain_data)

terrain_data <- cbind(
  station_triplet = stations[["station_triplet"]],
  terrain_data
)

terrain_data <- as.data.frame(terrain_data)

write_sf(
  terrain_data,
  dsn = gpkg,
  layer = "terrain"
)

remove(service, rr, terrain_data)

# gridmet ----------------------------------------------------------------

service <- file.path(
  "https://www.northwestknowledge.net",
  "metdata",
  "data"
)

variables <- c("pr", "srad", "tmmn", "tmmx")

years <- 2010:2020

features <- stations |>
  st_geometry() |>
  vect() |>
  project("epsg:4326")

gridmet_data <- vector(mode = "list", length = length(variables))

# takes nearly 40 minutes on my machine
for (i in seq_along(variables)) {
  layers <- file.path(
    service,
    paste0(variables[i], "_", years, ".nc")
  )

  r <- rast(layers, win = ext(features))

  # data matrix: nrows = n stations, ncols = n dates
  M <- terra::extract(r, features, ID = FALSE, raw = TRUE)

  # transpose and unfold M to get long vector with
  # all dates for each station as contiguous elements
  gridmet_data[[i]] <- c(t(M))

  cat("Data for", variables[i], "variable successfuly extracted.\n")

  remove(i, layers, r, M)
}

gridmet_data <- do.call(cbind, gridmet_data)
gridmet_data <- as.data.frame(gridmet_data)

names(gridmet_data) <- c(
  "precipitation",
  "radiation",
  "temperature_min",
  "temperature_max"
)

dates <- seq(as.Date("2010-01-01"), as.Date("2020-12-31"), by = 1)

keys <- data.frame(
  station_triplet = rep(stations[["station_triplet"]], each = length(dates)),
  date = rep(dates, each = nrow(stations))
)

gridmet_data <- cbind(keys, gridmet_data)

# subset to time period of interest
i <- which(gridmet_data[["date"]] >= as.Date(time_frame[["start"]]))

gridmet_data <- gridmet_data[i, ]
gridmet_data[["date"]] <- as.character(gridmet_data[["date"]])

write_sf(
  gridmet_data,
  dsn = gpkg,
  layer = "climate"
)

remove(service, variables, years, features, dates, keys, i, gridmet_data)
