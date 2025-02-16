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

# stations ---------------------------------------------------------------

# primary key = id:state:network

# get stations from USDA NWCC AWDB REST API that are in the aoi and have soil
# moisture data
stations <- western_us |>
  st_geometry() |>
  get_stations(elements = "SMS:*") |>
  select(station_triplet, name)

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

# soil moisture ----------------------------------------------------------

# SMS:* - soil moisture percent, all depths
soil_data <- western_us |>
  get_elements(
    elements = "SMS:*",
    awdb_options = set_options(
      begin_date = time_frame[["start"]],
      end_date = time_frame[["end"]],
      request_size = 35L
    )
  ) |>
  unnest(element_values) |>
  select(station_triplet, height_depth, date, value) |>
  mutate(id = row_number()) |>
  pivot_wider(
    id_cols = c(id, station_triplet, date),
    names_from = height_depth,
    values_from = value
  ) |>
  select(-id) |>
  rename_with(gsub, pattern = "-", replacement = "sm_")

write_sf(
  soil_data,
  dsn = gpkg,
  layer = "soil"
)

remove(soil_data)

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

terrain_data <- as_tibble(terrain_data)

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

# pr   = precipitation
# sph  = specific humidity
# srad = solar radiation
# tmmn = temperature minimum
# tmmx = temperature maximum
# vs   = wind speed
variables <- c("pr", "sph", "srad", "tmmn", "tmmx", "vs")

years <- 2010:2020

features <- stations |>
  st_geometry() |>
  vect() |>
  project("epsg:4326")

climate_data <- vector(mode = "list", length = length(variables))

# each variable takes ~15-20 minutes to process
for (i in seq_along(variables)) {
  start <- Sys.time()

  layers <- file.path(
    service,
    paste0(variables[i], "_", years, ".nc")
  )

  r <- rast(layers, win = ext(features))

  # data matrix: nrows = n stations, ncols = n dates
  M <- terra::extract(r, features, ID = FALSE, raw = TRUE)

  # transpose and unfold M to get long vector with
  # all dates for each station as contiguous elements
  climate_data[[i]] <- c(t(M))

  end <- Sys.time()

  sprintf(
    "Successfully downloaded: %s. Total elapsed time: %s minutes.\n",
    variables[i],
    round(difftime(end, start, units = "mins"), 2)
  ) |> cat()

  remove(i, start, layers, r, M, end)
}

climate_data <- do.call(cbind, climate_data)
climate_data <- as_tibble(climate_data)

names(climate_data) <- c(
  "precipitation",
  "humidity",
  "radiation",
  "temperature_min",
  "temperature_max",
  "wind_speed"
)

dates <- seq(as.Date("2010-01-01"), as.Date("2020-12-31"), by = 1)

keys <- data.frame(
  station_triplet = rep(stations[["station_triplet"]], each = length(dates)),
  date = rep(dates, nrow(stations))
)

climate_data <- cbind(keys, climate_data)

# subset to time period of interest
i <- which(climate_data[["date"]] >= as.Date(time_frame[["start"]]))

climate_data <- climate_data[i, ]
climate_data[["date"]] <- as.character(climate_data[["date"]])

write_sf(
  climate_data,
  dsn = gpkg,
  layer = "climate"
)

remove(service, variables, years, features, dates, keys, climate_data)

# land cover -------------------------------------------------------------

service <- "https://s3-us-west-2.amazonaws.com/mrlc"

years <- 2010:2020

layers <- paste0(service, "/Annual_NLCD_LndCov_", years, "_CU_C1V0.tif")

# degenerate crs
target_crs <- crs(rast(layers[1], vsi = TRUE))

aoi <- western_us |>
  st_transform(target_crs) |>
  vect() |>
  ext()

virtual_raster <- rast(layers, vsi = TRUE, win = aoi)

names(virtual_raster) <- years

features <- stations |>
  st_transform(target_crs) |>
  st_buffer(units::set_units(4000, m)) |>
  vect()

Mode <- function(x) {
  tbl <- table(x, useNA = "no")
  as.integer(names(tbl)[which.max(tbl)])
}

land_data <- terra::extract(
  virtual_raster,
  features,
  fun = Mode,
  ID = FALSE
)

land_data <- land_data |>
  mutate(
    station_triplet = stations[["station_triplet"]],
    .before = everything()
  ) |>
  pivot_longer(
    -station_triplet,
    names_to = "year",
    values_to = "land_cover"
  ) |>
  mutate(
    year = as.integer(year),
    land_cover = case_match(
      land_cover,
      11 ~ 1L,
      12 ~ 2L,
      21 ~ 3L,
      22 ~ 4L,
      23 ~ 5L,
      24 ~ 6L,
      31 ~ 7L,
      41 ~ 8L,
      42 ~ 9L,
      43 ~ 10L,
      52 ~ 11L,
      71 ~ 12L,
      81 ~ 13L,
      82 ~ 14L,
      90 ~ 15L,
      95 ~ 16L,
      .default = NA_integer_
    )
  )

land_cover_lookup <- data.frame(
  id = as.integer(c(1:16, NA_integer_)),
  nlcd_id = c(
    11, 12, 21, 22, 23, 24, 31, 41,
    42, 43, 52, 71, 81, 82, 90, 95, 250
  ),
  name = c(
    "Open Water",
    "Perennial Ice/Snow",
    "Developed, Open Space",
    "Developed, Low Intensity",
    "Developed, Medium Intensity",
    "Developed, High Intensity",
    "Barren Land (Rock/Sand/Clay)",
    "Deciduous Forest",
    "Evergreen Forest",
    "Mixed Forest",
    "Shrub/Scrub",
    "Grassland/Herbaceous",
    "Pasture/Hay",
    "Cultivated Crops",
    "Woody Wetlands",
    "Emergent Herbaceous Wetlands",
    "No Data"
  )
)

write_sf(
  land_cover_lookup,
  dsn = gpkg,
  layer = "land_cover_lookup"
)

write_sf(
  land_data,
  dsn = gpkg,
  layer = "land"
)

remove(
  service, years, layers, target_crs, aoi, virtual_raster,
  features, Mode, land_data, land_cover_lookup
)

# table ------------------------------------------------------------------

coordinates <- stations |>
  st_coordinates() |>
  as_tibble() |>
  rename_with(tolower) |>
  mutate(
    station_triplet = stations[["station_triplet"]],
    .before = everything()
  )

soil_data <- read_sf(gpkg, "soil")
terrain_data <- read_sf(gpkg, "terrain")
climate_data <- read_sf(gpkg, "climate")
land_data <- read_sf(gpkg, "land")

all_data <- soil_data |>
  left_join(coordinates, by = "station_triplet") |>
  full_join(
    climate_data |> filter(station_triplet %in% soil_data[["station_triplet"]]),
    by = c("station_triplet", "date")
  ) |>
  select(where(\(x) sum(!is.na(x)) / length(x) >= 0.25)) |>
  filter(!is.na(sm_8)) |>
  left_join(terrain_data, by = "station_triplet") |>
  mutate(year = as.integer(strsplit(date, "-")[[1]][1])) |>
  left_join(land_data, by = c("station_triplet", "year")) |>
  select(
    station_triplet,
    date,
    x, y,
    sm_8,
    precipitation,
    humidity,
    radiation,
    temperature_min,
    temperature_max,
    wind_speed,
    elevation,
    tri_rmsd,
    tpi,
    p_flat,
    land_cover
  ) |>
  filter(!is.na(sm_8))

all_data |> readr::write_csv("data/edhm-prototype-data.csv")
