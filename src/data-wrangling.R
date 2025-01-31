library(arcgislayers)
library(dplyr)
library(httr2)
library(sf)
library(terra)

# path to geopackage database
gpkg <- file.path("data", "edhm.gpkg")

# CRS
epsg <- 3857

# time frame
tf <- as.Date(c(start = "2011-01-01", end = "2020-12-31"))

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

# write to geopackage (sqlite) database
write_sf(
  western_us,
  dsn = gpkg,
  layer = "window"
)

remove(service, states, targets)

# soil moisture station metadata -----------------------------------------

# get SCAN and SNOTEL stations from USDA National Water and Climate Center
service <- file.path(
  "https://wcc.sc.egov.usda.gov",
  "awdbRestApi",
  "services/v1",
  "stations"
) |> request()

# api wants triplets of the form stationId:stateCode:networkCode
# '*' wildcard also allowed
targets <- sprintf(
  "*:%1s:SCAN, *:%1s:SNTL",
  c("AZ", "CA", "CO", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY"),
  c("AZ", "CA", "CO", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY")
) |> paste(collapse = ",")

stations <- service |>
  req_url_query(stationTriplets = targets) |>
  req_perform() |>
  resp_body_json(simplifyVector = TRUE) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(epsg) |>
  rename(
    "id" = stationId,
    "state" = stateCode,
    "network" = networkCode
  ) |>
  select(id, state, network)

# get stations from US Climate Reference Network run by NOAA
# - can find spatial metadata in the NOAA NCEI Historical Observing
#   Metadata Repository, but it's a poorly structured text file
# - so we grab the table put together by the National Soil Moisture Network and
#   subset to the USCRN stations in the project area
# - who knows how long this link will work...
path <- file.path(
  "http://nationalsoilmoisture.com",
  "test/VWC_QAQC",
  "Stationlist.csv"
)

uscrn <- read.csv(path) |>
  rename_with(tolower) |>
  filter(network == "USCRN") |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(epsg) |>
  select(id, state, network) |>
  st_filter(western_us)

# primary key = id:state:network
stations <- rbind(stations, uscrn)

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
  layer = "sm-station-metadata"
)

remove(service, targets, path, uscrn)

# soil moisture station data ---------------------------------------------
