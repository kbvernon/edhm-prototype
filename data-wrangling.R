library(arcgislayers)
library(dplyr)
library(httr2)
library(sf)
library(terra)

gpkg <- file.path("data", "edhm.gpkg")

# CRS
epsg <- 3857

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
# need to validate geometries to ensure that states with islands
# are represented as MULTIPOLYGON, not POLYGON
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
