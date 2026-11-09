library(tidyverse)
library(sf)
library(mapview)

# raster example ----------------------------------------------------------

wsd <- st_read("tools/hybas_na_lev07_v1c.shp") %>%
  filter(HYBAS_ID == "7070643030")

mapview(wsd, map.types = "OpenStreetMap")

file.path(normalizePath("."))

path0 <- here::here()
setwd("C://Users//a_terui//OneDrive - UNCG//Documents//shared_folder")

key <- wrapbox::get_key(wsd)
list_ras <- list.files(
  "gis",
  pattern = paste0(key, c("_upa", "_dir"), collapse = "|"),
  full.names = TRUE
) %>%
  lapply(terra::rast)

setwd(path0)

mask <- wsd %>%
  st_transform(crs = 32617) %>%
  st_buffer(dist = 100) %>%
  st_transform(crs = 4326)

lapply(list_ras,
       function(x) {
         y <- terra::crop(x, mask) %>%
           terra::mask(mask)

         terra::writeRaster(
           y,
           filename = paste0("inst/raster/",
                             str_extract(names(y), "upa|dir"),
                             ".tif"),
           overwrite = TRUE
         )
       })


# point example -----------------------------------------------------------

df_site <- tibble(lon = c(-80.51863, -80.80812), lat = c(35.78029, 36.02860))

sf_site <- st_as_sf(df_site, coords = c("lon", "lat"), crs = 4326) %>%
  mutate(id = row_number(), .before = "geometry")

saveRDS(sf_site, "inst/vector/outlet.rds")


