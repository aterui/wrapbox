pacman::p_load(terra,
               sf,
               tidyverse,
               mapview)

list_spr <- list.files("tests",
                       recursive = TRUE,
                       pattern = "\\.tif",
                       full.names = TRUE) %>%
  lapply(FUN = function(x) {
    tmp <- rast(x)
    r <- st_bbox(tmp)

    xmin <- r["xmin"]
    xmax <- r["xmin"] + abs(r["xmin"] - r["xmax"])/4
    ymin <- r["ymin"]
    ymax <- r["ymin"] + abs(r["ymin"] - r["ymax"])/4

    crop(tmp, ext(xmin, xmax, ymin, ymax))
  })

spr_dir <- arc2d8(list_spr[[1]])
spr_upa <- list_spr[[2]]

cell_id <- sample(which(values(spr_upa) > 50), 1)
cell_id <- c(cell_id, cell_id, cell_id + 1)

# convert to coordinates
outlet <- xyFromCell(spr_upa, cell_id) %>%
  data.frame() %>%
  st_as_sf(coords = c("x", "y"),
           crs = 4326)

# stream grid
spr_strg <- flow2grid(f_acc = spr_upa,
                      threshold = 50,
                      output = "tests/testthat/testdata/strg.tif")

# nested watershed
sf_wsd_n <- wsd_nested(outlet = outlet,
                       f_dir = spr_dir,
                       str_grid = spr_strg,
                       snap = TRUE)

# unnested watershed
sf_wsd_un <- wsd_unnested(outlet = outlet,
                          f_dir = spr_dir,
                          str_grid = spr_strg,
                          snap = TRUE)

mapview(spr_strg) +
  mapview(outlet) +
  mapview(sf_wsd_un[1,])
