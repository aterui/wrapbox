# tests/testthat/test-watershed.R

test_that("nested and unnested watershed generation works", {

  library(terra)
  library(sf)
  library(dplyr)

  # ------------------------------------------------------------------
  # Load and preprocess test rasters
  # ------------------------------------------------------------------

  list_spr <- list.files(
    "tests",
    recursive = TRUE,
    pattern = "\\.tif$",
    full.names = TRUE
  ) %>%
    lapply(function(x) {
      tmp <- rast(x)
      r <- st_bbox(tmp)

      xmin <- r["xmin"]
      xmax <- r["xmin"] + abs(r["xmin"] - r["xmax"]) / 4
      ymin <- r["ymin"]
      ymax <- r["ymin"] + abs(r["ymin"] - r["ymax"]) / 4

      crop(tmp, ext(xmin, xmax, ymin, ymax))
    })

  spr_dir <- arc2d8(list_spr[[1]])
  spr_upa <- list_spr[[2]]

  # ------------------------------------------------------------------
  # Create outlet points
  # ------------------------------------------------------------------

  set.seed(133)

  valid_cells <- which(values(spr_upa) > 50)
  cell_id <- sample(valid_cells, 1)
  cell_id <- c(cell_id, cell_id + 10, cell_id, cell_id + 1, cell_id + 2)

  outlet <- xyFromCell(spr_upa, cell_id) %>%
    as.data.frame() %>%
    st_as_sf(coords = c("x", "y"), crs = 4326) %>%
    mutate(site = paste0("s", row_number()))

  expect_s3_class(outlet, "sf")

  # ------------------------------------------------------------------
  # Create stream grid
  # ------------------------------------------------------------------

  tmp_strg <- tempfile(fileext = ".tif")

  spr_strg <- flow2grid(
    f_acc = spr_upa,
    threshold = 50,
    output = tmp_strg
  )

  expect_true(file.exists(tmp_strg))

  on.exit(
    unlink(tmp_strg),
    add = TRUE
  )

  # ------------------------------------------------------------------
  # Run watershed functions
  # ------------------------------------------------------------------

  sf_wsd_n <- wsd_nested(
    outlet = outlet,
    f_dir = spr_dir,
    str_grid = spr_strg,
    snap = TRUE,
    id_col = "site"
  )

  sf_wsd_un <- wsd_unnested(
    outlet = outlet,
    f_dir = spr_dir,
    str_grid = spr_strg,
    snap = TRUE,
    id_col = "site"
  )

  # ------------------------------------------------------------------
  # Assertions
  # ------------------------------------------------------------------

  expect_s3_class(sf_wsd_n, "sf")
  expect_s3_class(sf_wsd_un, "sf")

  expect_gt(nrow(sf_wsd_n), 0)
  expect_gt(nrow(sf_wsd_un), 0)

  expect_true("site" %in% names(sf_wsd_n))
  expect_true("site" %in% names(sf_wsd_un))

  # nested should be >= unnested in count OR at least not empty
  expect_gte(nrow(sf_wsd_n), nrow(sf_wsd_un))

})
