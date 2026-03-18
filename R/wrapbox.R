#' Return file name keys for MERIT Hydro rasters
#'
#' @inheritParams get_tf
#'
#' @section
#'  MERIT Hydro: see https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_Hydro/
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @importFrom dplyr %>%
#'
#' @export

get_key <- function(shape) {

  ## base data frame for lat/lon
  ## - longitude
  df_lon <- data.frame(min = c(-seq(5, 180, by = 5),
                               seq(0, 175, by = 5)),
                       max = c(-seq(0, 175, by = 5),
                               seq(5, 180, by = 5)))

  ## negative = west, positive = east
  ew <- with(df_lon, ifelse(min < 0, "w", "e"))
  numeric_x <- with(df_lon, stringr::str_pad(abs(min), 3, pad = "0"))

  df_lon <- df_lon %>%
    dplyr::mutate(code = paste0(ew, numeric_x)) %>%
    dplyr::arrange(min)

  ## - latitude
  df_lat <- data.frame(min = c(-seq(5, 90, by = 5),
                               seq(0, 85, by = 5)),
                       max = c(-seq(0, 85, by = 5),
                               seq(5, 90, by = 5)))

  ## negative = south, positive = north
  sn <- with(df_lat, ifelse(min < 0, "s", "n"))
  numeric_y <- with(df_lat, stringr::str_pad(abs(min), 2, pad = "0"))

  df_lat <- df_lat %>%
    dplyr::mutate(code = paste0(sn, numeric_y)) %>%
    dplyr::arrange(min)

  ## longitude
  xmin <- get_tf(data = df_lon,
                 shape = shape,
                 mode = "xmin")

  xmax <- get_tf(data = df_lon,
                 shape = shape,
                 mode = "xmax")

  ## latitude
  ymin <- get_tf(data = df_lat,
                 shape = shape,
                 mode = "ymin")

  ymax <- get_tf(data = df_lat,
                 shape = shape,
                 mode = "ymax")

  ## raster file key
  y <- df_lat$code[sort(which(ymin == 1):which(ymax == 1))]
  x <- df_lon$code[sort(which(xmin == 1):which(xmax == 1))]

  df_key <- expand.grid(y = y, x = x)
  z <- with(df_key, paste0(y, x))
  key <- paste(z, collapse = "|")

  return(key)
}

#' Convert ArcGIS flow direction to D8 flow direction
#'
#' @param x Flow direction raster of class \code{SpatRaster}
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

arc2d8 <- function(x) {
  ## whitebox uses flow pointer D8
  # 64, 128,  1
  # 32,   0,  2
  # 16,   8,  4

  ## ArcGIS uses D8 algorithm
  # 32, 64, 128
  # 16,  0,   1
  #  8,  4,   2

  # check class()
  if (!inherits(x, "SpatRaster")) stop("Provide data in class 'SpatRaster'")

  # begin with northeast through north
  fdir_arc <- as.integer(c(0, 2^(0:7), 247, 255))
  fdir_d8 <- as.integer(c(0, fdir_arc[3:9], fdir_arc[2], NA, NA))
  y0 <- terra::subst(x, from = fdir_arc, to = fdir_d8)

  # convert to integer just in case
  y <- terra::as.int(y0)

  return(y)
}

#' Convert raster to polygon(s)
#'
#' @param x Raster of class \code{SpatRaster}
#' @param simplify Logical.
#'  Whether output polygons are simplified or not.
#' @param keep Numeric.
#'  Proportion of vertices kept after polygon simplifications.
#'  Ignored if \code{simplify = FALSE}
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

rast2poly <- function(x,
                      simplify = TRUE,
                      keep = 0.05) {

  if (!(keep < 1 && keep > 0))
    stop("'keep' must be greater than 0 and less than 1")

  poly_raw <- stars::st_as_stars(x) %>%
    sf::st_as_sf(merge = TRUE,
                 as_points = FALSE) %>%
    sf::st_cast(to = "MULTIPOLYGON")

  if (simplify) {
    poly <- rmapshaper::ms_simplify(poly_raw,
                                    keep = keep)
  } else {
    poly <- poly_raw
  }

  return(poly)
}

#' Return UTM zone (EPSG code) based on geographic coordinates
#'
#' @param x Numeric scalar.
#'  Longitude coordinate.
#' @param y Numeric scalar.
#'  Latitude coordinate.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

get_utm <- function(x, y) {

  ## check inputs
  if (x < -180 || x >= 180)
    stop("Invalid value in x")
  if (y < -90 || y > 90)
    stop("Invalid value in y")

  ## UTM zone
  zone <- (floor((x + 180) / 6) %% 60) + 1

  ## Norway exception (zone 32V)
  if (y >= 56 && y < 64 && x >= 3 && x < 12)
    zone <- 32

  ## Svalbard exceptions (X band)
  if (y >= 72 && y < 84) {
    if (x >= 0  && x < 9)  zone <- 31
    if (x >= 9  && x < 21) zone <- 33
    if (x >= 21 && x < 33) zone <- 35
    if (x >= 33 && x < 42) zone <- 37
  }

  ## EPSG code
  if (y >= 0) {
    if (y <= 84) {
      epsg <- 32600 + zone
    } else {
      epsg <- 32661
    }
  } else {
    if (y >= -80) {
      epsg <- 32700 + zone
    } else {
      epsg <- 32761
    }
  }

  return(epsg)
}

#' Return UTM zone (EPSG code) based on a point layer
#'
#' @param point Point object of class \code{sf}.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

point2utm <- function(point) {

  if (!inherits(point, "sf"))
    stop("'point' must be class 'sf'")

  v_x <- sf::st_coordinates(point)[, 1]
  v_y <- sf::st_coordinates(point)[, 2]

  cout <- sapply(seq_len(length(v_x)), function(i) {
    get_utm(x = v_x[i], y = v_y[i])
  })

  return(cout)
}

#' Delineate unnested watersheds
#'
#' Delineate unnested watersheds for multiple outlet points using a
#' D8 flow direction raster.
#'
#' @param outlet Outlet point layer of class \code{sf}.
#' @param id_col Column name specifying outlet id.
#'  This column information will be appended to the output polygon layer.
#' @param f_dir Flow direction raster of class \code{SpatRaster}.
#' @param str_grid Stream grid raster of class \code{SpatRaster}.
#'  Required if \code{snap = TRUE}.
#' @param snap Logical. Whether snapping of outlet points should be performed.
#'  Defaults to \code{TRUE}.
#' @param snap_dist Numeric. Distance threshold for snapping points to stream grid.
#'  Measured in the unit of input raster files. Defaults to \code{5}.
#'
#' @return A polygon layer of class \code{sf} with one row per outlet point.
#'  If \code{id_col} is supplied, the output includes the identifier column and
#'  coordinates of the original (\code{x0}, \code{y0}) and snapped
#'  (\code{x}, \code{y}) outlet points.
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

wsd_unnested <- function(outlet,
                         id_col = NULL,
                         f_dir,
                         str_grid = NULL,
                         snap = TRUE,
                         snap_dist = 5) {


  # index outlet ------------------------------------------------------------

  ## this may not be necessary, but just in case
  outlet <- dplyr::mutate(outlet, idx = dplyr::row_number())

  # temporary files ---------------------------------------------------------

  message("Saving temporary files...")

  ## unique temp dir per call/worker
  temppath <- tempfile(
    pattern = paste0("wsd_", Sys.getpid(), "_")
  )
  dir.create(temppath)

  on.exit(
    unlink(temppath,
           recursive = TRUE,
           force = TRUE),
    add = TRUE,
    after = FALSE
  )

  ## setup temporary file names
  v_name <- file.path(temppath,
                      c("strg.tif",
                        "outlet.shp",
                        "outlet_snap.shp",
                        "dir.tif",
                        "wsd.tif")) %>%
    stats::setNames(c("strg",
                      "outlet",
                      "outlet_snap",
                      "dir",
                      "wsd"))

  ## write base raster files
  terra::writeRaster(f_dir,
                     filename = unname(v_name["dir"]),
                     overwrite = TRUE)

  sf::st_write(outlet,
               dsn = unname(v_name["outlet"]),
               append = FALSE,
               quiet = TRUE)

  # snapping outlets --------------------------------------------------------

  if (snap) {
    ## w/ snapping
    if (is.null(str_grid))
      stop("Stream grid 'str_grid' must be supplied if snap = TRUE")

    terra::writeRaster(str_grid,
                       filename = unname(v_name["strg"]),
                       overwrite = TRUE)

    message("Snap outlet points to the nearest stream grid...")
    whitebox::wbt_jenson_snap_pour_points(
      pour_pts = unname(v_name["outlet"]),
      streams = unname(v_name["strg"]),
      output = unname(v_name["outlet_snap"]),
      snap_dist = snap_dist,
      wd = temppath
    )

  } else {
    ## w/o snapping
    sf::st_write(outlet,
                 dsn = unname(v_name["outlet_snap"]),
                 quiet = TRUE,
                 append = FALSE)
  }

  # delineation -------------------------------------------------------------

  message("Delineate watersheds...")

  whitebox::wbt_unnest_basins(
    d8_pntr = unname(v_name["dir"]),
    pour_pts = unname(v_name["outlet_snap"]),
    output = unname(v_name["wsd"]),
    wd = temppath
  )

  # vectorize ---------------------------------------------------------------

  message("Vectorize raster watersheds...")

  ## read snapped outlets, then re-id with unique coordinates
  ## ordered as input outlets
  outlet_snap <- sf::st_read(dsn = unname(v_name["outlet_snap"])) %>%
    dplyr::select(.data$geometry) %>% # drop FID
    dplyr::mutate(idx = dplyr::row_number()) %>%
    dplyr::group_by(.data$geometry) %>%
    dplyr::mutate(pid = dplyr::cur_group_id()) %>%
    dplyr::ungroup()

  ## vectorize raster watersheds
  sf_wsd0 <- list.files(path = temppath,
                        pattern = "wsd.*\\.tif$",
                        full.names = TRUE) %>%
    lapply(terra::rast) %>%
    lapply(stars::st_as_stars) %>%
    lapply(sf::st_as_sf,
           merge = TRUE,
           as_points = FALSE) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(
      tifid = rowSums(dplyr::across(dplyr::ends_with("tif")),
                      na.rm = TRUE)
    ) %>%
    dplyr::select(.data$tifid) %>%
    dplyr::mutate(pid = outlet_snap$pid[.data$tifid])

  ## drop residual polygons, order by tifid (outlet id)
  sf_wsd <- sf_wsd0 %>%
    dplyr::mutate(area = units::set_units(sf::st_area(sf_wsd0), "km^2")) %>%
    dplyr::group_by(.data$pid) %>%
    dplyr::slice(which.max(.data$area)) %>% # remove duplicates by outlet
    dplyr::ungroup() %>%
    dplyr::relocate(.data$pid,
                    .data$tifid,
                    .data$area) %>%
    dplyr::arrange(.data$tifid)

  v_tifid <- sf_wsd$tifid

  ## subset by selected outlet, then extract coordinates
  ## - merging occurs when outlets are close to each other
  ## - 'idx' is original row ID, and 'tifid' should correspond to it
  xy0 <- outlet %>%
    dplyr::filter(.data$idx %in% v_tifid) %>%
    sf::st_coordinates()

  xy <- outlet_snap %>%
    dplyr::filter(.data$idx %in% v_tifid) %>%
    sf::st_coordinates()

  ## append outlet coordinates
  sf_wsd <- sf_wsd %>%
    dplyr::mutate(x = xy[, 1],
                  y = xy[, 2],
                  x0 = xy0[, 1],
                  y0 = xy0[, 2],
                  .before = .data$geometry)

  if (!is.null(id_col)) {

    v_sid <- outlet %>%
      dplyr::filter(.data$idx %in% v_tifid) %>%
      dplyr::pull(id_col)

    sf_wsd <- sf_wsd %>%
      dplyr::mutate(!!id_col := v_sid,
                    .before = .data$tifid)

  }

  return(sf_wsd)
}

#' Delineate nested watersheds
#'
#' Delineate nested watersheds for multiple outlet points using a D8 flow
#' direction raster.
#'
#' @inheritParams wsd_unnested
#'
#' @return A polygon layer of class \code{sf} with one row per outlet point.
#'  If \code{id_col} is supplied, the output includes the identifier column and
#'  coordinates of the original (\code{x0}, \code{y0}) and snapped
#'  (\code{x}, \code{y}) outlet points.
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

wsd_nested <- function(outlet,
                       id_col = NULL,
                       f_dir,
                       str_grid = NULL,
                       snap = TRUE,
                       snap_dist = 5) {

  # index outlet ------------------------------------------------------------

  ## this may not be necessary, but just in case
  outlet <- dplyr::mutate(outlet, idx = dplyr::row_number())

  # temporary files ---------------------------------------------------------

  message("Saving temporary files...")

  ## unique temp dir per call/worker
  temppath <- tempfile(pattern = paste0("wsd_", Sys.getpid(), "_"))
  dir.create(temppath)

  on.exit(
    unlink(temppath,
           recursive = TRUE,
           force = TRUE),
    add = TRUE,
    after = FALSE
  )

  ## setup temporary file names
  v_name <- file.path(temppath,
                      c("strg.tif",
                        "outlet.shp",
                        "outlet_snap.shp",
                        "dir.tif",
                        "wsd.tif")) %>%
    stats::setNames(c("strg",
                      "outlet",
                      "outlet_snap",
                      "dir",
                      "wsd"))

  ## write temporary files
  terra::writeRaster(f_dir,
                     filename = unname(v_name["dir"]),
                     overwrite = TRUE)

  # snapping ----------------------------------------------------------------

  if (snap) {
    ## w/ snapping
    message("Snap outlet points to the nearest stream grid...")

    if (is.null(str_grid))
      stop("Stream grid 'str_grid' must be supplied if snap = TRUE")

    terra::writeRaster(
      str_grid,
      filename = unname(v_name["strg"]),
      overwrite = TRUE
    )

    sf::st_write(
      outlet,
      dsn = unname(v_name["outlet"]),
      append = FALSE,
      quiet = TRUE
    )

    whitebox::wbt_jenson_snap_pour_points(
      pour_pts = unname(v_name["outlet"]),
      streams = unname(v_name["strg"]),
      output = unname(v_name["outlet_snap"]),
      snap_dist = snap_dist,
      wd = temppath
    )

  } else {

    ## w/o snapping
    sf::st_write(
      outlet,
      unname(v_name["outlet_snap"]),
      append = FALSE,
      quiet = TRUE
    )

  }

  # delineation -------------------------------------------------------------

  message("Delineate watersheds...")

  whitebox::wbt_watershed(
    d8_pntr = unname(v_name["dir"]),
    pour_pts = unname(v_name["outlet_snap"]),
    output = unname(v_name["wsd"]),
    wd = temppath
  )

  # vectorize ---------------------------------------------------------------

  message("Vectorize raster watersheds...")

  outlet_snap <- sf::st_read(dsn = unname(v_name["outlet_snap"])) %>%
    dplyr::select(.data$geometry) %>%
    dplyr::mutate(idx = dplyr::row_number()) %>%
    dplyr::group_by(.data$geometry) %>%
    dplyr::mutate(pid = dplyr::cur_group_id()) %>%
    dplyr::ungroup()

  sf_wsd0 <- terra::rast(unname(v_name["wsd"])) %>%
    stars::st_as_stars() %>%
    sf::st_as_sf(merge = TRUE,
                 as_points = FALSE) %>%
    dplyr::rename(tifid = .data$wsd.tif) %>%
    dplyr::mutate(pid = outlet_snap$pid[.data$tifid])

  sf_wsd <- sf_wsd0 %>%
    dplyr::mutate(area = units::set_units(sf::st_area(sf_wsd0), "km^2")) %>%
    dplyr::group_by(.data$tifid) %>%
    dplyr::slice(which.max(.data$area)) %>%
    dplyr::ungroup() %>%
    dplyr::relocate(.data$pid,
                    .data$tifid,
                    .data$area)

  v_tifid <- sf_wsd$tifid

  ## subset by selected outlet, then extract coordinates
  ## - merging occurs when outlets are close to each other
  ## - 'idx' is original row ID, and 'tifid' should correspond to it
  xy0 <- outlet %>%
    dplyr::filter(.data$idx %in% v_tifid) %>%
    sf::st_coordinates()

  xy <- outlet_snap %>%
    dplyr::filter(.data$idx %in% v_tifid) %>%
    sf::st_coordinates()

  ## append outlet coordinates
  sf_wsd <- sf_wsd %>%
    dplyr::mutate(x = xy[, 1],
                  y = xy[, 2],
                  x0 = xy0[, 1],
                  y0 = xy0[, 2],
                  .before = .data$geometry)

  if (!is.null(id_col)) {
    ## get unique outlet identifier
    v_sid <- outlet %>%
      dplyr::filter(.data$idx %in% v_tifid) %>%
      dplyr::pull(id_col)

    sf_wsd <- sf_wsd %>%
      dplyr::mutate(!!id_col := v_sid,
                    .before = .data$tifid)

  }

  return(sf_wsd)
}

#' Convert flow accumulation raster to stream grid
#'
#' Extracts a stream grid from a flow accumulation raster by applying a
#' minimum drainage area threshold via \code{whitebox::wbt_extract_streams()}.
#'
#' @param f_acc Flow accumulation raster of class \code{SpatRaster}.
#' @param threshold Numeric. Minimum drainage area threshold for stream
#'  initiation. The unit inherits from the flow accumulation layer.
#' @param output Character. File path for the output stream grid raster.
#'
#' @return A stream grid raster of class \code{SpatRaster}.
#'
#' @importFrom dplyr %>%
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

flow2grid <- function(f_acc,
                      threshold,
                      output) {

  # Use a truly unique temp dir per call (safe across workers)
  temppath <- tempfile(pattern = paste0("strg_", Sys.getpid(), "_"))
  dir.create(temppath)

  on.exit(
    unlink(temppath,
           recursive = TRUE,
           force = TRUE),
    add = TRUE,
    after = FALSE
  )

  fname <- file.path(temppath, "upa.tif")

  # Write raster from file path, not SpatRaster object (avoids fork issues)
  terra::writeRaster(f_acc,
                     filename = fname,
                     overwrite = TRUE)

  # Ensure output paths are unique per call (caller's responsibility)
  whitebox::wbt_extract_streams(
    flow_accum = fname,
    output     = output,
    threshold  = threshold,
    wd         = temppath
  )

  return(terra::rast(output))
}

#' Convert stream grid to vector stream
#'
#' Converts a raster stream grid to a vector line layer using
#' \code{whitebox::wbt_raster_streams_to_vector()}.
#'
#' @inheritParams wsd_unnested
#' @param output Character. Optional file path to save the output stream vector.
#'  If \code{NULL} (default), the result is returned in memory only.
#' @param set_crs Logical. Whether the output should inherit the CRS from
#'  \code{f_dir}. Defaults to \code{TRUE}.
#'
#' @return A vector line layer of class \code{sf}.
#'
#' @importFrom dplyr %>%
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

grid2stream <- function(f_dir,
                        str_grid,
                        output = NULL,
                        set_crs = TRUE) {

  ## Use a truly unique temp dir per call (safe across workers)
  temppath <- tempfile(pattern = paste0("strv_", Sys.getpid(), "_"))
  dir.create(temppath)

  fname <- file.path(temppath,
                     c("dir.tif",
                       "strg.tif",
                       "strv.shp")) %>%
    stats::setNames(c("dir",
                      "strg",
                      "strv"))

  on.exit(
    unlink(temppath,
           recursive = TRUE,
           force = TRUE),
    add = TRUE,
    after = FALSE
  )

  ## write raster input in temporary folder
  terra::writeRaster(f_dir,
                     filename = unname(fname["dir"]),
                     overwrite = TRUE)

  terra::writeRaster(str_grid,
                     filename = unname(fname["strg"]),
                     overwrite = TRUE)

  ## stream grids
  whitebox::wbt_raster_streams_to_vector(streams = unname(fname["strg"]),
                                         d8_pntr = unname(fname["dir"]),
                                         output = unname(fname["strv"]),
                                         wd = temppath)

  channel <- sf::st_read(unname(fname["strv"]), quiet = TRUE)

  ## inherit CRS from source
  if (set_crs) {
    channel <- sf::st_set_crs(
      channel,
      terra::crs(f_dir)
    )
  }

  if (!is.null(output))
    sf::st_write(channel,
                 dsn = output,
                 append = FALSE,
                 quiet = TRUE)

  return(channel)
}
