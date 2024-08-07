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
  if (x < -180 || x > 180)
    stop("Invalid value in x")

  if (y < -90 || y > 90)
    stop("Invalid value in y")

  ## x - longitude
  x0 <- seq(-180, 180, by = 6)
  x_med <- x0 + 3

  ## find the nearest center UTM
  zone <- which.min(abs(x - x_med))

  ## exception
  if (y >= 56 && y < 64 && x >= 0 && x < 6) {
    zone <- 32
  }

  if (y > 0) {

    ## epsg code for northern hemisphere
    if (y < 60) {
      ## non-pole zone
      epsg <- 32600 + zone
    } else {
      ## pole zone
      epsg <- 32661
    }

  } else {

    ## epsg code for southern hemisphere
    if (y >= -60) {
      ## non-pole zone
      epsg <- 32700 + zone
    } else {
      ## pole zone
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
#' @param outlet Outlet point layer of class \code{sf}.
#' @param id_col Column name specifying outlet id.
#'  This column information will be appended to the output polygon layer.
#' @param f_dir Flow direction raster of class \code{SpatRaster}.
#' @param str_grid Stream grid raster of class \code{SpatRaster}.
#' @param snap Logical.
#'  Whether snapping of outlet points should be performed.
#' @param snap_dist Numeric.
#'  Distance threshold for snapping points to stream grid.
#'  Measured in the unit of input raster files.
#' @param export Logical.
#'  Whether output polygons should be exported.
#' @param output_dir Character.
#'  Specify a directly name for output.
#'  Ignored if \code{export = FALSE}.
#' @param filename Character.
#'  Specify a file name of output watershed polygons.
#'  Ignored if \code{export = FALSE}.
#' @param file_ext Character.
#'  Specify a file extension of output files.
#'  Either \code{"gpkg"} or \code{"shp"}.
#'  Ignored if \code{export = FALSE}.
#' @param keep_outlet Logical.
#'  Whether a snapped outlet layer should be exported.
#'  Ignored if \code{export = FALSE}.
#'
#' @importFrom stringr str_detect
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

wsd_unnested <- function(outlet,
                         id_col,
                         f_dir,
                         str_grid = NULL,
                         snap = TRUE,
                         snap_dist = 5,
                         export = FALSE,
                         output_dir = "watershed",
                         filename = "watershed",
                         file_ext = "gpkg",
                         keep_outlet = FALSE) {

  ## temporary files
  message("Saving temporary files...")
  temppath <- tempdir()
  v_name <- paste(temppath,
                  c("strg.tif",
                    "outlet.shp",
                    "outlet_snap.shp",
                    "dir.tif",
                    "wsd.tif"),
                  sep = "\\")

  terra::writeRaster(f_dir,
                     filename = v_name[str_detect(v_name, "dir")],
                     overwrite = TRUE)

  sf::st_write(outlet,
               dsn = v_name[str_detect(v_name, "outlet.shp")],
               append = FALSE)

  ## snapping
  if (snap) {
    ## with snapping
    if (is.null(str_grid))
      stop("Stream grid 'str_grid' must be supplied if snap = TRUE")

    terra::writeRaster(str_grid,
                       filename = v_name[str_detect(v_name, "strg")],
                       overwrite = TRUE)

    message("Snap outlet points to the nearest stream grid...")
    whitebox::wbt_jenson_snap_pour_points(pour_pts = v_name[str_detect(v_name, "outlet\\.")],
                                          streams = v_name[str_detect(v_name, "strg")],
                                          output = v_name[str_detect(v_name, "outlet_snap")],
                                          snap_dist = snap_dist)
  } else {
    ## with no snapping
    sf::st_write(outlet,
                 dsn = v_name[str_detect(v_name, "outlet_snap.shp")],
                 append = FALSE)
  }

  ## delineation
  message("Delineate watersheds...")
  whitebox::wbt_unnest_basins(d8_pntr = v_name[str_detect(v_name, "dir")],
                              pour_pts = v_name[str_detect(v_name, "outlet_snap")],
                              output = v_name[str_detect(v_name, "wsd")])

  ## vectorize
  message("Vectorize raster watersheds...")

  sf_wsd0 <- list.files(path = temppath,
                        pattern = "wsd",
                        full.names = TRUE) %>%
    lapply(terra::rast) %>%
    lapply(stars::st_as_stars) %>%
    lapply(sf::st_as_sf,
           merge = TRUE,
           as_points = FALSE) %>%
    dplyr::bind_rows() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(site_id = sum(dplyr::c_across(cols = dplyr::ends_with("tif")),
                                na.rm = TRUE)) %>%
    dplyr::select(.data$site_id) %>%
    dplyr::ungroup()

  sf_wsd <- sf_wsd0 %>%
    dplyr::mutate(area = units::set_units(sf::st_area(sf_wsd0), "km^2")) %>%
    dplyr::group_by(.data$site_id) %>%
    dplyr::slice(which.max(.data$area)) %>% # remove duplicates by outlet
    dplyr::ungroup() %>%
    dplyr::relocate(.data$site_id, .data$area) %>%
    dplyr::arrange(.data$site_id)

  outlet_snap <- sf::st_read(dsn = v_name[str_detect(v_name, "outlet_snap")]) %>%
    dplyr::select(NULL)

  if (!missing(id_col)) {
    ## pull id_col as an identifier
    identifier <- outlet %>%
      dplyr::pull(id_col)

    ## pull xy coordinates from original and snapped points
    xy0 <- sf::st_coordinates(outlet)
    xy <- sf::st_coordinates(outlet_snap)

    ## append point coordinates to polygons
    site_id_num <- sf_wsd$site_id

    sf_wsd <- sf_wsd %>%
      dplyr::mutate(id_col = identifier[site_id_num],
                    x0 = xy0[site_id_num, 1],
                    y0 = xy0[site_id_num, 2],
                    x = xy[site_id_num, 1],
                    y = xy[site_id_num, 2]) %>%
      dplyr::relocate(.data$id_col,
                      .data$x,
                      .data$y,
                      .data$x0,
                      .data$y0)
  }

  if (export) {

    if (!any(str_detect(list.files(".", recursive = TRUE), output_dir)))
      dir.create(output_dir)

    if (!any(file_ext %in% c("gpkg", "shp")))
      stop("'file_ext' must be either 'gpkg' or 'shp'")

    sf::st_write(sf_wsd,
                 dsn = paste0(output_dir,
                              "/",
                              filename,
                              ".",
                              file_ext),
                 append = FALSE)

    if (keep_outlet) {

      sf::st_write(outlet_snap,
                   dsn = paste0(output_dir,
                                "/",
                                "outlet_snap",
                                ".",
                                file_ext),
                   append = FALSE)
    }

  }

  ## remove temporary files
  message("Removing temporary files...")

  files <- list.files(temppath, full.names = TRUE)
  cl <- call("file.remove", files)
  suppressWarnings(eval(cl, envir = parent.frame()))

  return(sf_wsd)
}

#' Delineate nested watersheds
#'
#' @inheritParams wsd_unnested
#' @param id_col Column name specifying outlet id.
#'  This column information will be appended to the snapped outlet layer.
#' @param simplify Logical.
#'  Whether simplify the output polygons.
#' @param keep Numeric.
#'  Proportion of vertices kept after polygon simplifications.
#'  Ignored if \code{simplify = FALSE}
#'
#' @importFrom stringr str_detect
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

wsd_nested <- function(outlet,
                       id_col,
                       f_dir,
                       str_grid = NULL,
                       snap = TRUE,
                       snap_dist = 5,
                       export = FALSE,
                       output_dir = "watershed",
                       filename = "watershed",
                       file_ext = "gpkg",
                       keep_outlet = FALSE,
                       simplify = FALSE,
                       keep = 0.5) {

  message("Saving temporary files...")
  ## temporary file names
  temppath <- tempdir()
  v_name <- paste(temppath,
                  c("strg.tif",
                    "outlet.shp",
                    "outlet_snap.shp",
                    "dir.tif",
                    "wsd.tif"),
                  sep = "\\")

  ## write temporary files
  terra::writeRaster(f_dir,
                     filename = v_name[str_detect(v_name, "dir")],
                     overwrite = TRUE)

  ## snapping
  if (snap) {
    ## with snapping
    if (is.null(str_grid))
      stop("Stream grid 'str_grid' must be supplied if snap = TRUE")

    terra::writeRaster(str_grid,
                       filename = v_name[str_detect(v_name, "strg")],
                       overwrite = TRUE)

    sf::st_write(outlet,
                 dsn = v_name[str_detect(v_name, "outlet.shp")],
                 append = FALSE)

    message("Snap outlet points to the nearest stream grid...")
    whitebox::wbt_jenson_snap_pour_points(pour_pts = v_name[str_detect(v_name, "outlet\\.")],
                                          streams = v_name[str_detect(v_name, "strg")],
                                          output = v_name[str_detect(v_name, "outlet_snap")],
                                          snap_dist = snap_dist)
  } else {
    ## with no snapping
    sf::st_write(outlet,
                 v_name[str_detect(v_name, "outlet_snap")],
                 append = FALSE)
  }

  ## watershed delineation: raster output
  message("Delineate watersheds...")
  whitebox::wbt_watershed(d8_pntr = v_name[str_detect(v_name, "dir")],
                          pour_pts = v_name[str_detect(v_name, "outlet_snap")],
                          output = v_name[str_detect(v_name, "wsd")])

  ## watershed delineation: vectorize
  message("Vectorize raster watersheds...")
  sf_wsd0 <- terra::rast(v_name[str_detect(v_name, "wsd")]) %>%
    stars::st_as_stars() %>%
    sf::st_as_sf(merge = TRUE,
                 as_point = FALSE) %>%
    dplyr::rename(tifid = .data$wsd.tif)

  sf_wsd <- sf_wsd0 %>%
    dplyr::mutate(area = sf::st_area(sf_wsd0)) %>%
    dplyr::group_by(.data$tifid) %>%
    dplyr::slice(which.max(.data$area)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fid = dplyr::row_number()) %>%
    dplyr::relocate(.data$fid) %>%
    dplyr::select(-area)

  outlet_id <- dplyr::pull(sf_wsd, .data$tifid)

  if (simplify) {
    ## with simplification
    if (!(keep < 1 && keep > 0))
      stop("'keep' must be greater than 0 and less than 1")
    sf_wsd <- rmapshaper::ms_simplify(sf_wsd,
                                      keep = keep) %>%
      sf::st_make_valid()
  } else {
    ## without simplification
    sf_wsd <- sf::st_make_valid(sf_wsd)
  }

  outlet_snap <- sf::st_read(dsn = v_name[str_detect(v_name, "outlet_snap")]) %>%
    dplyr::group_by(.data$geometry) %>%
    dplyr::mutate(pid = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::relocate(.data$pid)

  if (!missing(id_col)) {
    outlet_snap <- outlet_snap %>%
      dplyr::mutate(id_col = dplyr::pull(outlet, id_col)) %>%
      dplyr::relocate(.data$id_col) %>%
      dplyr::slice(outlet_id)
  } else {
    outlet_snap <- outlet_snap %>%
      dplyr::slice(outlet_id)
  }

  ## export
  if (export) {
    ### create export directory
    if (!(output_dir %in% list.files(".")))
      dir.create(output_dir)

    ### watershed polygon
    sf::st_write(sf_wsd,
                 dsn = paste0(output_dir,
                              "/",
                              filename,
                              ".",
                              file_ext),
                 append = FALSE)

    ### snapped outlet
    if (keep_outlet) {
      sf::st_write(outlet_snap,
                   dsn = paste0(output_dir,
                                "/",
                                "outlet_snap",
                                ".",
                                file_ext),
                   append = FALSE)
    }
  }

  ## remove temporary files
  message("Removing temporary files...")
  files <- list.files(temppath, full.names = TRUE)
  cl <- call("file.remove", files)
  suppressWarnings(eval(cl, envir = parent.frame()))

  return(list(watershed = sf_wsd,
              outlet = outlet_snap))
}

#' Convert flow accumulation raster to stream grid
#'
#' @param f_acc Flow accumulation raster of class \code{SpatRaster}.
#' @param threshold Numeric.
#'  Threshold value for minimum stream grid.
#'  The unit inherits from input the flow accumulation layer.
#' @param output Character.
#'  File path for output stream grid.
#'
#' @importFrom stringr str_detect
#' @importFrom dplyr %>%
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

flow2grid <- function(f_acc,
                      threshold,
                      output) {

  ## temporary file names
  temppath <- tempdir()
  fname <- paste0(temppath,
                  "\\",
                  "upa.tif")

  ## write raster input in temporary folder
  terra::writeRaster(f_acc,
                     filename = fname[str_detect(fname, "upa")],
                     overwrite = TRUE)

  ## stream grids
  whitebox::wbt_extract_streams(flow_accum = fname[str_detect(fname, "upa")],
                                output = output,
                                threshold = threshold)

  strg <- terra::rast(output)

  ## remove temporary files
  message("Removing temporary files...")
  files <- list.files(temppath, full.names = TRUE)
  cl <- call("file.remove", files)
  suppressWarnings(eval(cl, envir = parent.frame()))

  return(strg)
}

#' Convert stream grid to vector stream
#'
#' @inheritParams wsd_unnested
#' @param output Character.
#'  File path for output stream vector.
#' @param set_crs Logical.
#'  Whether output file should inherit CRS from the input flow direction raster.
#'
#' @importFrom stringr str_detect
#' @importFrom dplyr %>%
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

grid2stream <- function(f_dir,
                        str_grid,
                        output,
                        set_crs = TRUE) {

  ## temporary file names
  temppath <- tempdir()
  fname <- paste0(temppath,
                  "\\",
                  c("dir.tif", "strg.tif"))

  ## write raster input in temporary folder
  terra::writeRaster(f_dir,
                     filename = fname[str_detect(fname, "dir")],
                     overwrite = TRUE)

  terra::writeRaster(str_grid,
                     filename = fname[str_detect(fname, "strg")],
                     overwrite = TRUE)

  ## stream grids
  whitebox::wbt_raster_streams_to_vector(streams = fname[str_detect(fname, "strg")],
                                         d8_pntr = fname[str_detect(fname, "dir")],
                                         output = output)

  if (set_crs) {
    sf::st_read(output,
                quiet = TRUE) %>%
      sf::st_set_crs(terra::crs(f_dir)) %>%
      sf::st_write(dsn = output,
                   append = FALSE,
                   quiet = TRUE)
  }

  channel <- sf::st_read(output,
                         quiet = TRUE)

  ## remove temporary files
  message("Removing temporary files...")
  files <- list.files(temppath, full.names = TRUE)
  cl <- call("file.remove", files)
  suppressWarnings(eval(cl, envir = parent.frame()))

  return(channel)
}
