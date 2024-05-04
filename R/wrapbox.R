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

  if (keep >= 1)
    stop("'keep' must be less than 1.0")

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

