#' Return file name keys for MERITHYDRO rasters
#'
#' @inheritParams get_tf
#'
#' @section MERITHYDRO: https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_Hydro/
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
