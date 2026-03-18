#' Utility function: Return the terminal cell of a given polygon
#'
#' @param data Dataframe with columns \code{min} and \code{max} listing
#'  coordinate ranges (latitude or longitude).
#' @param shape Polygon object of class \code{sf} defining the extent.
#' @param mode Character. Which bounding box coordinate of \code{shape}
#'  to use. One of \code{"xmin"}, \code{"xmax"}, \code{"ymin"}, \code{"ymax"}.
#'
#' @return A logical vector of length \code{nrow(data)}, indicating whether
#'  the bounding box coordinate falls within each row's \code{min}/\code{max}
#'  range.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

get_tf <- function(data,
                   shape,
                   mode) {

  choice <- c("xmin", "ymin", "xmax", "ymax")

  mode <- match.arg(mode, choice)

  if (!all(c("min", "max") %in% colnames(data)))
    stop("'data' must have columns 'min' and 'max'")

  if (!inherits(shape, "sf"))
    stop("'shape' must be an sf object")

  bbox_val <- sf::st_bbox(shape)[mode]

  cout <- sapply(seq_len(nrow(data)), function(i) {
    dplyr::between(bbox_val,
                   min(data$min[i], data$max[i]),
                   max(data$min[i], data$max[i]))
  })

  return(cout)
}
