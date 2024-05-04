#' Return the terminal cell of a given polygon
#'
#' @param data Dataframe.
#'  This object must list the maximum and minimum of coordinates
#'  (latitude or longitude)
#' @param shape \code{sf} polygon object.
#'  This object defines the extent.
#' @param mode Character.
#'  This argument specifies which terminal coordinate of the shape polygon
#'  should be used. Either \code{xmin, xmax, ymin, ymax}
#'
#' @export

get_tf <- function(data, shape, mode) {

  choice <- c("xmin", "ymin", "xmax", "ymax")

  if (!(any(choice == mode)))
    stop("Invalid mode")

  if (!(any(colnames(data) == "min")) || !(any(colnames(data) == "max")))
    stop("Invalid column names")

  cout <- sapply(seq_len(nrow(data)),
                 function(i) {
                   dplyr::between(sf::st_bbox(shape)[mode],
                                  min(data$min[i], data$max[i]),
                                  max(data$min[i], data$max[i]))
                 })

  return(cout)
}
