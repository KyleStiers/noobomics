#' find_plot_range
#'
#' This function finds the lower and upper limits from two seperate vectors to find the right symmetrical boundaries for plot axes 
#'
#' @param d1 First data frame column 
#' @param d2 Second data frame column
#' @export
#' @examples
#' d1.x <- c(0:10)
#' d2.x <- c(sin(d1.x))
#' x_range <- find_plot_range(df1$var1, df2$var2)

find_plot_range <- function(d1, d2){
  low = 0
  hi = 0
  low <- min(c(range(d1)[1], range(d2)[1]))
  hi <- max(c(range(d1)[2], range(d2)[2]))
  return (c(low,hi))
}