#' Plot Frequency Histogram
#'
#' This function will calculate the frequency and plot the histogram
#'
#' @param data a numeric vector
#' @param break_by the range of each bar
#' @param barnum number of bars in total
#' @param main title
#' @param xlab x-axis label
#' @export
#' @details
#' This is build on hist(), with a easier way to set the breaks.
#' @examples
#'

function(data, break_by=1, barnum = 0,main="", xlab=""){
  if (barnum > 0) break_by <- (max(data)-min(data))/barnum
  hist(data, breaks = seq(min(data), max(data)+break_by, by=break_by), main = main, xlab = xlab)
}
