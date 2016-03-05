#' Plotting the partiion
#' 
#' @export
#' 
#' 
Plot.DSIpart <- function(Part){
  x <- factor(names(Part$OBS),levels(factor(names(Part$OBS)))[c(3,2,1)])
  plot.default(x[-4],Part$OBS[-4], pch=16, ylim=range(c(Part$OBS,Part$CI)), xaxt="n", ylab="Mean Squares", xlab="Component")
  axis(1, 1:3, x[-4])
  segments(x0=1:3,y0=Part$CI[1,-4], y1=Part$CI[2,-4])
}
