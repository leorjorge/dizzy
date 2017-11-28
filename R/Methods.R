#' @export
#' 

print.dsi <- function(object){
  ncon <- length(object$consumers)
  nres <- length(object$resources)
  sing <- sum(object$samp == 1)
  cat(paste("\nObject of class dsi, with the distance-based specialization index calculated for", 
            ncon, " consumer species, using", nres, " resource types"))
  if (sing > 0) {
    cat(paste("\n",sing, " consumers are singletons and specialization was not calculated for them"))
  }
  cat("\n For DSI* values, use summary or access the object directly")
}

#' @export
#' 

summary.dsi <- function(object){
  Res <- data.frame(N = object$samp, 
                    S = object$richness, 
                    MPD = object$MPD,
                    DSI = object$DSI,
                    DSIstar = object$DSIstar,
                    class = object$class,
                    row.names = object$consumers)
  print(Res)
}


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
