#' @export
#' 

print.dsi <- function(object){
  ncon <- length(object$consumers)
  nres <- length(object$resources)
  sing <- sum(object$samp == 1)
  cat(paste("Object of class dsi, with the distance-based specialization index calculated for", 
            ncon, "consumer species, \n using", nres, "resource types \n"))
  print(object$DSIstar)
  if (sing > 0) {
    cat(paste("\n",sing, " consumer species are singletons and specialization was not calculated for them"))
  }
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

#' @export
#' 

print.dsicom <- function(object){
  ncon <- length(object$consumers)
  nres <- length(object$resources)
  ncom <- length(object$communities)
  cat(paste("Object of class dsicom, with the distance-based specialization index calculated for", 
            ncon, " consumer species, \n using", nres, " resource types in", ncom, " communities"))
  cat(paste("\n DSICom values calculated for", ncom, " communities:"))
  print(object$dsicom)
  if (is.data.frame(object$part)) {
    cat("Partition of the variability in DSI* measured locally: \n ")
    print(object$part)
  }
  cat("\n For individual DSI* levels measured locally, access the DSIstar element directly")
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
