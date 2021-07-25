#' @export
#' 

print.dsi <- function(x, ...){
  ncon <- length(x$consumers)
  nres <- length(x$resources)
  sing <- sum(x$samp == 1)
  cat(paste("Object of class dsi, with the distance-based specialization index calculated for", 
            ncon, "consumer species, \n using", nres, "resource types \n"))
  print(x$DSIstar)
  if (sing > 0) {
    cat(paste("\n",sing, " consumer species are singletons and specialization was not calculated for them"))
  }
}

#' @export
#' 

summary.dsi <- function(object, ...){
  Res <- data.frame(N = object$samp, 
                    S = object$richness, 
                    MPD = object$MPD,
                    DSI = object$DSI,
                    DSIstar = object$DSIstar,
                    class = object$class,
                    row.names = object$consumers)
  Res
}

#' @export
#' 

print.dsicom <- function(x, ...){
  ncon <- length(x$consumers)
  nres <- length(x$resources)
  ncom <- length(x$communities)
  cat(paste("Object of class dsicom, with the distance-based specialization index calculated for", 
            ncon, " consumer species, \n using", nres, " resource types in", ncom, " communities"))
  cat(paste("\n DSICom values calculated for", ncom, " communities: \n "))
  print(x$dsicom)
  if (is.data.frame(x$part)) {
    cat("\n Partition of the variability in DSI* measured locally: \n ")
    print(x$part)
  }
  cat("\n For individual DSI* values measured locally, access the DSIstar element directly")
}
