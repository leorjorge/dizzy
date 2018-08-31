#' Distance-based specialisation index
#'
#' \code{dsi} measures the distance-based specialization indesx (DSI*) of a set of consumers. 
#' Relies on resource use information from an interaction matrix or array, 
#' resource similarities from a distance matrix or phylogeny, and local resource availability 
#' data for creating diet null models.
#'
#' @export
#' 
#' @param Int Is an interaction array for which Specialization is to be measured. 
#' Consumers are in rows, resources in the columns and different locations 
#' (or other "local" partition - e.g. time) are slices. 
#' When the matrix has two dimensions no local variation is assumed.
#' 
#' @param Dist Object of class phylo with the $tip.label information matching Abundance data. 
#' The phylogeny should have all resource classes present in \code{Abund} 
#' and will be pruned to remove resource classes absent from \code{Abund}. 
#' Alternatively, a distance matrix, with dissimilarities between resource classes. 
#' Dimension names must match \code{Abund}.
#' 
#' @param Abund Matrix of local abundances/sampling effort of the resources in the 
#' interaction array at different locations. Resource items are in rows and locations in columns,
#' with labels as row and column names. Abundances can be inputted as relative or absolute values.
#' 
#' @param Rep Number of iterations for the null model, defaults to 999
#' 
#' @return Returns an object of class \code{dsi}, with the following elements:
#' \item{data}{A list with the three objects used to calculate DSI - The interaction array, distance
#' matrix and resource availability matrix}
#' \item{consumers}{A character vector with the names of consumer species}
#' \item{resources}{A character vector with the names of resource classes (often species)}
#' \item{richness}{A numeric vector with the number of resource classes used by each consumer species}
#' \item{samp}{A numeric vector with the number of individuals or interactions recorded for each
#'  consumer species}
#'  \item{MPD}{A numeric vector with the abundance averaged mean pairwise distance 
#'  between resource items used by each consumer species}
#'  \item{null}{A matrix with \code{Rep} rows and \code{consumers} in columns, 
#'  with all \code{MPD} values obtained by the null model}
#'  \item{DSI}{A numeric vector with unstandardized DSI values. These are calculated as Z-scores 
#'  of the MPD values, compared to the null distribution of MPDs obtained from the null model. 
#'  It can be used to classify consumers as specialists, non-selective or generalists, 
#'  but is not appropriate for comparisons among species with different sampling effort}
#'  \item{class}{A character vector with the specialization class of consumers. 
#'  Singletons will return NA}
#'  \item{lim}{Theoretical maximum (for species with positive DSI) and minimum 
#'  (for species with negative DSI) DSI values attainable for the sample size and 
#'  resource similarity of each consumer species. Used to standardize DSI and generate comparable
#'  DSI* values. Maximum is obtained by assuming the species is a monophage, and minimum is
#'  calculated by using a Simulated annealing algorithm to distribute the recorded number
#'  of interactions among resource classes in such a manner that maximizes MPD}
#'  \item{DSIstar}{A numeric vector with standardized DSI values (DSI*). This is the index as 
#'  proposed in Jorge et al. (2017). DSI* values vary between -1 (extreme generalization) to 1 
#'  (extreme specialization). This should be the preferred value to be used when comparing species,
#'  as it has a very straightforward interpretation and is controlled for differences both in sampling 
#'  intensity, co-occurrence with resources and resource similarity differences.}

dsi <- function(Int, Dist, Abund, Rep=999){
  if (length(dim(Int)) == 2) {
    Dims <- dimnames(Int)
    Int <- as.matrix(Int)
    dim(Int) <- c(dim(Int),1)
    dimnames(Int) <- Dims
  }
  IntMat <- apply(Int,c(1,2),sum) #Matrix with interactions at the regional level
  if (class(Dist) == "phylo") {
    Dist <- cophenetic(Dist)
  }
  Dist <- as.matrix(Dist)
  Phy <- Dist[rownames(Dist) %in% rownames(Abund), colnames(Dist) %in% rownames(Abund)]
  if (length(setdiff(rownames(Abund), rownames(Phy))) > 0) {
    stop("One or more resources in Abund absent from Dist")
  }
  if (length(setdiff(colnames(Int), rownames(Phy))) > 0) {
    stop("One or more resources in Int absent from Dist")
  }
  Res <- structure(list(), class = "dsi")
  Res$data$int <- Int
  Res$data$dist <- Phy
  Res$data$abund <- Abund
  Res$consumers <- rownames(Int)
  Res$resources <- colnames(Phy)
  Res$richness <- vegan::specnumber(IntMat) #Number of host plant species
  names(Res$richness) <- Res$consumers
  Res$samp <- apply(Int,1,sum) #Number of samples in which each consumer was collected
  names(Res$samp) <- Res$consumers
  Res$MPD <- mpd2(IntMat,dis = Phy, abundance.weighted = T) #Raw Mean Phylogenetic distance among Regources of a given consumer at the regional level
  names(Res$MPD) <- Res$consumers
  Null.MPD <- matrix(NA,Rep,length(Res$consumers)) #Null model. Each cell is MPD of a given model iteration (rows) for one species (columns).
  LocSamp <- apply(Int,c(1,3),sum)
  for (j in 1:ncol(Null.MPD)) {
    if (sum(LocSamp[j,]) <= 1) {
      Null.MPD[,j] <- NA
    } else {
      if (ncol(LocSamp) > 1) {
        Pesos <- apply(sweep(Abund, MARGIN = 2, LocSamp[j,], "*"), 1, FUN = sum)
        for (i in 1:nrow(Null.MPD)) {
          amostras <- sample(rownames(Abund), prob = Pesos, size = Res$samp[j], replace = TRUE)
          Null.MPD[i,j] <- mpd2(as.matrix(t(table(amostras))), Phy, abundance.weighted = T)
        } 
      } else {
        for (i in 1:nrow(Null.MPD)) {
          amostras <- sample(rownames(Abund), prob =  Abund[,1],
                             size = Res$samp[j], replace = TRUE)
          Null.MPD[i,j] <- mpd2(as.matrix(t(table(amostras))), Phy, abundance.weighted = T)
        }
      }
    }
  }
  Res$null <- Null.MPD
  Null.MPD.mn <- apply(Null.MPD,2,mean, na.rm = T) #Mean MPD for all null model iterations for each species
  Null.MPD.sd <- apply(Null.MPD,2,sd, na.rm = T) #Standard deviation of MPD for all null model iterations
  Res$DSI <- -1*(Res$MPD - Null.MPD.mn)/Null.MPD.sd #DSI Result measured as a Z-score of MPD
  names(Res$DSI) <- Res$consumers
  Res$class[Res$DSI < -1.96] <- "generalist"
  Res$class[Res$DSI > 1.96] <- "specialist"
  Res$class[Res$DSI >= -1.96 & Res$DSI <= 1.96] <- "non-selective"
  names(Res$class) <- Res$consumers
  DSIPos <- Res$DSI >= 0 & !is.na(Res$DSI)
  DSINeg <- Res$DSI < 0 & !is.na(Res$DSI)
  Res$lim[DSIPos] <- -1*(0 - Null.MPD.mn[DSIPos])/Null.MPD.sd[DSIPos] #Maximum value of DSI, calculated by assuming all species are monophages
  if (sum(DSINeg) > 0) {
    Gen <- MaxGenReg(Phy, Int, Spps = which(DSINeg == T)) #Calculate the maximum possible MPD by simulated annealing optimization of individuals in resources
    Res$lim[DSINeg] <- -1*(Gen - Null.MPD.mn[DSINeg])/Null.MPD.sd[DSINeg] #Minimum value of DSI, using the maximum possible value calculated above
  }
  names(Res$lim) <- Res$consumers
  Res$DSIstar <- Res$DSI/abs(Res$lim) #Final DSI* value, obtained by standardizing DSI with the limit values obtained above
  names(Res$DSIstar) <- Res$consumers
  return(Res)
}

