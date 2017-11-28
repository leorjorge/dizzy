
#' Community level Distance-based specialisation index
#'
#'Measures the distance-based specialization indesx (DSI) of consumers in a set of communities. 
#'Relies on resource use information from an interaction matrix or array, 
#'resource similarities from a distance matrix or phylogeny, and local resource availability 
#'data standardizing the null models.
#'
#' @export
#' 
#' @param Int Is an interaction array for which Specialization is to be measured. 
#' Consumers are in rows, resources in the columns and different locations 
#' (or other "local" partition - e.g. time) are slices. More than one locality is necessary to
#' the function to work. With a single location \code{\link{dsi}} should be used instead
#' 
#' @param Dist Object of class phylo with the $tip.label information matching Abundance data. 
#' The phylogeny should have all resource classes present in \code{Abund} 
#' and will be pruned to remove resource classes absent from \code{Abund}. 
#' Alternatively, a distance matrix, with dissimilarities between resource classes. 
#' Dimension names must match \code{Abund}.
#' 
#' @param Abund Matrix of local abundances/sampling effort of the resources in the interaction 
#' array at different locations.
#' 
#' @param Rep Number of iterations for the null model, defaults to 999.
#' 
#' @param Part A logical indicating whether variability in DSI measured locally should be partitioned 
#' between species and localities. If set to \code{TRUE}, (the default), variation in DSI* measured 
#'  locally is orthogonally partitioned between species and communities into three components 
#'  of mean squared distances between DSI* values: differences within species, 
#'  differences within communities, residual variation between species and communities. 
#'  The matrix with DSI* values is then randomized \code{Rep} times, keeping the absences fixed,
#'  and in each iteration the same components of variation are calculated. A Z-score, measuring the
#'  effect size of the observed components relative to the null model is then measured.
#' 
#' @return Returns an object of class \code{dsicom}, with the following elements:
#' \item{data}{A list with the three objects used to calculate DSI - The interaction array, distance
#' matrix and resource availability matrix}
#' \item{consumers}{A character vector with the names of consumer species}
#' \item{resources}{A character vector with the names of resource classes (often species)}
#' \item{communities}{A character vector with the names of communities or other local unities where
#' specialization is measured}
#'  \item{MPD}{A numeric matrix with the abundance averaged mean pairwise distance 
#'  between resource items used by each consumer species in each community separately}
#'  \item{null}{An array with \code{Rep} rows, \code{consumers} in columns, and \code{communities}
#'  in the third dimension, with all \code{MPD} values obtained by the null model}
#'  \item{DSI}{A numeric matrix with unstandardized DSI values measured locally. 
#'  These are calculated as Z-scores of the MPD values, compared to the null distribution 
#'  of MPDs obtained from the null model. It can be used to classify consumers as specialists, 
#'  non-selective or generalists, but is not appropriate for comparisons among species 
#'  with different sampling effort}
#'  \item{lim}{Theoretical maximum (for species with positive DSI) and minimum 
#'  (for species with negative DSI) DSI values attainable for the sample size and 
#'  resource similarity of each consumer species. Used to standardize DSI and generate comparable
#'  DSI* values. Maximum is obtained by assuming the species is a monophage, and minimum is
#'  calculated by using a Simulated annealing algorithm to distribute the recorded number
#'  of interactions among resource classes in such a manner that maximizes MPD}
#'  \item{DSIstar}{A numeric matrix with standardized DSI values (DSI*) measured locally 
#'  in communities. This is the index as proposed in Jorge et al. (2017). 
#'  DSI* values vary between -1 (extreme generalization) to 1 (extreme specialization). 
#'  This should be the preferred value to be used when comparing species and communities,
#'  as it has a very straightforward interpretation and is controlled for differences both in sampling 
#'  intensity, co-occurrence with resources and resource similarity differences.}
#'  \item{dsicom}{A numeric vector with abundance weighted avarage DSI* (DSICom) of the species 
#'  that occur in each community.}
#'  \item{part}{If \code{Part} is set to \code{TRUE} (the default), the three components
#'  of variation in DSI* are presented. \code{IntraSP} amounts to differences within species, 
#'  \code{IntraCom} amounts to differences within communities, and \code{Residual} amounts 
#'  to residual variation between species and communities.}
#'  \item{nullpart}{If \code{Part} is set to \code{TRUE} (the default), the DSI* variation 
#'  components calculated from the null model are recorded in a matrix with components in 
#'  columns and null model iterations in rows.}
#'  \item{partZ}{If \code{Part} is set to \code{TRUE} (the default), the Z-scores with the 
#'  effect sizes of the observed variability components compared to the null model is shown, with
#'  the same nomenclature as \code{part} above}

dsicom <- function(Int, Dist, Abund, Rep=999, Part = TRUE){
  if (dim(Int)[3] == 1) stop("Interaction data for only one community provided")
  if (class(Dist) == "phylo") {
    Phy <- cophenetic(Dist)
  }
  Phy <- as.matrix(Phy)[rownames(Phy) %in% rownames(Abund), colnames(Phy) %in% rownames(Abund)]
  if (length(setdiff(rownames(Abund), rownames(Phy))) > 0) {
    stop("Some resources in Abund absent from Dist")
  }
  if (length(setdiff(colnames(Int), rownames(Phy))) > 0) {
    stop("Some resources in Int absent from Dist")
  }
  Res <- structure(list(), class = "dsicom")
  Res$data$int <- Int
  Res$data$dist <- Phy
  Res$data$abund <- Abund
  Res$consumers <- rownames(Int)
  Res$resources <- colnames(Phy)
  Res$communities <- dimnames(Int)[[3]]
  LocMPD <- apply(X = Int, MARGIN = 3, FUN = mpd2, dis = Phy, abundance.weighted = T)
  rownames(LocMPD) <- Res$consumers
  Res$MPD <- LocMPD
  LocSamp <- apply(Int,c(1,3),sum)
  LocMPD[LocSamp == 0] <- NA
  nullMPDLoc <- array(dim = c(Rep, nrow(LocMPD), ncol(LocMPD)))  
  for (k in 1:dim(nullMPDLoc)[3]) {
    for (j in 1:dim(nullMPDLoc)[2]) {
      if (LocSamp[j,k] <= 1) {
        nullMPDLoc[,j,k] <- NA
      } else {
        for (i in 1:dim(nullMPDLoc)[1]) {
          amostras <- sample(row.names(Abund), prob = Abund[,k], 
                             size = LocSamp[j,k], replace = TRUE)
          nullMPDLoc[i,j,k] <- mpd2(as.matrix(t(table(amostras))), Phy, abundance.weighted = T)
        }
      }
    }
  }
  Res$null <- nullMPDLoc
  nullMPDLoc.mean <- apply(nullMPDLoc,c(2,3),mean, na.rm = T)
  nullMPDLoc.sd <- apply(nullMPDLoc,c(2,3),sd, na.rm = T)
  LocDSI <- -1*(LocMPD - nullMPDLoc.mean)/nullMPDLoc.sd
  Res$DSI <- LocDSI
  LocLim <- matrix(ncol = ncol(LocDSI), nrow = nrow(LocDSI))
  LocDSIPos <- LocDSI >= 0 & !is.na(LocDSI)
  LocDSINeg <- LocDSI < 0 & !is.na(LocDSI)
  LocLim[LocDSIPos] <- -1*(0 - nullMPDLoc.mean[LocDSIPos])/nullMPDLoc.sd[LocDSIPos]#Maximum value of DSI, calculated by assuming species are monophages
  if (sum(LocDSINeg) > 0) {
    GenLoc <- MaxGenLoc(Phy, Int, SpLocs = which(LocDSINeg == T, arr.ind = T))
    LocLim[LocDSINeg] <- -1*(GenLoc - nullMPDLoc.mean[LocDSINeg])/nullMPDLoc.sd[LocDSINeg]#Minimum value of DSI, calculated by using optimized MaxMPD values
  }
  Res$lim <- LocLim
  LocDSI.st <- LocDSI/abs(LocLim)
  Res$DSIstar <- LocDSI.st
  Res$dsicom <- NA
  for (i in 1:ncol(LocDSI.st)) {
    Res$dsicom[i] <- weighted.mean(LocDSI.st[,i], LocSamp[,i], na.rm = T)
  }
  if (Part == T) {
    part <- DSIpart(Res$DSIstar, rep = Rep)
    Res$part <- part$OBS
    Res$nullpart <- part$Null
    Res$partZ <- part$Z
  }
  return(Res)
}

#' Partition of DSI components among locations and species
#' 
#' 
#' @export
#' 
DSImean <- function(DSILoc){
  PairSQDist <- outer(c(DSILoc),c(DSILoc),FUN = function(x,y) (x - y)^2)
  ArrayPairDist <- array(NA,dim = c(nrow(DSILoc),nrow(DSILoc),ncol(DSILoc),ncol(DSILoc)))
  for (i in 1:ncol(DSILoc)) {
    for (j in 1:ncol(DSILoc)) {
      ArrayPairDist[,,i,j] <- PairSQDist[c((i - 1) * nrow(DSILoc) + 1):c(i * nrow(DSILoc)),
                                         c((j - 1) * nrow(DSILoc) + 1):c(j * nrow(DSILoc))]
    }
  }
  IntraComDist <- array(NA,dim = c(nrow(DSILoc), nrow(DSILoc), ncol(DSILoc)))
  InterComDist <- ArrayPairDist
  for (i in 1:ncol(DSILoc)) {
    IntraComDist[,,i] <- ArrayPairDist[,,i,i]
    InterComDist[,,i,i] <- NA
  }
  IntraSPDist <- array(NA,dim = c(nrow(DSILoc),ncol(DSILoc),ncol(DSILoc)))
  for (i in 1:nrow(DSILoc)) {
    IntraSPDist[i,,] <- InterComDist[i,i,,]
    InterComDist[i,i,,] <- NA
  }
  IntraSP <- mean(IntraSPDist, na.rm = T)
  IntraCom <- mean(IntraComDist, na.rm = T)
  InterCom <- mean(InterComDist, na.rm = T)
  return(c(IntraSP = IntraSP, IntraCom = IntraCom, Residual = InterCom))
}

#' Null model to test for the different components calculated through DSImean
#' 

NullNA <- function(DSILoc,rep=999){
  NullDSI <- array(NA,dim = c(dim(DSILoc),rep))
  for (i in 1:rep) {
    NullDSI[,,i][!is.na(DSILoc)] <- sample(DSILoc[!is.na(DSILoc)])
  }
  NullPart <- t(apply(NullDSI,3,DSImean))
  return(NullPart)
}


#' Wrapper to calculate DSImean and the corresponding null model and provide proper results
#' 
#' @export
#' 
DSIpart <- function(Mat, rep=999){
  OBSPart <- DSImean(Mat)
  NPart <- NullNA(Mat,rep)
  NullCI <- rbind(apply(NPart,2,quantile,probs = c(0.025,0.975), na.rm = T))
  Zscore <- (OBSPart - apply(NPart, 2, mean, na.rm = T)) / apply(NPart, 2, sd, na.rm = T)
  Res <- list(Mat = Mat, OBS = OBSPart, Null = NPart, CI = NullCI, Z = Zscore)
  return(Res)
}
