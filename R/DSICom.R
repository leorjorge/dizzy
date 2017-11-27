
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
#' (or other "local" partition - e.g. time) are slices. 
#' 
#' @param Dist Object of class phylo with the $tip.label information matching Abundance data. 
#' The phylogeny should have all resource species present in the Abundance data matrix 
#' and will be pruned accordingly. 
#' Alternatively, a distance matrix, with dissimilarities between resource classes. 
#' Dimension names must match Abundance data.
#' 
#' @param Abund Matrix of local abundances/sampling effort of the resources in the interaction array at different locations
#' 
#' @param Rep Number of iterations for the null model, defaults to 999
#' 
#' @param Part a logical indicating whether variability in DSI measured locally should be partitioned 
#' between species and localities.

dsicom <- function(Int, Dist, Abund, Rep=999, Part = TRUE){
  if (dim(Int)[3] == 1) stop("Interaction data for only one community provided")
  if (class(Dist) == "phylo") {
    Phy <- cophenetic(Dist)
  }
  Phy <- as.matrix(Phy)[rownames(Phy) %in% rownames(Abund), colnames(Phy) %in% rownames(Abund)]
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
    part <- DSIpart(Res$DSIstar)
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

NullNA <- function(DSILoc,rep=1000){
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
