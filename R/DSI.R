#' Distance-based specialisation index
#'
#'Measures the distance-based specialization indesx (DSI) of a set of consumers. 
#'Relies on resource use information from an interaction matrix or array, 
#'resource similarities from a distance matrix or phylogeny, and local resource availability 
#'data standardizing the null models.
#'
#' @export
#' 
#' @param Int Is an interaction array for which Specialization is to be measured. 
#' Consumers are in rows, resources in the columns and different locations 
#' (or other "local" partition - e.g. time) are slices. 
#' When the matrix has two dimensions no local variation is assumed.
#' 
#' @param Dist Object of class phylo with the $tip.label information matching Abundance data. 
#' The phylogeny should have all resource species present in the interaction array 
#' and will be pruned accordingly. 
#' Alternatively, a distance matrix, with dissimilarities between resource classes. 
#' Dimension names must match Abundance data.
#' 
#' @param Abund Matrix of local abundances/sampling effort of the resources in the interaction array at different locations
#' 
#' @param Rep Number of iterations for the null model, defaults to 999
#' 

dsi <- function(Int, Dist, Abund, Rep=999){
  if (length(dim(Int)) == 2) {
    Dims <- dimnames(Int)
    Int <- as.matrix(Int)
    dim(Int) <- c(dim(Int),1)
    dimnames(Int) <- Dims
  }
  IntMat <- apply(Int,c(1,2),sum) #Matrix with interactions at the regional level
  if (class(Dist) == "phylo") {
    Phy <- cophenetic(Dist)
  }
  Phy <- as.matrix(Phy)[rownames(Phy) %in% rownames(Abund), colnames(Phy) %in% rownames(Abund)]
  Res <- structure(list(), class = "dsi")
  Res$data$int <- Int
  Res$data$dist <- Phy
  Res$data$abund <- Abund
  Res$consumers <- rownames(Int)
  Res$resources <- colnames(Phy)
  Res$richness <- vegan::specnumber(IntMat) #Number of host plant species
  Res$samp <- apply(Int,1,sum) #Number of samples in which each consumer was collected
  Res$MPD <- mpd2(IntMat,dis = Phy, abundance.weighted = T) #Raw Mean Phylogenetic distance among Regources of a given consumer at the regional level
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
  Res$class[Res$DSI < -1.96] <- "generalist"
  Res$class[Res$DSI > 1.96] <- "specialist"
  Res$class[Res$DSI >= -1.96 && Res$DSI <= 1.96] <- "non-selective"
  DSIPos <- Res$DSI >= 0 & !is.na(Res$DSI)
  DSINeg <- Res$DSI < 0 & !is.na(Res$DSI)
  Res$Lim[DSIPos] <- -1*(0 - Null.MPD.mn[DSIPos])/Null.MPD.sd[DSIPos] #Maximum value of DSI, calculated by assuming all species are monophages
  if (sum(DSINeg) > 0) {
    Gen <- MaxGenReg(Phy, Int, Spps = which(DSINeg == T)) #Calculate the maximum possible MPD by simulated annealing optimization of individuals in resources
    Res$Lim[DSINeg] <- -1*(Gen - Null.MPD.mn[DSINeg])/Null.MPD.sd[DSINeg] #Minimum value of DSI, using the maximum possible value calculated above
  }
  Res$DSIstar <- Res$DSI/abs(Res$Lim) #Final DSI* value, obtained by standardizing DSI with the limit values obtained above
  return(Res)
}


#' Wrapper to calculate DSI for several different guilds
#' 
#' @export
#' 
DSI.Guilds <- function(Int,Dist,Abund,Rep=200){
  Res <- vector("list",length(Int))
  names(Res) <- names(Int)
  for(i in names(Int)){ #This loop calculates results for each interaction matrix separately.
    Res[[i]] <- dsi(Int[[i]],Dist,Abund[[i]],Rep=Rep)
  }
  return(Res)
}
