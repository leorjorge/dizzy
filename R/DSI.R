#' Distance-based specialisation index
#'
#'Measures DSI from an Interaction matrix and a Phylogeny. Sampling effort data is necessary for the null model.
#'
#' @export
#' 
#' @param Int Is an interaction array for which Specialization is to be measured. Consumers are in rows, resources in the columns and different locations (or other "local" partition - e.g. time) are slices. When the matrix has two dimensions no local variation is assumed.
#' 
#' @param Dist Object of class phylo with the $tip.label information matching Abundance data. The phylogeny should have all resource species present in the interaction array and will be pruned accordingly. Alternatively, a distance matrix, with dissimilarities between resource classes. Dimension names must match Abundance data.
#' 
#' @param Abund Matrix of local abundances/sampling effort of the resources in the interaction array at different locations
#' 
#' @param Rep Number of iterations for the null model selected
#' 
#' @param DSICom Whether specialization is to be calculated locally.
#' 
dsi <- function(Int,Phylo,Abund,Rep=200, DSICom = T){
  if (length(dim(Int)) == 2) {
    Dims <- dimnames(Int)
    Int <- as.matrix(Int)
    dim(Int) <- c(dim(Int),1)
    dimnames(Int) <- Dims
  }
  IntMat <- apply(Int,c(1,2),sum) #Matrix with interactions at the regional level
  if (class(Phylo) == "phylo") {
    Phy <- cophenetic(Phylo)
  }
  Phy <- as.matrix(Phy)[rownames(Phy) %in% rownames(Abund), colnames(Phy) %in% rownames(Abund)]
  print("Calculating DSI* at the regional level")
  Reg <- data.frame(row.names = rownames(Int)) #Data frame to store the results. 
  Reg$Rich <- vegan::specnumber(IntMat) #Number of host plant species
  Reg$Samp <- apply(Int,1,sum) #Number of samples in which each consumer was collected
  Reg$MPD <- mpd2(IntMat,dis = Phy, abundance.weighted = T) #Raw Mean Phylogenetic distance among Regources of a given consumer at the regional level
  Null.MPD <- matrix(NA,Rep,nrow(Reg)) #Null model. Each cell is MPD of a given model iteration (rows) for one species (columns).
  LocSamp <- apply(Int,c(1,3),sum)
  for (j in 1:ncol(Null.MPD)) {
    if (sum(LocSamp[j,]) <= 1) {
      Null.MPD[,j] <- NA
    } else {
      if (ncol(LocSamp) > 1) {
        Pesos <- apply(sweep(Abund, MARGIN = 2, LocSamp[j,], "*"), 1, FUN = sum)
        for (i in 1:nrow(Null.MPD)) {
          amostras <- sample(rownames(Abund), prob = Pesos, size = Reg$Samp[j], replace = TRUE)
          Null.MPD[i,j] <- mpd2(as.matrix(t(table(amostras))), Phy, abundance.weighted = T)
        } 
      } else {
        for (i in 1:nrow(Null.MPD)) {
          amostras <- sample(rownames(Abund), prob =  Abund[,1],
                             size = Reg$Samp[j], replace = TRUE)
          Null.MPD[i,j] <- mpd2(as.matrix(t(table(amostras))), Phy, abundance.weighted = T)
        }
      }
    }
  }
  Null.MPD.mn <- apply(Null.MPD,2,mean, na.rm = T) #Mean MPD for all null model iterations for each species
  Null.MPD.sd <- apply(Null.MPD,2,sd, na.rm = T) #Standard deviation of MPD for all null model iterations
  Reg$DSI <- -1*(Reg$MPD - Null.MPD.mn)/Null.MPD.sd #DSI Result measured as a Z-score of MPD
  DSIPos <- Reg$DSI >= 0 & !is.na(Reg$DSI)
  DSINeg <- Reg$DSI < 0 & !is.na(Reg$DSI)
  Reg$Lim[DSIPos] <- -1*(0 - Null.MPD.mn[DSIPos])/Null.MPD.sd[DSIPos] #Maximum value of DSI, calculated by assuming all species are monophages
  if (sum(DSINeg) > 0) {
    Gen <- MaxGenReg(Phy, Int, Spps = which(DSINeg == T)) #Calculate the maximum possible MPD by simulated annealing optimization of individuals in resources
    Reg$Lim[DSINeg] <- -1*(Gen - Null.MPD.mn[DSINeg])/Null.MPD.sd[DSINeg] #Minimum value of DSI, using the maximum possible value calculated above
  }
  Reg$DSI.st <- Reg$DSI/abs(Reg$Lim) #Final DSI* value, obtained by standardizing DSI with the limit values obtained above
  if (dim(Int)[3] > 1 & DSICom == T) {
    #Local data 
    print("Calculating DSI* at the local level")
    LocMPD <- apply(X = Int, MARGIN = 3, FUN = mpd2, dis = Phy, abundance.weighted = T)
    row.names(LocMPD) <- row.names(Int)
    LocMPD[LocSamp == 0] <- NA
    nullMPDLoc <- array(dim = c(Rep, nrow(LocMPD), ncol(LocMPD)))  
    for (k in 1:dim(nullMPDLoc)[3]) {
      for (j in 1:dim(nullMPDLoc)[2]) {
        if (LocSamp[j,k] <= 1) {
          nullMPDLoc[,j,k] <- NA
        } else {
          Avail <- rep(x = row.names(Abund), times = Abund[,k])
          for (i in 1:dim(nullMPDLoc)[1]) {
            amostras <- sample(row.names(Abund), prob = Abund[,k], 
                               size = LocSamp[j,k], replace = TRUE)
            nullMPDLoc[i,j,k] <- mpd2(as.matrix(t(table(amostras))), Phy, abundance.weighted = T)
          }
        }
      }
    }
    nullMPDLoc.mean <- apply(nullMPDLoc,c(2,3),mean, na.rm = T)
    nullMPDLoc.sd <- apply(nullMPDLoc,c(2,3),sd, na.rm = T)
    LocDSI <- -1*(LocMPD - nullMPDLoc.mean)/nullMPDLoc.sd
    LocLim <- matrix(ncol = ncol(LocDSI), nrow = nrow(LocDSI))
    LocDSIPos <- LocDSI >= 0 & !is.na(LocDSI)
    LocDSINeg <- LocDSI < 0 & !is.na(LocDSI)
    LocLim[LocDSIPos] <- -1*(0 - nullMPDLoc.mean[LocDSIPos])/nullMPDLoc.sd[LocDSIPos]#Maximum value of DSI, calculated by assuming species are monophages
    if (sum(LocDSINeg) > 0) {
      GenLoc <- MaxGenLoc(Phy, Int, SpLocs = which(LocDSINeg == T, arr.ind = T))
      LocLim[LocDSINeg] <- -1*(GenLoc - nullMPDLoc.mean[LocDSINeg])/nullMPDLoc.sd[LocDSINeg]#Minimum value of DSI, calculated by using optimized MaxMPD values
    }
    LocDSI.st <- LocDSI/abs(LocLim)
    Loc <- list(MPD = LocMPD, DSI = LocDSI, Lim = LocLim, DSI.st = LocDSI.st)
    Loc$DSICOM <- NA
    for (i in 1:ncol(LocDSI.st)) {
      Loc$DSICOM[i] <- weighted.mean(LocDSI.st[,i], LocSamp[,i], na.rm = T)
    }
    return(list(Reg = Reg,Loc = Loc))
  } else {
    return(Reg)
  }
}


#' Wrapper to calculate DSI for several different guilds
#' 
#' @export
#' 
DSI.Guilds <- function(Int,Phylo,Abund,Rep=200){
  Res <- vector("list",length(Int))
  names(Res) <- names(Int)
  for(i in names(Int)){ #This loop calculates results for each interaction matrix separately.
    Res[[i]] <- dsi(Int[[i]],Phylo,Abund[[i]],Rep=Rep)
  }
  return(Res)
}
