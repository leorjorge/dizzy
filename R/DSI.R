###DSI Calculation####

#This function should be able to measure DSI from an Interaction matrix and a Phylogeny. Sampling effort data is necessary for the null model.

dsi <- function(Int,Phylo,Abund,Rep=200){
  #Int is an interaction array for which Specialization is to be measured. Consumers are in rows, resources in the columns and different locations (or other "local" partition - e.g. time) are slices. When the matrix has two dimensions no local variation is assumed.
  #Phylo is a phylo object with the $tip.label information matching column names. The phylogeny should have all resource species present in the interaction array and will be pruned accordingly
  #Abund is a matrix of local abundances/sampling effort of the resources in the interaction array at different locations
  #Rep is the number of iterations for the null model selected
  if (length(dim(Int))==2){
    Dims <- dimnames(Int)
    Int <- as.matrix(Int)
    dim(Int) <- c(dim(Int),1)
    dimnames(Int) <- Dims
  }
  IntMat <- apply(Int,c(1,2),sum) #Matrix with interactions at the regional level
  Phy <- prune.sample(IntMat,Phylo) #Prune phylogeny for the resource species in ith interaction matrix
  print("Regional data")
  Reg <- data.frame(row.names = rownames(Int)) #Data frame to store the results. 
  Reg$Rich <- specnumber(IntMat) #Number of host plant species
  Reg$Div <- diversity(IntMat,"simpson") #Simpson diversity of host plant species
  #  Reg$dPrime <- specieslevel(t(IntMat),index="d", level="higher")[,1] #Blüthgen's d Prime calculated from the interaction matrix
  Reg$Samp <- apply(Int,1,sum) #Number of samples in which each consumer was collected
  Reg$MPD <- mpd2(IntMat,dis=cophenetic(Phy),abundance.weighted=T) #Raw Mean Phylogenetic distance among Regources of a given consumer at the regional level
  Gen <- MaxGenReg(Phy, Int) #Regional Max Phylo Dispersion calculation. MPD is calculated
  Null.MPD <- matrix(NA,Rep,nrow(Reg)) #Null model. Should be re-written as a separate function, so that for different datasets (geographic variation, no sampling effort etc.) different null models are called. Each cell is MPD of a given model iteration (rows) for one species (columns).
  LocSamp <- apply(Int,c(1,3),sum)
  for(j in 1:ncol(Null.MPD)){
    print(j)
    if (ncol(LocSamp)>1){
      Locs <- dimnames(Int)[[3]][LocSamp[j,]>0]
      Avail <- rep.int(x = rownames(Abund),times = as.numeric(apply(as.matrix(Abund[,Locs]),1,sum)))
      Pesos <- rep(x=LocSamp[j,LocSamp[j,]>0], times=colSums(Abund)[Locs])
      for (i in 1:nrow(Null.MPD)){
        amostras <- sample(Avail, prob=Pesos, size=Reg$Samp[j], replace=FALSE)
        Null.MPD[i,j] <- mpd2(as.matrix(t(table(amostras))), cophenetic(Phy),abundance.weighted=T)
      } 
    } else {
      Avail <- rep.int(x=rownames(Abund),times=Abund[,1])
      for (i in 1:nrow(Null.MPD)){
        amostras <- sample(Avail, size=Reg$Samp[j], replace=TRUE)
        Null.MPD[i,j] <- mpd2(as.matrix(t(table(amostras))), cophenetic(Phy),abundance.weighted=T)
      }
    }
  }
  Null.MPD.mn <- apply(Null.MPD,2,mean, na.rm=T) #Mean MPD for all null model iterations for each species
  Null.MPD.sd <- apply(Null.MPD,2,sd, na.rm=T) #Standard deviation of MPD for all null model iterations
  Reg$Max <- -1*(0-Null.MPD.mn)/Null.MPD.sd #Maximum value of DSI, calculated by assuming all species are monophages
  Reg$Min <- (Gen-Null.MPD.mn)/Null.MPD.sd #Minimum value of DSI. Distributes individuals in the phylogeny according to uniqueness.
  Reg$DSI <- -1*(Reg$MPD-Null.MPD.mn)/Null.MPD.sd #DSI Regult measured as a Z-score of MPD
  DSIPos <- Reg$DSI>=0 & !is.na(Reg$DSI)
  DSINeg <- Reg$DSI<0 & !is.na(Reg$DSI)
  Reg$DSI.st[DSIPos] <- Reg$DSI[DSIPos]/Reg$Max[DSIPos]
  Reg$DSI.st[DSINeg] <- Reg$DSI[DSINeg]/Reg$Min[DSINeg]#DSI Regcaled to vary between Max and Min.
  if (dim(Int)[3]>1){
    #Local data 
    print("Local Data")
    LocMPD <- apply(X=Int, MARGIN=3, FUN=mpd2, dis=cophenetic(Phy), abundance.weighted=T)
    row.names(LocMPD) <- row.names(Int)
    GenLoc <- MaxGenLoc(Phy, Int)    
    nullMPDLoc <- array(dim=c(Rep,nrow(LocMPD),ncol(LocMPD)))  
    for(k in 1:dim(nullMPDLoc)[3]){
      print(k)
      for(j in 1:dim(nullMPDLoc)[2]){
        Avail <- rep(x=row.names(Abund),times=Abund[,k])
        for(i in 1:dim(nullMPDLoc)[1]){
          amostras <- sample(Avail, size=apply(Int,c(1,3),sum)[j,k])
          nullMPDLoc[i,j,k] <- mpd2(as.matrix(t(table(amostras))), cophenetic(Phy),abundance.weighted=T)
        }
      }
    }
    nullMPDLoc.mean <- apply(nullMPDLoc,c(2,3),mean, na.rm=T)
    nullMPDLoc.sd <- apply(nullMPDLoc,c(2,3),sd, na.rm=T)
    LocDSI <- -1*(LocMPD-nullMPDLoc.mean)/nullMPDLoc.sd
    LocMax <- -1*(0-nullMPDLoc.mean)/nullMPDLoc.sd
    LocMin <- -1*(GenLoc-nullMPDLoc.mean)/nullMPDLoc.sd
    LocDSIPos <- LocDSI>=0 & !is.na(LocDSI)
    LocDSINeg <- LocDSI<0 & !is.na(LocDSI)
    LocDSI.st <- array(dim = dim(LocDSI))
    LocDSI.st[LocDSIPos] <- LocDSI[LocDSIPos]/LocMax[LocDSIPos]
    LocDSI.st[LocDSINeg] <- LocDSI[LocDSINeg]/LocMin[LocDSINeg]
    Loc <- list(MPD=LocMPD, DSI=LocDSI, Max=LocMax, Min=LocMin, DSI.st=LocDSI.st)
    Loc$DSICOM <- NA
    for (i in 1:ncol(LocDSI.st)){
      Loc$DSICOM[i] <- weighted.mean(LocDSI.st[,i],LocSamp[,i], na.rm=T)
    }
    return (list(Reg=Reg,Loc=Loc))
  } else {
    return (Reg)
  }
}


####Wrapper to calculate DSI for several different guilds####
DSI.Guilds <- function(Int,Phylo,Abund,Rep=200){
  Res <- vector("list",length(Int))
  names(Res) <- names(Int)
  for(i in names(Int)){ #This loop calculates results for each interaction matrix separately.
    Res[[i]] <- dsi(Int[[i]],Phylo,Abund[[i]],Rep=Rep)
  }
  return(Res)
}
