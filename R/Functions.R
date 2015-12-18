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

#### MaxGen Functions, used to calculate the maximum possible MPD for a given number of individuals in the phylogeny.####
#The function uses as input a phylo object and the interaction data. It should be used within the DSI calculation function and measures the maximum values of MPDs for the consumers in a given Interaction array
#It returns the Maximum MPD for each herbivore species
require(adephylo, quietly=TRUE)
MaxGenReg <- function (Phy, Int){
  Gen.Mat <- matrix(nrow=nrow(Int), ncol=500)
  Samp <- apply(Int, 1, sum)
  ConsLoc <- apply(Int, c(1,3), sum)
  ConsLoc <- ConsLoc>0
  ResLoc <- apply(Int,c(2,3),sum)
  for(i in 1:nrow(Gen.Mat)) {
    if (Samp[i]<2){
      Gen.Mat[i,] <- NA
    } else {
    ResSums <- rowSums(as.matrix(ResLoc[,ConsLoc[i,]]))
    Phy.sp <- prune.sample(t(ResSums[ResSums>0]), Phy)
    Prob <- 1/(2^distRoot(Phy.sp,method="nNode"))
    #     if (min(Samp[i]*Prob)>1){
    #       Gen.Mat[i,] <- mpd2(matrix(Samp[i]*Prob,nrow=1,dimnames=list(NULL,Phy$tip.label)),
    #                           dis=cophenetic(Phy),abundance.weighted = T)
    #     } else {
    for (j in 1:ncol(Gen.Mat)){
      SPs <- factor(sample(Phy.sp$tip.label,prob=Prob,size = Samp[i],replace = T),levels=Phy.sp$tip.label)
      Gen.Mat[i,j] <- mpd2(t(table(SPs)),dis=cophenetic(Phy.sp),abundance.weighted = T)
          }
    }
  }
  MaxGenReg <- apply(Gen.Mat,1,max)
  return (MaxGenReg)
}

MaxGenLoc <- function (Phy, Int){
  Gen.Arr <- array(dim=c(dim(Int)[c(1,3)],500))
  for(k in 1:dim(Gen.Arr)[2]){
    IntLoc <- Int[,colSums(Int[,,k])>0,k]
    Samp <- rowSums(IntLoc)
    PhyLoc <- prune.sample(IntLoc,Phy)
    Prob <- 1/(2^distRoot(PhyLoc,method="nNode"))
    for(i in 1:dim(Gen.Arr)[1]) {
      if (Samp[i]<2){
        Gen.Arr[i,k,] <- NA
      } else {
        for (j in 1:dim(Gen.Arr)[3]){
          SPs <- factor(sample(PhyLoc$tip.label,prob=Prob,size = Samp[i],replace = T),levels=Phy$tip.label)
          Gen.Arr[i,k,j] <- mpd2(t(table(SPs)),dis=cophenetic(PhyLoc),abundance.weighted = T)
        }
      }
    }
  }
  MaxGen <- apply(Gen.Arr,c(1,2),max)
  return (MaxGen)
}

####Mean DSI.st for each guild and plant species####
GuildPart <- function(Int,DSI, weighted=TRUE) {
  Plant.list <- levels(as.factor(unlist(lapply(Int,colnames))))
  DSIsts <- Int
  DSIstPlant <- as.data.frame(matrix(nrow = length(Int), ncol=length(Plant.list),dimnames = list(names(Int),Plant.list)))
  for (i in names(Int)){
    for(j in 1:nrow(Int[[i]])){
      DSIsts[[i]][j,DSIsts[[i]][j,]>0] <- DSI[[i]]$DSI.st[j]
      DSIsts[[i]][j,DSIsts[[i]][j,]==0] <- NA
    }
    if (weighted==TRUE){
      DSIstsmean <- numeric(ncol(Int[[i]]))
      for (k in 1:length(DSIstsmean)){
        DSIstsmean[k] <- weighted.mean(DSIsts[[i]][,k], weighted=Int[[i]][,k][Int[[i]][,k]>0], na.rm=T)
      }
    } else {
      DSIstsmean <- as.data.frame(t(apply(DSIsts[[i]],2,mean, na.rm=T)))
    }
    DSIstPlant[i,colnames(DSIstPlant) %in% colnames(Int[[i]])] <- DSIstsmean
  }
  return(as.matrix(DSIstPlant))
}

#####MPD####
#Function to calculate MPD correctly. The rationale is similar to the MPD function in picante, but the diagonal in the sample weights matrix is corrected, using the number of combinations of diferent individuals in the same species at the diagonal. MPD for 1 single species is also changed to 0 instead of NA.

mpd2 <- function (samp, dis, abundance.weighted = FALSE) {
  N <- dim(samp)[1]
  mpd <- numeric(N)
  for (i in 1:N) {
    sppInSample <- names(samp[i, samp[i, ] > 0])
    if (length(sppInSample) > 1) {
      sample.dis <- dis[sppInSample, sppInSample]
      if (abundance.weighted) {
        sample.weights <- t(as.matrix(samp[i, sppInSample, 
                                           drop = FALSE])) %*% as.matrix(samp[i, sppInSample, 
                                                                              drop = FALSE])
        diag(sample.weights) <- as.numeric(samp[i,sppInSample]*(samp[i,sppInSample]-1)/2)
        mpd[i] <- weighted.mean(sample.dis[lower.tri(sample.dis, diag=TRUE)], sample.weights[lower.tri(sample.weights, diag=TRUE)])
      }
      else {
        mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
      }
    }
    else {
      mpd[i] <- 0
    }
  }
  mpd
}

####Wrapper to calculate DSI for several different guilds####
DSI.Guilds <- function(Int,Phylo,Abund,Rep=200){
  require(picante,quietly = T)
  require(bipartite,quietly = T)
  Res <- vector("list",length(Int))
  names(Res) <- names(Int)
  for(i in names(Int)){ #This loop calculates results for each interaction matrix separately.
    Res[[i]] <- dsi(Int[[i]],Phylo,Abund[[i]],Rep=Rep)
  }
  return(Res)
}

#### Partition of DSI components among locations and species ####
DSImean <- function(DSILoc){
  PairSQDist <- outer(c(DSILoc),c(DSILoc),FUN=function (x,y)(x-y)^2)
  ArrayPairDist <- array(NA,dim=c(nrow(DSILoc),nrow(DSILoc),ncol(DSILoc),ncol(DSILoc)))
  for(i in 1:ncol(DSILoc)){
    for(j in 1:ncol(DSILoc)){
      ArrayPairDist[,,i,j] <- PairSQDist[c((i-1)*nrow(DSILoc)+1):c(i*nrow(DSILoc)),c((j-1)*nrow(DSILoc)+1):c(j*nrow(DSILoc))]
    }
  }
  IntraComDist <- array(NA,dim=c(nrow(DSILoc),nrow(DSILoc),ncol(DSILoc)))
  InterComDist <- ArrayPairDist
  for(i in 1:ncol(DSILoc)){
    IntraComDist[,,i] <- ArrayPairDist[,,i,i]
    InterComDist[,,i,i] <- NA
  }
  IntraSPDist <- array(NA,dim=c(nrow(DSILoc),ncol(DSILoc),ncol(DSILoc)))
  for(i in 1:nrow(DSILoc))  {
    IntraSPDist[i,,] <- InterComDist[i,i,,]
    InterComDist[i,i,,] <- NA
  }
  IntraSP <- mean(IntraSPDist,na.rm=T)
  IntraCom <- mean(IntraComDist,na.rm=T)
  InterCom <- mean(InterComDist,na.rm=T)
  TotSQDist <- mean(PairSQDist,na.rm=T)
  return(c(IntraSP=IntraSP,IntraCom=IntraCom,Interaction=InterCom,Total=TotSQDist))
}

#### Null model to test for the different components calculated through DSImean ####
NullNA <- function(DSILoc,rep=1000){
  NullDSI <- array(NA,dim = c(dim(DSILoc),rep))
  for(i in 1:rep){
    NullDSI[,,i][!is.na(DSILoc)] <- sample(DSILoc[!is.na(DSILoc)])
  }
  NullPart <- t(apply(NullDSI,3,DSImean))
  #  NullDist <- matrix(unlist(NullPart),length(NullPart),4,byrow = T)
  #  colnames(NullDist) <- names(NullPart[[1]])
  #  return(NullDist)
  return(NullPart)
}

#### Wrapper to calculate DSImean and the corresponding null model and provide proper results ####
DSIpart <- function(Mat, rep=1000){
  OBSPart <- DSImean(Mat)
  NPart <- NullNA(Mat,rep)
  NullCI <- rbind(apply(NPart,2,quantile,probs=c(0.025,0.975), na.rm=T))
  Zscore <- (OBSPart-apply(NPart,2,mean, na.rm=T))/apply(NPart,2,sd, na.rm=T)
  Res <- list(Mat=Mat, OBS=OBSPart,Null=NPart, CI=NullCI,Z=Zscore)
  return(Res)
}

#### Plotting the partiion ####
Plot.DSIpart <- function(Part){
  library(ggplot2)
  x <- factor(names(Part$OBS),levels(factor(names(Part$OBS)))[c(3,2,1,4)])
  qplot(x[-4],Part$OBS[-4], xlab= "Component", ylab= "Mean Squares") + theme_classic() + geom_pointrange(aes(ymax=Part$CI[2,-4],ymin=Part$CI[1,-4]))
}

