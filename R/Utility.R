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
