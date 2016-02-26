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


#### Wrapper to calculate DSImean and the corresponding null model and provide proper results ####
DSIpart <- function(Mat, rep=1000){
  OBSPart <- DSImean(Mat)
  NPart <- NullNA(Mat,rep)
  NullCI <- rbind(apply(NPart,2,quantile,probs=c(0.025,0.975), na.rm=T))
  Zscore <- (OBSPart-apply(NPart,2,mean, na.rm=T))/apply(NPart,2,sd, na.rm=T)
  Res <- list(Mat=Mat, OBS=OBSPart,Null=NPart, CI=NullCI,Z=Zscore)
  return(Res)
}
