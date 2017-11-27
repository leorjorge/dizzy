#### MaxGen Functions, used to calculate the maximum possible MPD for a given number of individuals in the phylogeny.####
#The function uses as input a distance matrix and the interaction data. It should be used within the DSI calculation function and measures the maximum values of MPDs for a given list of consumers in a given Interaction array. This function depends on the other functions below, that implement the simulated annealing optimization
#It returns the Maximum MPD for each herbivore species

MaxGenReg <- function(Phy, Int, Spps){
  MaxMPDReg <- vector(mode = "numeric", length = length(Spps))
  Samp <- apply(Int, 1, sum)
  ConsLoc <- apply(Int, c(1,3), sum)
  ConsLoc <- ConsLoc > 0
  ResLoc <- apply(Int,c(2,3),sum)
  for(i in 1:length(MaxMPDReg)) {
    if (Samp[Spps[i]]<2){
      MaxMPDReg[i] <- NA
    } else {
      ResSums <- rowSums(as.matrix(ResLoc[,ConsLoc[Spps[i],]]))
      Phy.sp <- Phy[rownames(Phy) %in% names(ResSums[ResSums > 0]), 
                    colnames(Phy) %in% names(ResSums[ResSums > 0])]
      MaxMPDReg[i] <- MaxMPD(Samp[Spps[i]],Phy.sp)$value*-1
    }
  }
  return(MaxMPDReg)
}


MaxGenLoc <- function(Phy, Int, SpLocs){
  MaxMPDLoc <- vector(mode = "numeric", length = nrow(SpLocs))
  Samp <- apply(Int, c(1,3), sum)
  Abund <- apply(Int, c(2,3), sum)
  for (i in 1:length(MaxMPDLoc)) {
    if (Samp[SpLocs[i]] < 2) {
      MaxMPDLoc[i] <- NA
    } else {
      LocAbund <- Abund[,SpLocs[i,2]]
      Phy.sp <- Phy[rownames(Phy) %in% names(LocAbund[LocAbund > 0]), 
                    colnames(Phy) %in% names(LocAbund[LocAbund > 0])]
      MaxMPDLoc[i] <- MaxMPD(Samp[SpLocs[i,1],SpLocs[i,2]],Phy.sp)$value*-1
    }
  }
  return(MaxMPDLoc)
}

### The MPD maximization function, employing simulated annealing. It is called by both MaxGenReg and MaxGenLoc, and depends on genAbund, which is the transition function between different distributions of resource items in the phylogeny, and MPD, which is the function that calculates MPD among these individuals. This very modular design is demanded by the simulated annealing implementation in the optim function.

MaxMPD <- function(N, Phylo){
  Prob <- 1/(2^adephylo::distRoot(ape::compute.brlen(ape::as.phylo(hclust(as.dist(Phylo)))), 
                                  method = "nNode"))
  Start <- t(table(factor(sample(rownames(Phylo), prob = Prob, size = N, replace = T), 
                          levels = rownames(Phylo))))
  res <- optim(par = Start, fn = MPD, gr = genAbund, Dist = Phylo,
               method = "SANN", control = list(maxit = 150000, temp = 300, trace = 0))
  return(res)
}


genAbund <- function(Start, Dist) { 
  idx <- 1:length(Start)
  LossAbund <- sample(x = idx[Start > 0], size = 1, replace = FALSE)
  GainAbund <- sample(x = idx[-LossAbund], size = 1, replace = FALSE)
  Start[LossAbund] <- Start[LossAbund] - 1
  Start[GainAbund] <- Start[GainAbund] + 1
  Start
}


#####MPD####
#Function to calculate MPD correctly. The rationale is similar to the MPD function in picante, but the diagonal in the sample weights matrix is corrected, using the number of combinations of diferent individuals in the same species at the diagonal. MPD for 1 single species is also changed to 0 instead of NA. (mpd2)
#In addition, a function using only numeric arguments is also available, for use in the optimization of MaxMPD (MPD)
#I'm considering changing everything in the package to use MPD.

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
      } else {
        mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
      }
    } else {
      mpd[i] <- 0
    }
  }
  mpd
}


MPD <- function(Start, Dist){
  Dist <- matrix(Dist, nrow = length(Start))
  Dist.samp <- Dist[which(Start>0),which(Start>0)]
  Start <- Start[Start>0]
  sample.weights <- outer(Start,Start,"*")
  diag(sample.weights) <- as.numeric(Start*(Start-1)/2)
  weighted.mean(Dist.samp[lower.tri(Dist.samp, diag=TRUE)], sample.weights[lower.tri(sample.weights, diag=TRUE)])*-1
}
