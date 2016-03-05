#' Convert interaction records to array
#' 
#' This function converts records for the occurrence of interactions registered 
#' in 'long' format, with columns for consumer, resource, location and frequency
#' into an interaction array.
#' 
#' @export
#' 
#' @param data A dataframe of recorded interactions. It needs at least one 
#'   column for consumer identity and another for resource identity. If the
#'   study spans different locations, one column with location identity must be
#'   provided. The frequency of a given interaction is recorded in a separate
#'   column or as different records.
#'   
#' @param cons,res Consumer and resource column identity. Can be provided as 
#'   the column number or name. If those are not provided, it is assumed that 
#'   consumers are in the first column and resources in the second
#'   
#' @param loc Location column identity, also column number or name. If not 
#'   provided, a warning is given and it is assumed that interactions occur in a
#'   single location
#'   
#' @param freq Interaction frequency column identity, also column number or
#'   name. If the frequency of a given interaction is recorded as a separate
#'   column instead of being placed as different records in the dataframe, this 
#'   information must be provided.
#'   
#' @return An object of class interaction array. To calculate DSI, resource 
#'   names must match the phylogeny tip labels and resource availability names.
#'   

RecordToArray <- function (data, cons=1, res=2, loc, freq=NULL){
  if (is.null(loc)) {
    warning ("No location column - assuming one single location")
    loc <- "Location"
    data$loc <- "Location"
  }
  if (!is.null(freq)) {
    raw <- data.frame(data[,cons],data[,res],data[,loc],data[,freq])
    data <- raw[rep(1:nrow(raw), times=raw$freq),]
  }
int.array <- table(data[,cons], data[,res], data[,loc]) 
return(int.array)
}

#' Organise a set of interaction matrices into one array
#' 
#' Organises a set of local interaction matrices into one array for interactions
#' in all locations
#' 
#' @export
#' 
#' @param mats A vector with matrix object names.
#' 

MatToArray <- function (mats){
  Mat.list <- lapply(X = mats, FUN = get)
  Rows <- levels(as.factor(unlist(lapply(Mat.list,rownames))))
  Cols <- levels(as.factor(unlist(lapply(Mat.list,colnames))))
  int.array <- array(0, dim = c(length(Rows), length(Cols), length(Mat.list)))
  dimnames(int.array) <- list(Rows,Cols,mats)
  for (i in 1:length(Mat.list)) {
    int.array[match(rownames(Mat.list[[i]]), dimnames(int.array)[[1]]),
              match(colnames(Mat.list[[i]]), dimnames(int.array)[[2]]), i] <- 
      Mat.list[[i]]
  }
  return(int.array)
}