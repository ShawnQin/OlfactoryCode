RandomizedMatrx <- function(OldMatrix,method="switch",...){
  # this function return randomized interaction network
  # matrix     old interaction network
  # method     mathod used to randomize the matrix, can be "switch","bnOdorWise" and "bngloble"
  # ====================================
  if(method=='switch'){
    MaxTry <- 100*sum(OldMatrix != 0) # maximum number of try
    # all interaction pairs index
    allInx <- which(OldMatrix !=0,arr.ind = TRUE)
    LEN <- dim(allInx)[1]  #total number of edges
    # ValueInter <- OldMatrix[allInx]  # store the interaction type or value
    # RandomizedMatrx <- matrix(0,nrow = nrow(OldMatrix),ncol = ncol(OldMatrix))
    # select pairs of interaction
    count <- 0  # how many times has been tried
    UpdateMat <- OldMatrix
    while (count <= MaxTry) {
      sampleFlag <- 1
      while(sampleFlag){
      temp <- sample(LEN,2)  # which two rows
      InxTry <- allInx[temp,]
      if(InxTry[1,1] != InxTry[2,1] && InxTry[1,2] != InxTry[2,2]){
        sampleFlag <- 0 
      }
      }
      
      # switch or not
      if(UpdateMat[InxTry[1,1],InxTry[2,2]] == 0 && UpdateMat[InxTry[2,1],InxTry[1,2]] == 0){
        allInx[temp[1],] <- c(InxTry[1,1],InxTry[2,2])
        allInx[temp[2],] <- c(InxTry[2,1],InxTry[1,2])
        UpdateMat[allInx[temp,]] <- UpdateMat[InxTry]
        UpdateMat[InxTry] <- 0
        # ValueInter[temp] <- ValueInter[c(temp[2],temp[1])]
      }
      count <- count + 1
      # RandomizedMatrx <- matrix(0,nrow = nrow(OldMatrix),ncol = ncol(OldMatrix))
      # RandomizedMatrx[allInx] <- ValueInter
      # if(any(colSums(abs(UpdateMat)) != colSums(abs(OldMatrix))) | any(rowSums(abs(UpdateMat)) != rowSums(abs(OldMatrix)))){
      #   browser()
      # }
    }
  }
    return(UpdateMat)
}