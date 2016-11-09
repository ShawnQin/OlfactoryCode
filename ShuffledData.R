# This function shuffle a matrix or data frame in a column wise, row wise or both
# There are serveral method can be used to do the shuffle
# Simplest way is to use sample

ShuffledData <- function(Data,  #matrix or data frame to be shuffled
                         WhichWay = 'col', #default, column wise
                         Method = 1, #using simplest sample
                         ...
){
  #data type of input
  if(WhichWay=='col'){
    #ShuffledData <- sapply (1:nrow(Data), function (row) Data[row,]<-sample(Data[row,]))
    ShuffledMatrix <- t(apply(Data, 1, sample))
  }
  if(WhichWay=='row'){
    #ShuffledData <- sapply (1:ncol(Data), function (column) Data[,column]<-sample(Data[,column]))
    ShuffledMatrix <- t(apply(Data, 2, sample))
  }
  if(WhichWay=='both'){
    Shuffled <- matrix(Data[sample(1:(dim(Data)[1]*dim(Data)[2]))],nrow=nrow(Data),ncol = ncol(Data))
    if(is.data.frame(Data)){
      ShuffledMatrix <- as.data.frame(Shuffled)
    }
    else{
      ShuffledMatrix <- Shuffled
    }
  }
  return(ShuffledMatrix)
}