#partion the response matrix into digitalized matrix
#this function tranform the response matrix into digitalized response matrix
digitMatx <- function(X,    #response matrix
                      DigitSp, #digits responses
                      ORNmax,
                      ...
){
    newMatx <- matrix(0,dim(X)[1],dim(X)[2])   #initialize
    allRS <- c(DigitSp,ORNmax)  #be care of the odor
    newMatx[X<=DigitSp[1]] <- -1#inhibition
    for(i in 4:length(allRS)){
        newMatx[X>=allRS[i-1]&X <allRS[i]] <- i-3       
    }
    return(newMatx)
}
