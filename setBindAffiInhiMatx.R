#this function set the equilibrium constants matrix so as to balance the excitatory and inhibitory response
#the methods are revised following Prof. Tu's suggestion
#last revised on May 16,2016

setBindAffiInhiMatx <- function(digitMatx,  #digital states matrix
                          DigitSp,      #spiking rates to threshold states
                          ORNmax,      #maximum spiking rate
                          alpha,       #free engergy difference without odorant
                          Beta = 1,        #set the balance degree, default is 1
                          method = 1, #balance method, 1 for colum balance,2 for row balance, 3 for global balance
                          ResType = "multi",
                          ...
){
  
  #set the equilibrium of excitatory response
  Kexcit<- alpha*DigitSp[3:6]/(ORNmax -DigitSp[3:6])-1    
  equilMatx <- matrix(0,dim(digitMatx)[1],dim(digitMatx)[2]) #initialization
  colnames(equilMatx) <- colnames(digitMatx)
  rownames(equilMatx) <- rownames(digitMatx)
  for(i in 1:4){
    equilMatx[which(digitMatx== i)] <- Kexcit[i]
  }
  
  #set the inhibitory response equilibrium so as to balance the excitatory response
  #Kinhi <- (ORNmax-DigitSp[1])/alpha/DigitSp[1]-1 #if the responses of an ORN  are inhibitory
  if(method ==1 & ResType == "multi"){
    #statistic of different states
    #stateNum <-table(digitMatx) #which gives the number of each state in the state matrix
    #ORNstatis <- apply(digitMatx,2,table)  #return a matrix or a list of statitics of each ORN
    for(i in 1:dim(digitMatx)[2]){
        ORNstatis <- table(digitMatx[,i])
        if(prod(range(digitMatx[,i]))<0){
            exitKprod <- sum(log(1 + Kexcit[1:(length(ORNstatis) - 2)])*ORNstatis[3:length(ORNstatis)])
            inhiKprod <- exp((exitKprod-log(Beta))/ORNstatis[1]) -1
            equilMatx[which(digitMatx[,i] == -1),i] <-inhiKprod
      }
    }
  }
  
  # OSN-wise, competitive binding
  if(method ==1 & ResType == "comp"){
    for(i in 1:dim(digitMatx)[2]){
      ORNstatis <- table(digitMatx[,i])
      if(prod(range(digitMatx[,i]))<0){
        exitKprod <- sum(Kexcit[1:(length(ORNstatis) - 2)]*ORNstatis[3:length(ORNstatis)]) + 1
        inhiKprod <- (exitKprod/Beta - 1)/ORNstatis[1]
        equilMatx[which(digitMatx[,i] == -1),i] <-inhiKprod
      }
    }
  }
  
  #multiple binding, global balance
  else if(method ==2 & ResType == "multi"){
    stateNum <- table(digitMatx)
    KinhiAll <- exp((sum(log(1 + Kexcit[1:(length(stateNum)-2)])*stateNum[3:length(stateNum)]) - dim(digitMatx)[2]*log(Beta))/stateNum[1]) -1
    equilMatx[digitMatx< 0] <-KinhiAll
  }
  
  #competitive binding, global balance
  else if(method ==2 & ResType == "comp"){
    stateNum <- table(digitMatx)
    KinhiAll <- ((sum(Kexcit[1:(length(stateNum)-2)]*stateNum[3:length(stateNum)]) + 1)/dim(digitMatx)[2]/Beta - 1)/stateNum[1]
    equilMatx[digitMatx< 0] <-KinhiAll
  }
  # else{
  #   stop("method has to be 1, 2 please try again!")
  # }
  return(equilMatx)
}