#this function set the equilibrium constants matrix so as to balance the excitatory and inhibitory response
equiConstMatx <- function(digitMatx, #digital states matrix
                          DigitSp,   #spiking rates to threshold states
                          ORNmax,    #maximum spiking rate
                          alpha,       #free engergy difference without odorant
                          method = 1, #balance method, 1 for colum balance,2 for row balance, 3 for global balance
                          ...
){
    #set the equilibrium of excitatory response
    Kexcit<- alpha*DigitSp[3:6]/(ORNmax -DigitSp[3:6])-1    
    equilMatx <- matrix(0,dim(digitMatx)[1],dim(digitMatx)[2]) #initialization
    for(i in 1:4){
        equilMatx[which(digitMatx== i)] <- Kexcit[i]
    }
    
    #set the inhibitory response equilibrium so as to balance the excitatory response
    Kinhi <- (ORNmax-DigitSp[1])/alpha/DigitSp[1]-1 #if the responses of an ORN  are inhibitory
    if(method ==1){
        #statistic of different states
        stateNum <-table(digitMatx) #which gives the number of each state in the state matrix
        ORNstatis <- apply(digitMatx,2,table)  #return a list of statitics of eahc ORN
        for(i in 1:dim(digitMatx)[2]){
            if(prod(range(digitMatx[,i]))<0){
                exitKprod<-sum((1+Kexcit[1:(length(ORNstatis[[i]])-2)])*ORNstatis[[i]][3:length(ORNstatis[[i]])])
                inhiKprod <- exitKprod/ORNstatis[[i]][1] -1
                if(inhiKprod > 0){
                    equilMatx[which(digitMatx[,i]== -1),i] <-inhiKprod
                }
                else{
                    equilMatx[which(digitMatx[,i]== -1),i] <-Kinhi
                }
            }
        else{
            equilMatx[which(digitMatx[,i]== -1),i] <- Kinhi
        }
        }
    }
    else if(method == 2){
        ORNstatis <- apply(digitMatx,1,table)  #return a list of statitics of eahc ORN
        for(i in 1:dim(digitMatx)[1]){
            if(prod(range(digitMatx[i,]))<0){
                exitKprod<-sum((1+Kexcit[1:(length(ORNstatis[[i]])-2)])*ORNstatis[[i]][3:length(ORNstatis[[i]])])
                inhiKprod <- exitKprod/ORNstatis[[i]][1] -1
                if(inhiKprod > 0){
                    equilMatx[i,which(digitMatx[i,]== -1)] <-inhiKprod
                }
                else{
                    equilMatx[i,which(digitMatx[i,]== -1)] <-Kinhi
                }
            }
            else{
                equilMatx[i,which(digitMatx[i,]== -1)] <- Kinhi
            }
        }
        
    }
    else if(method ==3){
        stateNum <-table(digitMatx)
        KinhiAll<-sum((1+Kexcit[1:(length(stateNum)-2)])*stateNum[3:length(stateNum)])/stateNum[1] -1
        equilMatx[digitMatx<0] <-KinhiAll
    }
    else if(method ==4){
        equilMatx[digitMatx<0] <-Kexcit[2]
    }
    else{
        stop("method has to be 1, 2, or 3, please try again!")
    }
    return(equilMatx)
}