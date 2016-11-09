#this function calculate the mutual information between selected ORN to odor mixtures
#using smapling method
#last revised on May 17

mutualInfo <-function(p0=0.2,     #the probability of any given odorant apprearing in mixture
                      weightMatx,  #digital states matrix
                      equilMatx,    #spiking rate of typic response
                      ORNmax,      #maximum spiking rate
                      alpha,       #free engergy difference without odorant
                      beta = 1,        #set the balance degree, default is 1
                      concentration = 1, #concentration of odorant in mixture
                      method = 1, #balance method, 1 for colum balance,2 for row balance, 3 for global balance
                      interType="multi",  #default interaction among odors
                      numSamp=10000, #default sampling times
                          ...
){
  library(infotheo)
  #set the binding affinity of 5 different response
  basline <- ORNmax/(1+alpha);
  
  #mixture response fucntion 
  mixResFunc <-function(weightMatx,selectInx,inx,iny,interType){
    K <- equilMatx[,c(inx,iny)]
    W <- weightMatx[,c(inx,iny)]
    if(interType=="multi"){
      #multiple binding sites
      interac <- c(1,1)  #initialize the response
      # browser()
      for(i in 1:length(selectInx)){
        interac <- interac*(1 + concentration*K[selectInx[i],])^(-W[selectInx[i],])
      }
      response<- ORNmax/(1+alpha*interac)
      return(response)
    }
    else if(ResType=="comp"){
      reac_ini<- c(0,0)  #initialize the response
      reac_exc<- c(0,0)  #initialize the response
      for(i in 1:2){
        for(j in 1:length(selectInx)){
          if(W[selectInx[j],i] > 0){
            reac_exc[i] = reac_exc[i] + concentration*K[selectInx[j],i]
          }
          else if(W[selectInx[j],i] < 0){
            reac_ini[i] = reac_ini[i] + conc*K[selectInx[j],i]
          }
        }
      }
      response <- ORNmax/(1+alpha*(1+reac_ini)/(1+reac_exc))
      return(response)
    }
    
  }

  #calculate the mutual information using package entropy  or infotheo
  muFunc <- function(allResponse,numBin,withInhi=1){
    #library(entropy)
    
    # numBin <- 250
    if(withInhi==1){
      #with inhibiton 
      rangRes <- c(0,ORNmax)
    }
    else if(withInhi==2){
      #without inhibition
      rangRes <- c(29.9,ORNmax)
    }
    disc <- discretize(allResponse,disc="globalequalwidth", nbins=numBin)
    #browser()
    MutualInformation <- mutinformation(disc[,1],disc[,2])
    
    #disc<-discretize2d(allResponse[,1], allResponse[,2], numBins1=numBin, numBins2=numBin,r1=rangRes, r2=rangRes)
    #mutualInfo<-mi.empirical(disc,unit = "log2")
    return(MutualInformation)
  }
  
  #only sampling among the odors that can elicit either response of ORN
  #set the interaction matrix for without inihbition
  weightMatxNoInhi <-weightMatx
  weightMatxNoInhi[weightMatxNoInhi<0] <- 0

allMutualInfo <- matrix(NaN,dim(weightMatx)[2],dim(weightMatx)[2])
allMutualInfoNoInhi <- matrix(NaN,dim(weightMatx)[2],dim(weightMatx)[2])
numBin = 25
for(l in 1:(dim(weightMatx)[2]-1)){
    for(m in (l+1):dim(weightMatx)[2] ){
        inx <-2
        iny <-9
        
        #find index that have at leat one response
        temp <-apply(abs(weightMatx[,c(inx,iny)]),1,sum)
        INX <-which(temp>0)
        numOdor <- length(INX)
        allResponse <- matrix(0,numSamp,2)
        
        for(i in 1:numSamp){
            samp <- rbinom(numOdor,1,p0)
            if(all(samp==0)){
                allResponse[i,] <- ORNmax/(1+alpha)*c(1,1)
            }
            else{
                selectInx <- INX[as.logical(samp)]
                allResponse[i,] <- mixResFunc(weightMatx,selectInx,inx,iny,interType)
            }
            
        }
        
        allMutualInfo[l,m] <- muFunc(allResponse,numBin,withInhi = 1)
        
        #without inhibition
        tempNo <-apply(abs(weightMatxNoInhi[,c(inx,iny)]),1,sum)
        INXno<- which(tempNo>0)
        numOdorNo <- length(INXno)
        allResponseNoInhi <- matrix(0,numSamp,2)
        for(i in 1:numSamp){
            samp <- rbinom(numOdorNo,1,p0)
        if(all(samp==0)){
            allResponseNoInhi[i,] <- ORNmax/(1+alpha)*c(1,1)
            }
        else{
            selectInx <- INXno[as.logical(samp)]
            allResponseNoInhi[i,] <- mixResFunc(weightMatxNoInhi,selectInx,inx,iny,interType)
            }
        }
        browser()
        allMutualInfoNoInhi[l,m] <- muFunc(allResponseNoInhi,numBin,withInhi = 2)
    }
}
browser()
write.xlsx(allMutualInfo,"mutualInfoMatx.xlsx")
write.xlsx(allMutualInfoNoInhi,"mutualInfoMatxNoInhi.xlsx")
  
  
#without inhibition


#   gplot <- ggplot(allEntropy,aes(x = OrnNames,y =entropyORN,fill= system))+geom_bar(stat="identity", position=position_dodge())+labs(x='ORN',y="entropy (bit)") + theme(axis.text.x=element_text(size = 16,angle = 60,hjust = 1))
#   if(interType==1){
#     intTy<-"multi"
#   }
#   else if(interType==2){
#     intTy<-"compete"
#   }
#   entropyFile <- paste("newMethod",as.character(method),"_sparse",as.character(p0),"_",intTy,"_entropy.pdf",sep = "")
#   # browser()
#   ggsave(gplot, file=entropyFile, width=10, height=5)
#   return(allEntr)
}