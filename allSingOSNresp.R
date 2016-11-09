allSingOSNresp <-function(p0=0.2,     
                          digitMatx,  
                          DigitSp,   
                          ORNmax,      
                          alpha,    
                          beta = 1,        
                          concentration = 1, 
                          method = 1,  
                          interType="multi",
                          plotOrNot = TRUE,
                          ...
){
  ############################################################
  # this function use enumerate method to calculate the response of all mixture for a certain OSN
  # assuming the existence probability of every odorant is te same, suggested by Prof. Tu
  # p0              the probability of any given odorant apprearing in mixture
  # digitMatx       discrized interaction matrix, withe elements to be -1 ~ 4
  # DigitSp         threshold used to decide a level of spike
  # ORNmax          Maximum spiking rate
  # alpha           free engergy difference without odorant
  # beta,           set the balance of excitatory and inhibitory, default is 1
  # concentration   concentration of odorant in mixture
  # method          balance method, 1 for colum balance,2 for global
  # set the binding  affinity of 5 different response
  # interType        default interaction among odors, multiple or competitive binding
  # plotOrNot       whether plot the figure
  
  # last revised on 10/22/2016
  ################################################################################
  
  #set baseline activity and excitatory binding affinity
  basline <- ORNmax/(1+alpha);
  Kexcit<- alpha*DigitSp[3:6]/(ORNmax -DigitSp[3:6])-1 
  #stateNum <-table(digitMatx)   #which gives the number of each state in the state matrix
  # ORNstatis <- apply(digitMatx,2,table)   #return a list of statitics of eahc ORN, the number of different response 
  
  #ONS-wise balance and multiple binding sites
  if(method ==1 & interType == "multi"){
    inhiK = rep(0,times =dim(digitMatx)[2]) #initialize the inhibtory binding constant vector
    for(i in 1:dim(digitMatx)[2]){
        ORNstatis <- table(digitMatx[,i])
        exitKprod<-sum(log(1+Kexcit[1:(length(ORNstatis)-2)])*ORNstatis[3:length(ORNstatis)])
        inhiK[i] <- exp((exitKprod-log(beta))/ORNstatis[1]) -1
    }
  }
  
  #ONS-wise balance and competitive binding sites
  if(method ==1 & interType == "comp"){
    inhiK = rep(0,times =dim(digitMatx)[2]) #initialize the inhibtory binding constant vector
    for(i in 1:dim(digitMatx)[2]){
      ORNstatis <- table(digitMatx[,i])
      exitKprod<-sum(Kexcit[1:(length(ORNstatis)-2)]*ORNstatis[3:length(ORNstatis)])+1
      inhiK[i] <- ((exitKprod+1)/beta - 1)/ORNstatis[1]
    }
  }
  
  #global balance and multiple binding site model
  else if(method ==2 & interType == "multi"){
    stateNum <-table(digitMatx)
    numPari <- dim(digitMatx)[2]  #modified 10/22/2016
    KinhiAll<-exp((sum(log(1+Kexcit[1:(length(stateNum)-2)])*stateNum[3:length(stateNum)])-numPari*log(beta))/stateNum[1]) -1
    inhiK <-rep(KinhiAll,times =dim(digitMatx)[2])
  }
  
  #global balance and competitive binding model
  else if(method ==2 & interType == "comp"){
    stateNum <-table(digitMatx)
    numPari <- dim(digitMatx)[2]  #modified 10/22/2016
    KinhiAll<- ((sum(Kexcit[1:(length(stateNum)-2)]*stateNum[3:length(stateNum)])+1)/numPari/beta - 1)/stateNum[1]
    inhiK <-rep(KinhiAll,times =dim(digitMatx)[2])
  }
  
  # browser()
  # define several functions used in enumeration method
  #probability function of mixture
  pVector <-function(odor_vec,p0,N_odor){
    pVec <- prod(choose(N_odor+1,odor_vec)*p0^odor_vec*(1-p0)^(1+N_odor-odor_vec))
    return(pVec)
  }
  
  #response function of mixture
  respFunc <- function(odor_vec,ORNmax,Klist,alpha,concentration,exciType,method=1)
  {
    #default, multiple binding sites
    if(interType=="multi"){
        if(length(Klist)==5){
          resp <- ORNmax/(1+alpha*(1+concentration*Klist[1])^odor_vec[1]/prod((1+concentration*Klist[1+exciType])^odor_vec[2:length(odor_vec)]))
      # browser()
        }
    else if(length(Klist)==4){
        resp <-ORNmax/(1+alpha/prod((1+concentration*Klist[exciType])^odor_vec))
    }
    else{
      stop("K list is wrong!")
        }
    }
    #competitive binding
    else if(interType=="comp"){
         if(length(Klist)==5){
          resp <- ORNmax/(1+alpha*(1+concentration*Klist[1]*odor_vec[1])/(1+sum(concentration*Klist[1+exciType]*odor_vec[2:length(odor_vec)])))
        # browser()
          }
        else if(length(Klist)==4){
          resp <-ORNmax/(1+alpha/(1+sum(concentration*Klist[exciType]*odor_vec)))
        }
    }
    
    
    return(resp)
  }

  # entropy fucntion, calcuate the entropy of each OSN
  entropyRes <- function(allResp,#all response
                         allProb,#probability of each response
                         type = 1, #1 for original and 2 for inhibiton
                         numBin = 250, #number of bins used to calcuate the entropy
                         ...
  ){
    # rangeRes <- range(allResp)
    # cutY <- floor(rangeRes[1]):ceiling(rangeRes[2])
    if(type==1){
      cutY <- (0:numBin)*ORNmax/numBin
      # countY <-table(cut(allResp,breaks = cutY))
      countY <- factor(unclass(cut(allResp,breaks = cutY)))
      # sel_count<-countY[countY>0]
    }
    else if(type==2){
      cutY<- basline-0.1 + (0:numBin)*(ORNmax-basline+0.1)/numBin
      # countY <-table(cut(allResp,breaks = cutY))
      countY <- factor(unclass(cut(allResp,breaks = cutY)))
      # sel_count<-countY[countY>0]
    }
    # countY <- hist(allResp,cutY,plot=FALSE)
    # sel_count <- countY$counts[countY$counts>0]
     # browser()
    probFrame <-data.frame(bins=countY,
                           prob=allProb)
    newPro<-tapply(probFrame$prob, probFrame$bins,FUN=sum)
    ep <- -sum(log2(newPro[newPro>0])*newPro[newPro>0])
    return(ep)
  }  
 
  #statistics of each OSN,number of different odoratns
  library(combinat)
  # ORNstatis <- apply(digitMatx,2,table)
  entropyOri <-rep(0,times=dim(digitMatx)[2])
  entropyNoInhi<-rep(0,times=dim(digitMatx)[2])
  
  for(i in 1:dim(digitMatx)[2]){
    ORNstatis <- table(digitMatx[,i])
    if(any(as.numeric(names(ORNstatis))<0)){
      N_odor <- ORNstatis[c(1,3:length(ORNstatis))]  #number of odorants for each response type
      Klist <- c(inhiK[i],Kexcit)  #set the binidng affinity of this receptor
      N_odor_noInhi <- ORNstatis[3:length(ORNstatis)] #for without inhibition situation
    }
    else{
      N_odor <-ORNstatis[2:length(ORNstatis)]  #only excitatory response
      Klist <- Kexcit  #set the binidng affinity of this receptor
      N_odor_noInhi <- ORNstatis[2:length(ORNstatis)] #for without inhibition situation
    }
    
    # browser()
    allType <- as.numeric(names(N_odor))
    exciType <- allType[allType>0]
    # if(length(exciType) !=length(N_odor)){
    #   browser()
    # }
    n1all <- hcube(N_odor+1,1)-1    #all the combination of odorants for a certain OSN
    N_tot <- prod(N_odor+1)     #number of total sample for this OSN
    
    
    allProb <-apply(n1all,MARGIN = 1, FUN=function(x) pVector(x,p0,N_odor))
    # browser()
    # debug(respFunc)
    if(length(exciType)==0){
      print(i)
      stop("wrong!")
    }
    allResp <- apply(n1all, MARGIN=1, FUN=function(y) respFunc(y,ORNmax,Klist,alpha,concentration,exciType,method=interType))
  # debug(entropyRes)
    entropyOri[i] <-entropyRes(allResp,allProb,type = 1)
    # browser()
    
    # N_odor_noInhi <- ORNstatis[[i]][3:length(ORNstatis[[i]])] #for without inhibition situation
    n1all_noInhi <- hcube(N_odor_noInhi+1,1)-1    #all the combination of odorants for a certain OSN
    # browser()
    Klist_noInhi <- Kexcit    #no inhibition
    allProbNoInhi <-apply(n1all_noInhi, MARGIN=1, FUN=function(x) pVector(x,p0,N_odor_noInhi))
    allRespNoInhi <- apply(n1all_noInhi, MARGIN=1, FUN=function(y) respFunc(y,ORNmax,Klist_noInhi,alpha,concentration,exciType,method=interType))
    # browser()
    entropyNoInhi[i]<-entropyRes(allRespNoInhi,allProbNoInhi,type = 2)
  }

  allEntr <- list(entrOri=entropyOri,entroNoInhi=entropyNoInhi)
  # browser()
  if (is.null(colnames(digitMatx))){
    ColNames <- c(1:dim(digitMatx)[2])
  }
  else{
    ColNames <- colnames(digitMatx)
  }
  allEntropy <- data.frame(
    OrnNames = rep(factor(ColNames,as.character(ColNames)),times=2),
    system = factor(rep(c("Ori", "NoI"), each= dim(digitMatx)[2])),
    entropyORN = c(entropyOri,entropyNoInhi)
  )
  
  # plot entropy of each OSN
  if(plotOrNot){
      gplot <- ggplot(allEntropy,aes(x = OrnNames,y = entropyORN,fill= system)) + geom_bar(stat="identity",  position = position_dodge()) + labs(x='ORN',y="entropy (bit)") + theme(axis.text.x = element_text(size = 14,angle = 60, hjust = 1))
      entropyFile <- paste("newMethod",as.character(method),"_sparse",as.character(p0),"_",interType,"_entropy.pdf",sep = "")
      ggsave(gplot, file=entropyFile, width=10, height=5)
  }
  
  return(allEntr)
}