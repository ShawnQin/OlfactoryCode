#this grogram try to organized all the results, especially the PCA for different number of mixtures
compAll <- function(method = 1,#the balancing method
                    compNumMix = FALSE,#number of odorants in mixture
                    compCon = 1, #whether compare different concentration
                    ResType = "multi",  #interaction type of odor mixture on ORN
                    plotOrNot = TRUE,  #whether plot hisogram and correlation matrix
                    baseline = 30,  #basal activity
                    ORNmax = 250,
                    numBin=250,#bins used to calculate the entropy of OSN
                    interMatx = newMatx, #interaction matrix, default
                    ...                    
){
  ##=============================================###
  # this function use sampling method to do the calcuation and analysis, inculuding the PCA and entropy
  # method           used to balance the excitatory and inhibitory response, 1 for OSN-wise, 2 for global
  # compCon          whether use different concentration
  # ResType          interaction type, multiple or competitive binding
  # plotOrNot        whether plot the figure
  # baseline         basal activity
  # ORNmax           maximum response
  # interMatx        interaction matrix, digitalized
  #
  # last revised 10/22/2016
  #####################################################
  
  
    sampleSize <- 100000 #sampling size, defaulte
    if(compCon == FALSE){
        concentration =1  #concentration of odorant
    }
    else{
        concentration = compCon
    }
    
    #data frame to store the returned PCA results
    PCAallVarEigen <- data.frame(system=NULL,numMix=NULL,component=NULL,concen=NULL,variance_expl=NULL,eigenvalue=NULL)
    
    # set the binding affinity matrix for the orignial and without inhibition
    # browser()
    equilMatx <- setBindAffiInhiMatx(interMatx,DigitSp,ORNmax,alpha,beta,method,ResType)
    equilMatxNoInhi <-equilMatx
    equilMatxNoInhi[interMatx<0] <- 0
    
    #initialized entropy matrix
    sumEntropy = matrix(0,length(compNumMix),2)
    colnames(sumEntropy)<-c("Ori","No")
    
    #cycle through different number of odorants in mixture
    for(i in 1:length(compNumMix)){ 
        allResp<- mixResp(compNumMix[i],sampleSize,ORNmax,alpha,weightMatx,equilMatx,concentration,ResType=ResType)
        allRespNoInhi <- mixResp(compNumMix[i],sampleSize,ORNmax,alpha,weightMatxNoInhi,equilMatxNoInhi,concentration,ResType=ResType)
        
        digiRespMatx <- digitMatx(allResp,DigitSp,ORNmax)
        digiRespMatxNoInhi <- digitMatx(allRespNoInhi,DigitSp,ORNmax)
        fname <- paste(as.character(compNumMix[i]),"-con",as.character(concentration),"-method",as.character(method),"_",ResType,sep = "")
      
        
        #pca analysis is carried out to the digit response matrix
        pcaVarEig <- pcaSpectrum(digiRespMatx,digiRespMatxNoInhi,numMix,scale=FALSE,fname,plot=FALSE)
        
        PEPtimes<-c(length(pcaVarEig[[1]]),length(pcaVarEig[[3]]))
        
        newRows <- data.frame(
                          system=factor(rep(c("Ori", "Noi"), times= PEPtimes)),
                          numMix=factor(rep(c(compNumMix[i],compNumMix[i]), times = PEPtimes)),
                          component=c(1:PEPtimes[1],1:PEPtimes[2]),
                          concen=factor(rep(c(concentration,concentration), times = PEPtimes)),
                          variance_expl = c(pcaVarEig$variance,pcaVarEig$varianceNoInhi),
                          eigenvalue = c(pcaVarEig$eig,pcaVarEig$eigNoInhi)
                          )
        PCAallVarEigen <-rbind(PCAallVarEigen,newRows)
                          
        #plot the results
        if(plotOrNot==TRUE){
          # histgram of the response of each OSN
          saveName <- paste("histoCompareDigt",as.character(compNumMix[i]),"-con",as.character(concentration),"-method",as.character(method),"_",ResType,".pdf",sep = "")
          plotHistCompORN(allResp,allRespNoInhi,fileName=saveName)
          
          #plot the sum of entropy
           # pdf(paste("compare_sum_entropy_method",as.character(method),"_",as.character(ResType),".pdf",sep = ""),width = 6,height = 6)
           # plot(compNumMix,sumEntropy)
                          
        }
        
        #plot correlation matrix and entropy of each OSN
        fileStore <- paste("odorMix",as.character(compNumMix[i]),"-con",as.character(concentration),"-method",as.character(method),"_",ResType,sep = "")
        if (is.null(colnames(interMatx))){
          ORNname <- c(1:dim(interMatx)[2])
        }
        else{
          ORNname <- colnames(interMatx)
        }
        allEntro <- plotMeanCorr(allResp,allRespNoInhi,ORNname,fileStore,plotOrNot,baseline,ORNmax,numBin)
        # browser()
        sumEntropy[i,] <-apply(allEntro, 2, sum)

    }
    
    returnList <-list(PCAresult=PCAallVarEigen,sumEnt=sumEntropy)
    return(returnList)
}