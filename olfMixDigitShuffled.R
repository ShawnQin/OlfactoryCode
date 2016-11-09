#*************************************************************
#GROGRAM: olfMixDigitShuffled.R
#DESCRIPTION:
#this program calculate the ORN response patterns of odor mixtures
#for the sytem with inhibitory response and without inhibitoyr response
#the integration rule at ORN is proposed by Prof. Tu
#it differs from olfMixDsic by
# 1) same baseline of all ORN
# 2) responses are digitalized from -1 to 5

#modifed from  "olfMixDigitDis.R", here the interaction matrix was randomly shuffled before 
#doing the calculation
#last revised on 10/20/2016
#*******************************************************************

#load liberies
library(xlsx)
library(RColorBrewer) #this package is used for color settings
library(ggplot2)

#source the functions
source("setBindAffiInhiMatx.R")
source("mixRespFunc.R")
source("mixResp.R")
source("pcaSpectrum.R")
source("digitMatx.R")
source("plotHistCompORN.R")
source("plotMeanCorr.R") 
source("compAll.R")
source("plotPCAEntropy.R")#plot the summary of entropy at each OSN
source("histCount.R")     #hist gram of response pattern on OSN
source("ShuffledData.R")  #shuffling the interaction matrix
source("PlotComparePermutation.R")  #plot statistics of interaction matrix
source("allSingOSNresp.R")

#load data file
digiFile<-"data/CarslonORNdigit.xlsx"
rawData1<-read.xlsx(digiFile,1)
rawMatx <- as.matrix(rawData1[,3:length(rawData1)])
rownames(rawMatx)<- rawData1$odorantName[1:110]

#delete odorants that don't elicite any response
adjustSpikingMatx<-rawMatx     #estimated absolute ORNs spiking rate
temp<-apply(abs(adjustSpikingMatx),1,sum)
newMatx <- adjustSpikingMatx[which(temp >0),]

#delete the ORNs that only have inhibitory response or no repoonse
tempORN <- apply(apply(newMatx,2,range),2,sum)
newMatx <-newMatx[,-which(tempORN==-1)]


#basic parameters in the 
ORNmax <- 250      #inorder to get the value of Kij
spontaneous <- 30  #assuming 10% basal level response
alpha <- ORNmax/spontaneous - 1
OdorNum <- dim(newMatx)[1] #number of odorants

#weight matrix before permuation
# oldWeightMatrix <- newMatx
# oldWeightMatrix[oldWeightMatrix>0] <- 1
# oldWeightMatrix[oldWeightMatrix<0] <- -1

##compare the interaction matrix before and after permutation
#PlotComparePermutation(oldWeightMatrix,weightMatx)
#assuming the concentrations are 1
DigitSp<- c(15,30,50,100,150,200)
beta <-1
balanceMethod <-1 #method used to balance
cCon <-  1 #concnetration
Rtype <- "comp"  #competive or multiple binding
compNumMix <- c(5,20,35)
# compNumMix <- 35
interactionType = "comp"  # 1 for multiple binding sites and 2 for competitive binding
sparsePo = (1:12)*5/107

NumTry <- 5  #try 10 times
# AllWeightMatrix <- array(0,c(dim(newMatx)[1],dim(newMatx)[2],NumTry))
# AllWeightMatrixNoInhi <- array(0,c(dim(newMatx)[1],dim(newMatx)[2],NumTry))
AllSumEntropy <- array(0,c(length(sparsePo),2,NumTry))
AllEffectiveDim <- data.frame(system=character(),numMix=numeric(),x=numeric())
startTime <- proc.time()
for (i in 1:NumTry){
  
    #make sure the shuffled matrix has at leat one positive interaction
    flag =1
    while(flag==1){
      shuffedMatx <- ShuffledData(newMatx)
      if (all(apply(shuffedMatx, 2, max) > 0)){
        flag =0
      }
    }
    
    weightMatx <- shuffedMatx
    weightMatx[weightMatx>0] <- 1
    weightMatx[weightMatx<0] <- -1
    # AllWeightMatrix[,,i] <- weightMatx

    #weightMatxNoInhi
    weightMatxNoInhi <-weightMatx
    weightMatxNoInhi[weightMatxNoInhi<0] <- 0
    # AllWeightMatrixNoInhi[,,i] <- weightMatxNoInhi

    # bindig constants matrix
    # equilMatx <- setBindAffiInhiMatx(shuffedMatx,DigitSp,ORNmax,alpha,beta,balanceMethod)
    # equilMatxNoInhi <-equilMatx
    # equilMatxNoInhi[shuffedMatx<0] <- 0

    #pca and entropy

    pcaEntr <-compAll(method = balanceMethod,compNumMix,compCon = cCon,ResType = Rtype,plotOrNot = FALSE,spontaneous,ORNmax,ORNmax,interMatx =shuffedMatx)

    # pcacmp<-subset(pcaEntr$PCAresult,numMix %in% c(5,20,35))
    pcacmp<- pcaEntr$PCAresult[c("system","numMix","eigenvalue")]
     effectivDim <- aggregate(pcacmp$eigenvalue,by=list(system = pcacmp$system,numMix=pcacmp$numMix),function(x) sum(x>=1))
    AllEffectiveDim <- rbind(AllEffectiveDim,effectivDim)
    sumEntropy<-pcaEntr$sumEnt
    plotPCAEntropy(compNumMix,sumEntropy,pcaEntr$PCAresult,balanceMethod,cCon,Rtype)


    ## using enumeration method to calculate total entropy
    sumEntropy = matrix(0,length(sparsePo),2)
    startTime <- proc.time()
    for(j in 1:length(sparsePo)){

        allEntropy <-allSingOSNresp(p0=sparsePo[j],shuffedMatx,DigitSp,ORNmax,alpha,beta = 1,concentration = 1,method = 1,interType=interactionType, plotOrNot = FALSE)

        sumEntropy[j,]<-c(sum(allEntropy$entrOri),sum(allEntropy$entroNoInhi))
    }
    proc.time() - startTime
    # plotPCAEntropy((1:12)*5,sumEntropy)
    AllSumEntropy[,,i] <- sumEntropy
    
}

myColor <- brewer.pal(11,"Spectral")
## plot effective coding dimension
SummaryDim <- aggregate(AllEffectiveDim$x,by=list(sys = AllEffectiveDim$system,num = AllEffectiveDim$numMix),function(x) c(dimen = mean(x),std = sd(x)))

#plot
gEffDim <-ggplot(SummaryDim,aes(x=num,y=x[,1],group=sys,color = sys)) + geom_pointrange(aes(ymin=x[,1]-x[,2], ymax=x[,1]+x[,2])) + geom_line() + geom_point()+scale_color_manual(values=c(myColor[10],myColor[1]))
gEffDim <- gEffDim + labs(x="number of odor", y = "effective coding dimension")
effectiveDimFileName <- paste("EffectiveCodingDimShuffled",as.character(balanceMethod),"_",Rtype,".pdf",sep = "")
ggsave(effectiveDimFileName,width = 4,height = 3)

## plot average total entropy
AverageTotalEntroy <- apply(AllSumEntropy, c(1,2), mean)  
StdTotalEntropy <- apply(AllSumEntropy, c(1,2), sd)  #standard deviation

#data frame used to plot
TotEntrStat <- data.frame(
  sys = factor(c(rep("Ori", times= length(sparsePo)),rep("Noi",times = length(sparsePo)))),
  con = factor(rep((1:12)*5, times = 2)),
  average = c(AverageTotalEntroy[,1],AverageTotalEntroy[,2]),
  std = c(StdTotalEntropy[,1],StdTotalEntropy[,2])
    )
gTotEntr <- ggplot(TotEntrStat,aes(x=con,y=average,group=sys,color=sys)) + geom_pointrange(aes(ymin=average-std, ymax=average+std)) + geom_line() + geom_point()+scale_color_manual(values=c(myColor[10],myColor[1]))
gTotEntr + labs(title="comparing total entropy", x="average number of odor", y = "total entropy(bit)")
CompareTotEntropyFileName <- paste("TotalEntropyShuffled_bm",as.character(balanceMethod),"_","ntry",as.character(NumTry),Rtype,".pdf",sep = "")
ggsave(gTotEntr,width = 4,height = 3)


