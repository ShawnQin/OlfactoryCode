#this program plot the odor-OSN interaction distribution
#we are interested to find the motif and the sparsity of these interaction
#we use the discreted response matrix from Hallem et al., Cell, 2006


#load liberies
library(xlsx)
library(RColorBrewer) #this package is used for color settings
library(ggplot2)

#load the data
digiFile<-"/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/CarslonORNdigit.xlsx" 
digitResponse<-read.xlsx(digiFile,1)
numMatrix <- digitResponse[,3:dim(digitResponse)[2]]

InhibitoryRatio = sum(numMatrix<0)/dim(numMatrix)[1]/dim(numMatrix)[2]
ActivationRatio = sum(numMatrix>0)/dim(numMatrix)[1]/dim(numMatrix)[2]

#how many OSN a given odor can ativate or inhibit
ActiNumVec <- apply(numMatrix,MARGIN=1,function(x) sum(x>0))
InhiNumVec <- apply(numMatrix,MARGIN=1,function(x) sum(x<0))

#how many odors a given OSN response
ActiOdorNum <- apply(numMatrix,MARGIN=2,function(x) sum(x>0))
InhiOdorNum <- apply(numMatrix,MARGIN=2,function(x) sum(x<0))

#generate binomial distribution with the same probability of excitatory and inhibitory interactions
ROW <- dim(numMatrix)[1]
COL <- dim(numMatrix)[2]
RandomBinoMatrixExci <- matrix(rbinom(ROW*COL,1,ActivationRatio), ROW,COL)
RandomBioVecExci <- apply(RandomBinoMatrixExci,MARGIN=1,function(x) sum(x>0))

RandomBinoMatrixInhi <- matrix(rbinom(ROW*COL,1,InhibitoryRatio), ROW,COL)
RandomBioVecInhi <-apply(RandomBinoMatrixInhi,MARGIN=1,function(x) sum(x>0))

#plot the distribution and compare it with random interaction
#first put all the frequence in a data frame
AllFreq <- data.frame(
        OSNActi = hist(ActiNumVec,plot = FALSE,breaks = 0:20)$counts,
        OSNActiRandom = hist(RandomBioVecExci,plot = FALSE,breaks = 0:20)$counts,
        OSNInhi = hist(InhiNumVec,plot = FALSE,breaks = 0:20)$counts,
        OSNInhiRandom = hist(RandomBioVecInhi,plot = FALSE,breaks = 0:20)$counts
    )

myColor <- brewer.pal(11,"Spectral")
plot(0:19,AllFreq$OSNActi,type = "b",pch=16,lty=1,lwd=2,cex=2,col=myColor[1],ylim = yRange,xlab = "Number of OSN",ylab = "Frequency")
lines(0:19,AllFreq$OSNActiRandom,type = "b",pch=16,lty=1,lwd=2,cex=2,col=myColor[10])