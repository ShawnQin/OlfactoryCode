#*****************************************************************
#GROGRAM: olfMixDisc.R
#DESCRIPTION:
#this program calculate the ORN response patterns of odor mixtures
#for the sytem with inhibitory response and without inhibitoyr response
#the integration rule at ORN is proposed by Prof. Tu
#writen by ssqinpku@163.com
#last revised on Apr. 25,2016
#*******************************************************************

#load liberies
library(xlsx)
library(RColorBrewer) #this package is used for color settings
library(ggplot2)

#load data file
filePath<-"/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/Carlson20016TableS1.xlsx"
digiFile<-"/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/CarslonORNdigit.xlsx"
rawData1<-read.xlsx(filePath,1)   #with the last line as the 
digitResponse<-read.xlsx(digiFile,1)
rawMatx <- as.matrix(rawData1[-nrow(rawData1),3:length(rawData1)])
rownames(rawMatx)<- rawData1$odorantName[1:110]

minSpike<-apply(rawMatx,2,min)
basalSpon <- rawData1[111,3:length(rawData1)]
realSpontaneous<-apply(rbind(abs(minSpike),basalSpon),2,max)

adjustSpikingMatx<-rawMatx     #estimated absolute ORNs spiking rate
ORNwithoutInhiMatx<- rawMatx   #assuming no inhibition response
for(i in 1:nrow(rawMatx))
{
    adjustSpikingMatx[i,]<-rawMatx[i,]+realSpontaneous
    ORNwithoutInhiMatx[i,ORNwithoutInhiMatx[i,]<0] <- 0
    ORNwithoutInhiMatx[i,]<-ORNwithoutInhiMatx[i,] + realSpontaneous
}

#maximum spiking rate
ORNmax <- max(c(max(adjustSpikingMatx),300))  #inorder to get the value of Kij
alpha <- ORNmax/realSpontaneous - 1
weightMatx <- rawMatx
weightMatx[weightMatx>0] <- 1
weightMatx[weightMatx<0] <- -1
#weightMatxNoInhi
weightMatxNoInhi <-weightMatx
weightMatxNoInhi[weightMatxNoInhi<0] <- 0

#assuming the concentrations are 1
equilMatx <- adjustSpikingMatx
adjustSpikingMatx[adjustSpikingMatx==0] <- 0.1 # to get rid of infinity of Kij
for(i in 1:110){
    for(j in 1:24){
        if(weightMatx[i,j]>0){
            equilMatx[i,j] <- adjustSpikingMatx[i,j]*alpha[j]/(ORNmax - adjustSpikingMatx[i,j]) - 1
        }else if(weightMatx[i,j] < 0){        
            equilMatx[i,j] <- (ORNmax - adjustSpikingMatx[i,j])/(adjustSpikingMatx[i,j]*alpha[j]) -1
        }else{        
            equilMatx[i,j] <- 0
        }
    }
}

#response of mixture
#define the respone fucntion of mixture, provided by Yuhai Tu
source("mixRespFunc.R")
source("mixResp.R")
source("pcaSpectrum.R")
#source("digitMatx.R")
odorMix <- c(1)
debug(mixRespFunc)
undebug(mixRespFunc)
twoMixRes<- mixRespFunc(odorMix,ORNmax,alpha,weightMatxNoInhi,equilMatx)

numMix <- 5 #number of mixture
sampleSize <- 100000 #sampling size
conc <- 5   #concentration of odorant
#debug(mixResp)
#undebug(mixResp)
#allResp<- mixResp(numMix,sampleSize,ORNmax,alpha,weightMatx,equilMatx,conc)
allResp<- mixResp(numMix,sampleSize,ORNmax,alpha,weightMatx,equilMatx,conc,ResType="comp")
allRespNoInhi <- mixResp(numMix,sampleSize,ORNmax,alpha,weightMatxNoInhi,equilMatx,conc,ResType="comp")
#debug(pcaSpectrum)
#undebug(pcaSpectrum)
fname <- paste("continuosuas_".character(numMix),"-conc",as.character(concentration),"-method",as.character(method),sep = "")
pcaSpectrum(digiRespMatx,digiRespMatxNoInhi,numMix,scale=FALSE,fname)
#pcaSpectrum(allResp,allRespNoInhi,numMix,scale=TRUE)

#plot the histgram of response






#**************************************PN response*********************************
source("PNRatesFun.R")
Rmax<-165
n<-1.5
sigma<-12
m<-10.63

PNspiking <-allResp
PNspikingNoInhi <-allRespNoInhi
for(i in 1:dim(allResp)[1]){
    PNspiking[i,] <- PNRatesFun(allResp[i,],Rmax,n,sigam,m)
    PNspikingNoInhi[i,] <- PNRatesFun(allRespNoInhi[i,],Rmax,n,sigam,m)
}
pcaSpectrum(PNspiking,PNspikingNoInhi,numMix)

# #all the two-mixture pairs
# pairInx <- rep(0,times = 6545)
# twoMixResp <- matrix(0,6545,24)
# twoMixRespNoInhi <-matrix(0,6545,24)  #response matrix without inhibition
# count <-1
# for(i in 1:109){
#     for(j in (i+1):110){
#         mixInx <- c(i,j)
#         twoMixResp[count,]<-mixRespFunc(mixInx,ORNmax,alpha,weightMatx,equilMatx)
#         twoMixRespNoInhi[count,] <- mixRespFunc(mixInx,ORNmax,alpha,weightMatxNoInhi,equilMatx)
#         count <- count +1
#     }
# }


#priciple component analysis and comparing
#for the original system
pcaTwoMix <- princomp(twoMixResp,cor=TRUE, scores=TRUE)
loadPC<-pcaTwoMix$loadings[]
p.vaiance.explained <- pcaTwoMix$sdev^2/sum(pcaTwoMix$sdev^2)

#plot the contribution of variance of each component
par(mai = c(1.5,1.5,0.2,0.5))
barplot(100*p.vaiance.explained,width = 0.5,ylim = c(0,50),names.arg=1:24,xlab = "principle components",ylab = "explained varinace(%)",axisnames = FALSE,col="steelblue",cex.lab = 1.5,cex.axis = 1.5)

#for the system without inhibition
pcaTwoMixNoInhi <- princomp(twoMixRespNoInhi,cor=TRUE, scores=TRUE)
loadPC<-pcaTwoMixNoInhi$loadings[]
p.vaiance.explained <- pcaTwoMixNoInhi$sdev^2/sum(pcaTwoMixNoInhi$sdev^2)

#plot the contribution of variance of each component
par(mai = c(1.5,1.5,0.2,0.5))
barplot(100*p.vaiance.explained,width = 0.5,ylim = c(0,50),names.arg=1:24,xlab = "principle components",ylab = "explained varinace(%)",axisnames = FALSE,col="steelblue",cex.lab = 1.5,cex.axis = 1.5)

#3-D visulization
library("scatterplot3d")
s3d <- scatterplot3d(pcaTwoMixNoInhi$scores[,1:3],xlab = "Comp.1",ylab = "Comp.2",zlab = "Comp.3",pch = 20)

#histogram of response
hist(twoMixResp,30)
hist(twoMixRespNoInhi,30)

pdf("histORNchannelsTwoMixComp.pdf")
par(mfrow = c(2,1))
plotRange<-c(0,max(apply(twoMixResp,1,mean)))
barplot(apply(twoMixResp,2,mean),width = 0.1,ylim <-plotRange, main = "ORN",xlab = "odor",ylab = "mean spiking rate",axisnames = FALSE,col="gray",border = NULL,cex=2,cex.lab = 2,cex.axis = 2,cex.sub = 4)
barplot(apply(twoMixRespNoInhi,2,mean),width = 0.1,ylim <-plotRange, main = "PN",xlab = "odor",ylab = "mean spiking rate",axisnames = FALSE,col="gray",cex = 2,cex.lab = 2,cex.axis = 2,cex.sub=4)
dev.off()

#distance distribution
#distribution of distance pairs
par(pin<-c(12,10))
par(mai<-c(0.5,1,0.5,1))
hist(dist(twoMixResp,method = "euclidean"),50,xlab = "distance between odorants mixture",
     cex.lab = 1.5,cex.axis = 1.5,col="steelblue")
hist(dist(twoMixRespNoInhi,method = "euclidean"),50,xlab = "distance between odorants",cex.lab = 1.5,cex.axis = 1.5,col="steelblue")
