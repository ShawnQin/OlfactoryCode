#*************************************************************
#GROGRAM: olfMixDigitDis.R
#DESCRIPTION:
#this program calculate the ORN response patterns of odor mixtures
#for the sytem with inhibitory response and without inhibitoyr response
#the integration rule at ORN is proposed by Prof. Tu
#it differs from olfMixDsic by
# 1) same baseline of all ORN
# 2) responses are digitalized from -1 to 5
#writen by qinsspku@163.com
#last revised on 10/07/2016
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
source("plotPCAEntropy.R")
source("histCount.R")

#load data file
digiFile<-"data/CarslonORNdigit.xlsx"
#rawData1<-read.xlsx(filePath,1)   #with the last line as the 
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

#weight matrix
weightMatx <- newMatx
weightMatx[weightMatx>0] <- 1
weightMatx[weightMatx<0] <- -1

#weightMatxNoInhi
weightMatxNoInhi <-weightMatx
weightMatxNoInhi[weightMatxNoInhi<0] <- 0

#assuming the concentrations are 1
DigitSp<- c(15,30,50,100,150,200)


method <- 1
beta <- 0.1
# debug(setBindAffiInhiMatx)
# equilMatx <- setBindAffiInhiMatx(newMatx,DigitSp,ORNmax,alpha,beta,method,ResType="multi")
# equilMatxNoInhi <-equilMatx
# equilMatxNoInhi[newMatx<0] <- 0

# #debug(compAll)
balanceMethod <-1 #method used to balance
cCon <- 1 #concnetration
Rtype <- "comp"
compNumMix <- c(5,20,35)
#compNumMix <- 5
# debugonce(compAll)
# startTime <- proc.time()
pcaEntr <-compAll(method = balanceMethod,compNumMix,compCon = cCon,ResType = Rtype,plotOrNot = FALSE,spontaneous,ORNmax,ORNmax)
# proc.time() - startTime
#pcacmp<-subset(pcaEntr$PCAresult,numMix %in% c(5,20,35))
# browser()
pcacmp <- pcaEntr$PCAresult
sumEntropy<-pcaEntr$sumEnt
effectivDim <- aggregate(pcacmp$eigenvalue,by=list(system = pcacmp$system,numMix=pcacmp$numMix),function(x) sum(x>=1))
#
plotPCAEntropy(compNumMix,pcaEntr$sumEnt,pcacmp,balanceMethod,cCon,Rtype)

# myColor <- brewer.pal(11,"Spectral")
# yRange <-c(20,90)
# 
# par(bty="l")
# plot(compNumMix,sumEntropy[,1],type = "b",pch=16,lty=1,lwd=2,cex=2,col=myColor[1],ylim = yRange,xlab = "Number of odor",ylab = "Total Entropy (bit)")
# lines(compNumMix,sumEntropy[,2],type = "b",pch=16,lty=1,lwd=2,cex=2,col=myColor[10])
# # dev.off()
# sumEntro<data.frame(sys=factor(rep(c("Ori","No"),each=dim())))
# peig <- ggplot(subset(pcacmp,pcacmp$system=="Ori"),aes(x=component, y=eigenvalus, group=numMix)) + geom_line(data=subset(pcacmp,pcacmp$system=="Ori"),linetype="solid",color=myColor[1])+geom_point(data=subset(pcacmp,pcacmp$system=="Ori"),aes(shape=numMix),color=myColor[1],size=4)+theme(legend.position="top")



#new method 
source("allSingOSNresp.R")
# debug(allSingOSNresp)
interactionType = "comp"  #  multiple binding sites and competitive binding
sparsePo = (1:12)*5/107
 # sparsePo = c(5,20)/107
 sumEntropy = matrix(0,length(sparsePo),2)
 for(i in 1:length(sparsePo)){
   allEntropy <-allSingOSNresp(p0=sparsePo[i],newMatx,DigitSp,ORNmax,alpha,beta = 1,concentration = 1,method = 1,interType=interactionType)
   sumEntropy[i,]<-c(sum(allEntropy$entrOri),sum(allEntropy$entroNoInhi))
 }
 plotPCAEntropy((1:12)*5,sumEntropy)

# #mutual information among differnt OSNs
# source("mutualInfo.R")
# # debug(mutualInfo)
# mutualInfo(p0=0.2,weightMatx, equilMatx,ORNmax,alpha,beta = 1, concentration = 1, method = 1,interType="multi",numSamp=100000)


