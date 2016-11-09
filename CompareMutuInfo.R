#this script compares the mutual information between any two OSN
#we only sample the odors that stimulates either of the selected two OSN

library(xlsx)


source("mutualInfo.R")
source("setBindAffiInhiMatx.R")

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

#weight matrix
weightMatx <- newMatx
weightMatx[weightMatx>0] <- 1
weightMatx[weightMatx<0] <- -1


DigitSp<- c(15,30,50,100,150,200)

equilMatx <- setBindAffiInhiMatx(newMatx,DigitSp,ORNmax,alpha,beta=1,method=1)
#browser()
mutualInfo(p0=0.1,weightMatx, equilMatx,ORNmax,alpha,beta = 1, concentration = 1, method = 1,interType="multi",numSamp=10000)