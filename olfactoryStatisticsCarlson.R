#*****************************************************************************
#this program reanalyze the Olfactory Sensory Neurons in Carlson's paper,2006
#in a recent baper writen by Charles Steven at UCSD, the author reanalyzed this 
# data set and argued that the distribution is exponentail
#writen by Shanshan Qin @ Tang Lab,Peking University
#last revised on Mar 04/01/2016
#******************************************************************************

#load the packages
library(xlsx)
library(RColorBrewer) #this package is used for color settings
library(ggplot2)
library(fitdistrplus) #statistical test
library(vcd)
library(lattice)
library(latticeExtra) #plot cumulative density functions
# library(MASS)
# library(sm)
# library(vioplot)     # violin plot of data
# library(fheatmap)    #package used to plot awsome heatmap
# library(corrgram)    #package used to plot correlation matrix
# library(corrplot)    #package that draw elegant correlation matrix
# library(rgl)         #to plot 3D figure of PCA result

# source other R scripts
source("testDistriOSNOr.R")   #plot all single OSN or Odor
source("testDistAll.R")     #plot overall response distribtuion

#change the workspace
FileFolder<-"/home/shan/GoogleDrive/olfactoryCoding"
setwd(FileFolder)

#load data file
filePath<-"/home/shan/GoogleDrive/olfactoryCoding/data/Carlson20016TableS1.xlsx"
digiFile<-"/home/shan/GoogleDrive/olfactoryCoding/data/CarslonORNdigit.xlsx"
fruitFile <- "/home/shan/GoogleDrive/olfactoryCoding/data/Carlson2006TableS2.xlsx"
#filePath<-"Carlson20016TableS1.xlsx"
rawData1<-read.xlsx(filePath,1)   #with the last line as the 
concenData <- read.xlsx(fruitFile,1) # different concentration and fruit odor mixtures
digitResponse<-read.xlsx(digiFile,1) 

#load the function for correlaogram
source("rquery.cormat.R")

#the data type of rawData1 is a list with the first column as the name
#first, calculate the distribution of firing rate regardless of ORNs and odorants
#convert the list into a vector
rawMatx<-as.matrix(subset(rawData1,select = 3:length(rawData1)))
rownames(rawMatx)<-unlist(rawData1[1]) #set the name of each row, need to convert list to vector
rawList<-as.numeric(as.list(rawMatx[-1,]))   #first, converted to a list, without spontaneous firing rate

#convert this table into absolute spiking rate
spRate <- rbind(rawMatx[111,],abs(apply(rawMatx, 2, min)))
setSponRate <- apply(spRate, 2, max)  # set the spontaneous spiking rate of each Or
resetPosiMatr <- rawMatx[-1,] + matrix(rep(setSponRate,times = 110),nrow = 110,ncol = 24,byrow = TRUE)
# another way to set the spontaneous response is just use  the rate identified
resetPosiMatr2 <- rawMatx[-1,] + matrix(rep(rawMatx[111,],times = 110),nrow = 110,ncol = 24,byrow = TRUE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#overall response data
# histogram
dev.new()
figure1<-hist(rawList,seq(min(rawList),max(rawList),20),main = NA,xlab = "Spiking Rate(Hz)",ylab = "Frequency",col = "gray",cex.lab = 1.5, cex.axis = 1.5)
box() #add box

#test if they follows exponential or lognornal distribution
figureName <- "allStatistics"
testDistAll(resetPosiMatr,"lnorm",figureName)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#all excitatory response
exciList <- rawList[rawList>0]
figureName <- "allStatExci"
testDistAll(exciList,"lnorm",figureName)

#all inhibitory response
inhiList <- abs(rawList[rawList<0])
figureName <- "allStatInhi"
testDistAll(inhiList,"exp",figureName)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# cummulative distribution of each OSN, test if it is exponential 
fileName <- "testDistriOSN"
testDistriOSNOr(resetPosiMatr,"lnorm",normalizedRate = 48,fileName)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# cumulative distribution of odorants
fileName <- "testDistriOdorant"
testDistriOSNOr(t(resetPosiMatr),"exp",normalizedRate = 48,fileName)
# pdf("testDistriOdorant.pdf",width = 8,height = 6)
# par(mfrow=c(2,2),mai=c(0.7,0.8,0.3,0.5))
# 
# # original data
# plot.ecdf(rawMatx[1,],do.points=FALSE,xlim = c(-50,290),main = "",xlab = "spiking rate (HZ)", ylab = "cdf",cex.lab = 1.5,cex.axis = 1.5)
# for (i0 in 2:dim(rawMatx)[1]){
#   temp <- ecdf(rawMatx[i0,])
#   plot(temp,verticals=TRUE, do.points=FALSE,add=TRUE)
# }
# 
# #normalized to 100 Hz 
# normVal <- 50  # normalized average sipking rate
# averageSpikesOdor <- apply(resetPosiMatr, 1, mean)
# afterScaleMatrOdor <- resetPosiMatr*matrix(rep(normVal/averageSpikesOdor,24),nrow = 110, ncol = 24, byrow = FALSE)  #normalized to average spking rate as 100 Hz
# 
# 
# index <- seq(1,300,1)
# allCdf <- matrix(0,nrow = length(index),ncol = dim(afterScaleMatrOdor)[1])
# 
# fit1 <- ecdf(afterScaleMatrOdor[1,])
# allCdf[,1] <- fit1(index)
# plot(fit1,do.points=FALSE,xlim = c(0,300),main ="",ylab = "cdf",xlab = "spiking rate (Hz)",cex.lab = 1.5,cex.axis = 1.5)
# for (i0 in 2:dim(afterScaleMatrOdor)[1]){
#   temp <- ecdf(afterScaleMatrOdor[i0,])
#   plot(temp,verticals=TRUE, do.points=FALSE,add=TRUE)
#   allCdf[,i0] <- temp(index)
# }
# 
# # generated from exponential distributions
# plot(ecdf(rexp(24,rate = 1/normVal)),verticals=TRUE, do.points=FALSE,main = "",xlim = c(0,300),xlab = "spiking rate (Hz)",ylab = "cdf",cex.lab = 1.5,cex.axis = 1.5)
# for (i0 in 2:dim(afterScaleMatrOdor)[1]){
#   plot(ecdf(rexp(24,rate = 1/normVal)),verticals=TRUE, do.points=FALSE,add=TRUE)
# }
# 
# #plot the average cdf
# plot(index,apply(allCdf, 1,mean),main = "",xlab = "spiking rate (Hz)",ylab = "cdf",cex.lab = 1.5,cex.axis = 1.5)
# curve(pexp(x, rate = 1/normVal), col = "red", lwd =3,add = TRUE) #overlay an exponential distribution
# dev.off()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# response under different concentration and fruit odor mixture
# first, monomolecular at different concentration
odorNames <- concenData[,1]
dataMat <- as.matrix(concenData[,3:26])
allConPure <- concenData[1:40,2]
sp <- abs(apply(dataMat, 2, min))   # spontaneous rate
dataReset <- dataMat + matrix(rep(sp,dim(dataMat)[1]), nrow = dim(dataMat)[1],ncol = dim(dataMat)[2],byrow = T)
pureOdor <- dataReset[1:40,]  #pure odorants
fruit <- dataReset[41:76,]    #fruit odor

#pure odor under diffferent concentrations
saveFile <- "disPureOdorCon"
testDistAll(pureOdor,"lnorm",saveFile)

#fruit odor mixture
saveFile <- "disFruitOdorCon"
testDistAll(fruit,"lnorm",saveFile)

#all pure odor and fruit under different concentrations
saveFile <- "disAllOdorCon"
testDistAll(dataReset,"lnorm",saveFile)

con <- 10^(c(-8,-6,-4,-2))
# aggregate(pureOdor,by = list(pureOdor$dilution),mean)
meanStatisc <- matrix(0,nrow = length(con),2)  #mean and std
for (i0 in 1:4){
   srow <- allConPure == 2*i0
   meanStatisc[i0,]<- c(mean(pureOdor[srow,]),sd(pureOdor[srow,]))
}
plot(sort(con,decreasing = T),meanStatisc[,1], log = "xy")


#average spiking rate under different concentration

pureOdor <- grep("-[0-9]$",odorNames)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#function that plot all the basic statistical figures in the paper

plotBasic<-function(rawMatx)
{
    #x is the data frame that contain the original data
    rawList<-as.numeric(as.list(rawMatx))   #first, converted to a list, then convert data type
    dev.new()
    figure1<-hist(rawList,seq(min(rawList),max(rawList),20),xlab = "Spiking Rate(Hz)",ylab = "Frequency",col = "gray")
    
    maxSpike<-apply(rawMatx,2,max)  #the largest spiking rate of each ORN
    minSpike<-apply(rawMatx,2,min)  #strongest inhibition of each ORN
    #open a graphicd device and store the file
    pdf("allTheORNsResponse.pdf",width = 14,height = 14) #set the size of the graphics
    op<-par(mfrow = c(5,5))
    for(i in 1:length(rawMatx))
    {       
        barplot(rawMatx[order(rawMatx[,i]),i],xlab = NULL,ylim = c(minSpike,maxSpike))
        par(pin=c(2,2),mai=c(0.3,0.3,0.3,0.3)) #set the figure size
    }
    par(op)    #At the end of ploting, recet to previous settings
    dev.off()  #close the graphic device
    
    #strong and inhibitory response
    #the respones that are larger than 100 spikes per second
    strongActThreshold <-100
    inhiTheshold <- -10
    strongSp<-rep(0,length(rawMatx)) #initialization
    inhibiSp<-rep(0,length(rawMatx))
    for(i in 1:length(rawMatx))
    {
        strongSp[i]<-length(rawMatx[rawMatx[,i]>strongActThreshold,i])/nrow(rawMatx)
        inhibiSp[i]<-length(rawMatx[rawMatx[,i]< 0.5*min(rawMatx[,i]),i])/nrow(rawMatx)
    }
    dev.new()
    par(mfrow = c(2,1))
    with(actInhiSp,{
        barplot(strongSp[order(strongSp)],names.arg = ORNnames[order(strongSp)],ylab = "Percentage of strong activation",space = 0.5,col=colorRampPalette(brewer.pal(9,"Blues"))(24),axes = TRUE)
        barplot(inhibiSp[order(inhibiSp)],names.arg = ORNnames[order(inhibiSp)],ylab = "Percentage of inhibition",space = 0.5,col=colorRampPalette(brewer.pal(9,"Blues"))(24))
    })
    
    #percentage of odorant that elicits strong activation and inhibition response
    pdf("odorantsRecpActInhi.pdf",width = 12,height = 9)
    par(mfrow = c(2,3))
    percentActInbiAll<-matrix(0,100,5) #this data is used to produced fig 2C
    for(i in 1:nrow(rawMatx))
    {
        percentActInbiAll[i,1]<-length(rawMatx[i,rawMatx[i,]>50])
        percentActInbiAll[i,2]<-length(rawMatx[i,rawMatx[i,]>100])
        percentActInbiAll[i,3]<-length(rawMatx[i,rawMatx[i,]>150])
        percentActInbiAll[i,4]<-length(rawMatx[i,rawMatx[i,]>200])
        percentActInbiAll[i,5]<-length(rawMatx[i,rawMatx[i,]<0])              
    }
    for(i in 1:5)
    {
        barplot(percentActInbiAll[order(percentActInbiAll[,i]),i],names.arg = NULL,ylim = c(0,18),main = "excitation to >= 50spikes/s",ylab = "Receptors")
        par(pin=c(3.5,3),mai=c(0.5,0.5,0.3,0.3)) #set the figure size
    }
    dev.off()
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function test if the response follows a certain distribution




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxSpike<-apply(rawMatx,2,max)  #the largest spiking rate of each ORN
minSpike<-apply(rawMatx,2,min)  #strongest inhibition of each ORN

#we estimate the real spontaneous spiking rate of a ORN as the max between inhibitory spiking rate and basal spiking rate
realSpontaneous<-apply(cbind(abs(minSpike),rawMatx[nrow(rawMatx),]),1,max)
adjustSpikingMatx<-rawMatx[-nrow(rawMatx),]     #estimated absolute ORNs spiking rate
ORNwithoutInhiMatx<- rawMatx[-nrow(rawMatx),]   #assuming no inhibition response
for(i in 1:(nrow(rawMatx)-1))
{
    adjustSpikingMatx[i,]<-rawMatx[i,]+realSpontaneous
    ORNwithoutInhiMatx[i,ORNwithoutInhiMatx[i,]<0] <- 0
    ORNwithoutInhiMatx[i,]<-ORNwithoutInhiMatx[i,] + realSpontaneous
}


#the respones that are larger than 100 spikes per second
strongActThreshold <-100
inhiTheshold <- -10
strongSp<-rep(0,length(rawMatx)) #initialization
inhibiSp<-rep(0,length(rawMatx))
for(i in 1:length(rawMatx))
{
    strongSp[i]<-length(rawMatx[rawMatx[,i]>strongActThreshold,i])/nrow(rawMatx)
    inhibiSp[i]<-length(rawMatx[rawMatx[,i]< 0.5*min(rawMatx[,i]),i])/nrow(rawMatx)
}
ORNnames<-colnames(rawMatx)
actInhiSp<-data.frame(ORNnames,strongSp,inhibiSp)


#correlation graphics
cg1<-rquery.cormat(t(rawMatx[1:110,]),type = "full")
rquery.cormat(t(adjustSpikingMatx),type = "full")

#use spearman correlation
corrCoefMatx<-cor(rawMatx,method = "spearman")
corrgram(corrCoefMatx,order = "TRUE",lower.panel = panel.shade,upper.panel = panel.pie,text.panel = panel.txt,main="Corrlogram of spiking rate among ORNs")

#function that compares the correlation change
#debug(corrComp)
corrComp<-function(X,Y,
                   plot="TRUE")
{
    #input are two matrix has the same dimension matirx
    if(dim(X)[1] != dim(Y)[1]){
        stop("The lengths of X and Y are different!")}
    else{
        #correlation marix of X and Y
        c1<-cor(X)
        c2<-cor(Y)
        #because the matrix is symmetrical, only the upper or lower triangle is usefull
        cTri1<-c1[lower.tri(c1)]
        cTri2 <-c2[lower.tri(c2)]
        totChg<-sum(abs(cTri2)) - sum(abs(cTri1))  #totally change of correlation
        
        #plot the correlation value of each pair
        par(mai=c(2,2,1,2))
        plotlim<-range(c(abs(cTri1),abs(cTri2)))
        fig1<- plot(abs(cTri1),abs(cTri2),xlim = plotlim, ylim = plotlim,pch = 20,main = "Comparing the change of correlation of each pair",xlab = "orignal",ylab = "without inhibition",cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5)
        abline(a = 0,b=1,lwd = 2,col = "gray")  #add a line with 
        final<-list(totChg,fig1)
        return(final)
    }
}
f1<-corrComp(t(adjustSpikingMatx),t(ORNwithoutInhiMatx))
#strongest inhibition response of all the ORNs
maxInhibition<-lapply(rawMatx,1,min)

#cluster based on the response pattern
#heatmap of all the response
library("gplots")
pdf("heatMapORNnoinhi.pdf",width = 10,height = 16)
#heatmap.2(ORNwithoutInhiMatx, col=redblue(100), scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
fheatmap(ORNwithoutInhiMatx, cluster_rows = T,mat_color =colorRampPalette(c("navy","white","firebrick3"))(50),cell_border = FALSE, fontsize = 4)
#pheatmap(ORNwithoutInhiMatx, color =colorRampPalette(c("navy","white","firebrick3"))(50),fontsize = 12,fontsize_row = 13)
dev.off()
#heatmap.2(as.matrix(rawMatx[1:110,]), col=redblue(100), scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)

#the original data
pdf("heatMapORNrawAdjust.pdf",width = 10,height = 16)
fheatmap(adjustSpikingMatx, cluster_rows = T,cluster_distance_cols = "manhattan",mat_color =colorRampPalette(c("navy","white","firebrick3"))(50),cell_border = FALSE, fontsize = 4)
dev.off()


cluster<-function(X, #matrix contain observes
                  disMeth = c("euclidean","manhattan","cosine"), #distance method
                  clusMeth = c("average","minimum")#clustering method
                  )
{
    hcEuc<-hclust(dist(X,method = disMeth),method = clusMeth)
    plot(hcEuc,hang=-1,labels = rownames(rawMatx))
    return(hcEuc)   #return figure handle
}
#distance can be Euclidean, Cosine,Manhatten or others
#Euclean
OrnRawc1<-cluster(rawMatx)
#Manhattan
OrnRawc2<-cluster(rawMatx,disMeth="manhattan")


hcEuc<-hclust(dist(rawMatx,method = "euclidean"),method = "average")
plot(hcEuc,hang=-1,labels = rownames(rawMatx))
dev.off()
#Manhattan
dev.new()
hcEuc<-hclust(dist(rawMatx,method = "manhattan"),method = "average")
plot(hcEuc,hang=-1,labels = rownames(rawMatx))  #a minus hang means that label are 
dev.off()


#****************************Projection Neurons Response***************************
#Projection nerons and tranfrorm of previous pattern
#define the transform function,which is an empirical funcition from Wilson's 2010 paper
#PN=Rmax*ORN^n/(ORN^n + sigma^n + (m*sumofORNs/190)^n)
PNRates<-function(ORNrates,Rmax,n,sigam,m){
    #ORNrates[ORNrates<0] <- 0 
    sumOfRats<-sum(ORNrates)
    result<-Rmax*ORNrates^n/(ORNrates^n + sigma^n + (m*sumOfRats/190)^n)
    return(result)
}
Rmax<-165
n<-1.5
sigma<-12
m<-10.63

#initialization of PN response data frame
#PNspiking<-data.frame(matrix(NA, ncol = 24, nrow = nrow(rawMatx))-1) 
#names(PNspiking)<-colnames(rawMatx)
#row.names(PNspiking)<-subset(rownames(rawData1),select = 1:(nrow(rawData1)-1))
PNspiking<-rawMatx[-nrow(rawMatx),]  #the last line is deleted due to spontaneous response
PNspikingNoInhi <- rawMatx[-nrow(rawMatx),] #the spiking rate if no inhibition of ORNs
for(i in 1:nrow(PNspiking)){
    PNspiking[i,]<-PNRates(adjustSpikingMatx[i,],Rmax,n,sigam,m)
    PNspikingNoInhi[i,]<- PNRates(ORNwithoutInhiMatx[i,],Rmax,n,sigam,m)
}

#distribution of spiking rate in PNs, does it follow exponential distribution?
dev.new()
hist(as.matrix(PNspiking),30,xlab = "PN spiking rate",col = "steelblue",cex=1.5)
hist(as.matrix(PNspikingNoInhi),30,xlab = "PN spiking rate",col = "steelblue",cex=1.5)

#calculate the gain of tranformation, from ORN to PN
maxORN<-max(adjustSpikingMatx)  #maximum spiking rate of ORN
gainRaw<-(as.vector(PNspiking)/Rmax)/(as.vector(adjustSpikingMatx)/maxORN)
gainNoInhi<-(as.vector(PNspikingNoInhi)/Rmax)/(as.vector(ORNwithoutInhiMatx)/maxORN)

library(scales)
hist(gainRaw,25,xlim=c(0,6),col='skyblue',border=F,xlab = "relatively gain")
hist(gainNoInhi,25,add=T,col=scales::alpha('red',.5),border=F)
plot(h1,col=rgb(0,0,1,1/4),xlim = c(0,6))
plot(h2,col=rgb(0,0,1,1/4),xlim = c(0,6),add=T)
plot(as.vector(adjustSpikingMatx),as.vector(PNspiking),xlab = "ORN spiking rate", ylab = "PN spiking rate")
plot(as.vector(ORNwithoutInhiMatx),as.vector(PNspikingNoInhi),xlab = "ORN spiking rate", ylab = "PN spiking rate")


#mean spiking rate of ORNs and PNs
pdf("compareTranformationAverResp.pdf",width = 13,height =10)
par(mfrow = c(2,1))
plotRange<-c(0,max(apply(adjustSpikingMatx,1,mean)))
barplot(apply(adjustSpikingMatx,1,mean),width = 0.1,ylim <-plotRange, main = "ORN",xlab = "odor",ylab = "mean spiking rate",axisnames = FALSE,col="gray",border = NULL,cex=2,cex.lab = 2,cex.axis = 2,cex.sub = 4,mai = c(0.2,2,0.2,0.5))
barplot(apply(as.matrix(PNspiking),1,mean),width = 0.1,ylim <-plotRange, main = "PN",xlab = "odor",ylab = "mean spiking rate",axisnames = FALSE,col="gray",cex = 2,cex.lab = 2,cex.axis = 2,cex.sub=4,mai = c(0.2,2,0.2,0.5))
dev.off()

pdf("PNresponseAdjust.pdf",width = 12,height = 20)
heatmap.2(PNspiking, col=redblue(100), scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()

#test whether it follows exponential distribution
#cumulative emperical distribution
FitExp <- fitdist(as.vector(PNspiking),"exp")
plot(FitExp)
PNcdf<-ecdf(as.vector(PNspiking))
plot(PNcdf,do.points=FALSE, verticals=TRUE)

#**************************************************************************

#PCA of data
library(psych)
fa.parallel(adjustSpikingMatx,fa = "pc",n.iter = 100,show.legend = FALSE)
pcORNadj<-principal(adjustSpikingMatx,nfactors = 24,scores = TRUE)
pcPNadj<-principal(PNspiking,nfactors = 24,scores = TRUE)

PNcorrCoefMatx<-cor(PNspiking)
ORNcorrCoefMatxAdj<-cor(adjustSpikingMatx)

#corrlogram
corrgram(PNcorrCoefMatx,order = "TRUE",lower.panel = panel.shade,upper.panel = panel.pie,text.panel = panel.txt,main="Corrlogram of spiking rate among PNs")
adM<-rquery.cormat(PNspikingNoInhi,type = "full")

rquery.cormat(t(PNspikingNoInhi),type = "full")
rquery.cormat(PNspiking,type = "full")
rquery.cormat(t(PNspiking),type = "full")
corrgram(ORNcorrCoefMatxAdj,order = "TRUE",lower.panel = panel.shade,upper.panel = panel.pie,text.panel = panel.txt,main="Corrlogram of spiking rate among ORNs")

#compare the change of correlation matrix
corrComp(PNspiking,PNspikingNoInhi)
#********************without inhibiton ************************************
#we assume that the system can only have basal level response or activation
#we replace all the inhibitory response with the basal level activity
#**************************************************************************
PCAvar3dplot<-function(x,            #matrix or data frame
                       cateName,     #name for each observes
                       figure=TRUE
    )
{
    #library("rgl")           #for 3-d plot interactive
    library("scatterplot3d") #for 3-d plot
    
    #combine x and name into a new data frame
    if(is.matrix(x)){
        y<-cbind(as.data.frame(x),cateName)
    }
    else if(is.data.frame(x)){
        y<-cbind(x,cateName)
    }
    pc<-princomp(y[,1:(ncol(y)-1)],cor=TRUE, scores=TRUE)
    #load the components  ORN
    loadPC<-pc$loadings[]
    chemicalCat <- factor(y$cateName) #name of observe values
    cols <-as.double(chemicalCat)     #color for different categories
    vaiance.explained <- pc$sdev^2/sum(pc$sdev^2)  #variance explained
    #plot a histgram of the explained variance
    par(mai = c(1.5,1.5,0.2,0.5))
    barplot(100*vaiance.explained,width = 0.5,ylim = c(0,50),xlab = "principle components",ylab = "explained varinace(%)",axisnames = FALSE,col="steelblue",cex.lab = 1.5,cex.axis = 1.5)
    dev.new()
    s3d <- scatterplot3d(pc$scores[,1:3],xlab = "Comp.1",ylab = "Comp.2",zlab = "Comp.3",color = cols,pch=20)
    #legend('right', cex=.5,legend = levels(chemicalCat), fill = 1:nlevels(chemicalCat), merge = F, bty = 'n')
    
}

#calculate the the survival function

PCAvar3dplot(ORNwithoutInhiMatx,categoryName)

#pcORNnoInhi<-principal(ORNwithoutInhiMatx,nfactors = 4,scores = TRUE)
categoryName<-rawData1$category[-nrow(rawData1)]

newORNwithoutInhi<- cbind(as.data.frame(ORNwithoutInhiMatx),categoryName)
newPNnoInhiSpi<- cbind(as.data.frame(ORNwithoutInhiMatx),categoryName)
newRawMatx<-cbind(as.data.frame(rawMatx[1:110,]),categoryName)

pcORNnoInhi<-princomp(newORNwithoutInhi[,1:24],cor=TRUE, scores=TRUE)
pcPNnoInhi <-princomp(PNspikingNoInhi[,1:24],cor=TRUE, scores=TRUE)
pcORNraw <- princomp(newRawMatx[,1:24],cor=TRUE, scores=TRUE)

#load the components  ORN
loadPC<-pcORNnoInhi$loadings[]
chemicalCat <- factor(newORNwithoutInhi$categoryName)
cols = as.double(chemicalCat)

#variance explained by each conponent
p.vaiance.explained <- pcORNnoInhi$sdev^2/sum(pcORNnoInhi$sdev^2)

#plot a histgram of the explained variance
par(mai = c(1.5,1.5,0.2,0.5))
barplot(100*p.vaiance.explained,width = 0.5,ylim = c(0,50),names.arg=1:24,xlab = "principle components",ylab = "explained varinace(%)",axisnames = FALSE,col="steelblue",cex.lab = 1.5,cex.axis = 1.5)

#load the components  PN
loadPC<-pcPNnoInhi$loadings[]
chemicalCat <- factor(newPNnoInhiSpi$categoryName)
cols = as.double(chemicalCat)

#variance explained by each conponent
p.vaiance.explained <- pcPNnoInhi$sdev^2/sum(pcPNnoInhi$sdev^2)

#plot a histgram of the explained variance
par(mai = c(1.5,1.5,0.2,0.5))
barplot(100*p.vaiance.explained,width = 0.5,ylim = c(0,50),names.arg=1:24,xlab = "principle components",ylab = "explained varinace(%)",axisnames = FALSE,col="steelblue",cex.lab = 1.5,cex.axis = 1.5)

#show the first three components in 3-D
library('scatterplot3d')
s3d <- scatterplot3d(pcORNraw$scores[,1:3],xlab = "Comp.1",ylab = "Comp.2",zlab = "Comp.3",color = cols,pch=20)
legend('right', cex=.5,legend = levels(chemicalCat), fill = 1:nlevels(chemicalCat), merge = F, bty = 'n') 

plot3d(pcORNnoInhi$scores[,1:3],col=cols)
legend('topright', cex=.8,  legend = levels(sectors), fill = 1:nlevels(sectors), merge = F, bty = 'n') 
plot3d(pcORNnoInhi$scores[,1:3],box = TRUE,axis = TRUE,col = newORNwithoutInhi$categoryName)

ORNcorrMatxNoInhi<-cor(ORNwithoutInhiMatx)
PNcorrMatxNoInhi<-cor(PNspikingNoInhi)
rquery.cormat(ORNwithoutInhiMatx[1:110,],type = "full")
corrgram(ORNcorrMatxNoInhi,order = "TRUE",lower.panel = panel.shade,upper.panel = panel.pie,text.panel = panel.txt,main="Corrlogram of spiking rate among PNs")
corrgram(PNcorrMatxNoInhi,order = "TRUE",lower.panel = panel.shade,upper.panel = panel.pie,text.panel = panel.txt,main="Corrlogram of spiking rate among PNs")

#plot the spiking rate
pdf("histogramOfORNsSpikingRateNoInhi.pdf",width = 8,height = 6)
hist(ORNwithoutInhiMatx,30,xlab = "Spiking Rate(Hz)",ylab = "Frequency",col = "gray",cex.lab = 1.5,cex.axis=1.5,lwd = 1.5)
dev.off()

pdf("histogramOfPNsSpikingRateNoInhi.pdf",width = 8,height = 6)
hist(as.matrix(PNspikingNoInhi),30,xlab = "Spiking Rate(Hz)",ylab = "Frequency",col = "gray",cex.lab = 1.5,cex.axis=1.5,lwd = 1.5)
dev.off()

#comparing the average spiking rate of ORNs and PNs with the assumption that no inhibitory ORN response
pdf("compareORN-PN-AverRespNoInhi.pdf",width = 13,height =10)
par(mfrow = c(2,1))
plotRange<-c(0,max(apply(ORNwithoutInhiMatx,1,mean)))
barplot(apply(ORNwithoutInhiMatx,1,mean),width = 0.1,ylim <-plotRange, main = "ORN",xlab = "odor",ylab = "mean spiking rate",axisnames = FALSE,col="gray",border = NULL,cex=2,cex.lab = 2,cex.axis = 2,cex.sub = 4,mai = c(0.2,2,0.2,0.5))
barplot(apply(as.matrix(PNspikingNoInhi),1,mean),width = 0.1,ylim <-plotRange, main = "PN",xlab = "odor",ylab = "mean spiking rate",axisnames = FALSE,col="gray",cex = 2,cex.lab = 2,cex.axis = 2,cex.sub=4,mai = c(0.2,2,0.2,0.5))
dev.off()

#distance among different odor, defined in the 24-ORNs space
#Euclidean for ORN response
dev.new()
hcEucORNnoInhi<-hclust(dist(ORNwithoutInhiMatx,method = "euclidean"),method = "average")
plot(hcEucORNnoInhi,hang=-1,labels = rownames(ORNwithoutInhiMatx))
dev.off()

#Euclidean for PN response
dev.new()
hcEucPNnoInhi<-hclust(dist(PNspikingNoInhi,method = "euclidean"),method = "average")
plot(hcEucPNnoInhi,hang=-1,labels = rownames(PNspikingNoInhi))
dev.off()

#Manhattan distance for ORN response
dev.new()
hcManhORN<-hclust(dist(ORNwithoutInhiMatx,method = "manhattan"),method = "average")
plot(hcManhORN,hang=-1,labels = rownames(ORNwithoutInhiMatx))  #a minus hang means that label are 
dev.off()

#Manhattan distance for PN response
dev.new()
hcManhPN<-hclust(dist(PNspikingNoInhi,method = "manhattan"),method = "average")
plot(hcManhPN,hang=-1,labels = rownames(PNspikingNoInhi))  #a minus hang means that label are 
dev.off()

#Cosine distance/simiarity
cosineDist <- function(x){
    cosDis<-as.matrix(as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2)))),diag = T,upper = T))
    allAngle<-acos(pmax(pmin(x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2)))),1.0),-1.0))
    final<- list(cosDis,allAngle)
    return(final)
}
cosPN <- cosineDist(PNspiking)
cosPNnoInhi <- cosineDist(PNspikingNoInhi)
#rquery.cormat(cosPNnoInhi,graphType = "heatmap")
#rquery.cormat(cosPN,graphType = "heatmap")

#plot the change
plotlim<-range(c(cosPN[[1]],cosPNnoInhi[[1]]))
par(mai = c(1,1,1,1))
plot(cosPN[[1]][lower.tri(cosPN[[1]])],cosPNnoInhi[[1]][lower.tri(cosPNnoInhi[[1]])],pch = 20,main = "Comparing the change of correlation of each pair",xlab = "orignal",ylab = "without inhibition",cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5)
abline(a=0,b=1,lwd=2,col="gray")

#histgram of angles
limXlab<-range(c(cosPN[[2]]),cosPNnoInhi[[2]])
hist(cosPN[[2]],25,xlim=limXlab,col='skyblue',border=F,xlab = "angle between pattern")
hist(cosPNnoInhi[[2]],25,add=T,col=scales::alpha('red',.5),border=F)

#plot the corrematrix
corrplot(cosPNnoInhi[[1]],method = "circle",order="hclust")


dev.new()
hcManhORN<-hclust(cosPNnoInhi,method = "average")
plot(hcManhORN,hang=-1,labels = rownames(cosPNnoInhi))  #a minus hang means that label are 
dev.off()

#distribution of distance pairs
par(pin<-c(12,10))
par(mai<-c(0.5,1,0.5,1))
hist(dist(PNspikingNoInhi,method = "manhattan"),50,xlab = "distance between odorants",
     cex.lab = 1.5,cex.axis = 1.5,col="steelblue")
hist(dist(ORNwithoutInhiMatx,method = "manhattan"),50,xlab = "distance between odorants",cex.lab = 1.5,cex.axis = 1.5,col="steelblue")
hist(dist(PNspikingNoInhi,method = "euclidean"),50,xlab = "distance between odorants",cex.lab = 1.5,cex.axis = 1.5,col="steelblue")
hist(dist(rawMatx[1:110,],method = "euclidean"),50,xlab = "distance between odorants",cex.lab = 1.5,cex.axis = 1.5,col="steelblue")

#*******************************digital response********************************************
#correlogram 
digitOri<-as.matrix(digitResponse[,3:26])
rownames(digitOri) <- digitResponse$odorantName
digitNoInhi <- digitOri
digitNoInhi[digitNoInhi==-1] <- -0.01 #can't be set as 0
rquery.cormat(digitOri,type = "full")
rquery.cormat(digitNoInhi,type = "full")

#distance matrix
digitOriDist<-dist(digitOri)
digitNoInhiDist<-dist(digitNoInhi)

plotlim <- range(c(digitOriDist,digitNoInhiDist))
par(mai = c(2,2,1,0.5))
plot(digitOriDist,digitNoInhiDist,xlim = plotlim, ylim = plotlim,pch = 20,main = "Comparing the change of correlation of each pair",xlab = "orignal",ylab = "without inhibition",cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5)
abline(a=0,b=1,lwd=2,col="gray")

#cosine angle
digCosOri <- cosineDist(digitOri)
digCosNoIn<- cosineDist(digitNoInhi)

#plot cosine similarity
#plotlim<-range(c(digCosOri[[1]],digCosNoIn[[1]]))
plotlim <- c(0,2)
par(mai = c(2,2,1,1))
plot(digCosOri[[1]][lower.tri(digCosOri[[1]])],digCosNoIn[[1]][lower.tri(digCosNoIn[[1]])],xlim = plotlim, ylim = plotlim,pch = 20,main = "Comparing the change of correlation of each pair",xlab = "orignal",ylab = "without inhibition",cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5)
abline(a=0,b=1,lwd=2,col="gray")

#histgram of angles
limXlab<-range(c(digCosOri[[2]]),digCosNoIn[[2]])
limXlab <- c(0,2.5)
hist(digCosOri[[2]],25,xlim=limXlab,col='skyblue',border=F,xlab = "angle between pattern")
hist(digCosNoIn[[2]],25,add=T,col=scales::alpha('red',.5),border=F)


#*******************************attraction-inversion*****************************
#this part try to figure out whether there is a relation between odorant-evoked pattern and 
#behavior difference

#load file
pattBehFil<- "Knaden2012S1new.xlsx"
pBehaMatx<-read.xlsx(pattBehFil,1)

#pattern distances at ORN level
pattDist<- dist(pBehaMatx[,4:length(pBehaMatx)])
behaDist<- dist(pBehaMatx[,2])

#plot
#plotlim<-range(c(cosPN[[1]],cosPNnoInhi[[1]]))
par(mai = c(1.5,1.5,1,1))
plot(behaDist,pattDist,pch = 20,main = NULL,xlab = "attractive index distance",ylab = "pattern distance",cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5)
abline(a=0,b=1,lwd=2,col="gray")

#PN pattern
pbMatx<- as.matrix(pBehaMatx[,4:length(pBehaMatx)])
rownames(pbMatx) <- pBehaMatx$Odorant
minSpike<-apply(pbMatx,2,min)  #strongest inhibition of each ORN
realSpontaneous<-apply(cbind(abs(minSpike),pbMatx[nrow(pbMatx),]),1,max)
adjBPMatx<-pbMatx
bpPNspiking <- pbMatx
for(i in 1:(nrow(adjBPMatx)-1))
{
    adjBPMatx[i,]<-adjBPMatx[i,]+realSpontaneous
    bpPNspiking[i,] <- PNRates(adjBPMatx[i,],Rmax,n,sigam,m)
}
bpPNdist<- dist(bpPNspiking)
par(mai = c(1.5,1.5,1,1))
plot(behaDist,bpPNdist,pch = 20,main = NULL,xlab = "attractive index distance",ylab = "pattern distance",cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5)
abline(a=0,b=1,lwd=2,col="gray")

#cosine distance
cosBPPNspikeDist <- cosineDist(bpPNspiking)
par(mai = c(1.5,1.5,1,1))
plot(behaDist,cosBPPNspikeDist[[1]][lower.tri(cosBPPNspikeDist[[1]])],pch = 20,main = NULL,xlab = "attractive index distance",ylab = "pattern distance",cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.5)

#****************************Related to Prof.Luo experiments***************************
# in this part, we try to relate the distance between spiking pattern with Prof. Luo's expriments
# the test ORN is 85a, the WT, 85a-only and lose of function  response are from Carlson's 2006 paper
# while the corresponding PN response are different, for the 85a-only, no lateral inhibition, while
# both the WT and lose of function mutant use the transformation function above