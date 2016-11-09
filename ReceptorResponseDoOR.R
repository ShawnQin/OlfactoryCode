#This program extract olfactory response data from DoOR.2 data base

library(DoOR.functions)
library(DoOR.data)
library(ggplot2)
loadData()

#plot the response matrix(partly)
library(ggplot2)
SaveFile <- 'AllOrResponse.pdf'
plotFlag <-dplot_responseMatrix(resetSFR(response.matrix_non.normalized, "SFR")[2:50,], tag = "Name", base_size = 8)
ggsave(plotFlag, file=SaveFile, width=10, height=8)


#reset the sponteneous firing rate in the response matrix 
ResponseMatrixSFR <-within(resetSFR(response.matrix,"SFR"),rm("Or22c","Or24a","Or67d","pb2A"))
#compute the activation and inhibiton ratio use the globally normalized response matrix
#activation criterior: >= 0.1
#inhibition criterior: <= 0.5*minium response && minimum response <= -0.05
ResponseStatistics <- matrix(0,3,dim(ResponseMatrixSFR)[2])
colnames(ResponseStatistics) <- names(ResponseMatrixSFR)
for(i in 1:dim(ResponseMatrixSFR)[2]){
    TempData <- ResponseMatrixSFR[!is.na(ResponseMatrixSFR[,i]),i]
    #SpontenousRate <- TempData[1]
    ResponseStatistics[1,i] <- length(TempData)
    #browser()
    ResponseStatistics[2,i] <- length(TempData[TempData>=0.15])
    
    if(min(TempData) <= -0.05)
        ResponseStatistics[3,i] <- length(TempData[TempData <= 0.5*min(TempData)])
    else
        ResponseStatistics[3,i] <- 0
}
rownames(ResponseStatistics) <- c('NumOdorants','NumActivation','NumInhibition')

#excitatory probability
ExcitatoryProb <- sum(ResponseStatistics[2,]/ResponseStatistics[1,])/dim(ResponseStatistics)[2]
InhitatoryProb <- sum(ResponseStatistics[3,]/ResponseStatistics[1,])/dim(ResponseStatistics)[2]
#histgram of the probability
hist(ResponseStatistics[2,]/ResponseStatistics[1,],20,xlim = c(0,1),ylim = c(0,25),main = 'Excitatory Response',xlab = 'percetage of Or-OSN response',ylab = 'counts',col="#104E8B",cex=2,cex.lab=2,cex.axis=1.5)
hist(ResponseStatistics[3,]/ResponseStatistics[1,],20,ylim = c(0,25),main = 'Inhibitatory Response',xlab = 'percetage of Or-OSN response',ylab = 'counts',col="#104E8B",cex=2,cex.lab=2,cex.axis=1.5)

#only consider the Ors that have at lease 50 odors
SelectedResponseStatistics <-ResponseStatistics[,ResponseStatistics[1,]>=50]
SelectedExciProb <- sum(SelectedResponseStatistics[2,]/SelectedResponseStatistics[1,])/dim(SelectedResponseStatistics)[2]
SelectedInhiProb <- sum(SelectedResponseStatistics[3,]/SelectedResponseStatistics[1,])/dim(SelectedResponseStatistics)[2]


#from the odorant persepective
#only consider odorants that has been tested on at least 20 OSN
#get rid of SRF, it is all 0
NewResp <- ResponseMatrixSFR[-1,]
RawOdorantRespStat <- apply(NewResp, 1, function(x) length(x[!is.na(x)]))
OdorantRespStat <- NewResp[RawOdorantRespStat>=20,]
NonZero <- apply(OdorantRespStat, 1, function(x) length(x[!is.na(x)]))
PosiOdorRespNum <- apply(OdorantRespStat, 1, function(x) sum(x>=0.2,na.rm=TRUE))
NegaOdorRespNum <- apply(OdorantRespStat, 1, function(x) sum(x<=-0.05,na.rm=TRUE))

#histgram of number of OSN activated by odorant
hist(PosiOdorRespNum/NonZero,20,main = 'Number of Or an odorant activate',xlab = 'Percentage of OSN response',ylab = 'Frequency',col="#104E8B",cex=2,cex.main = 2,cex.lab=2,cex.axis=1.5)
hist(NegaOdorRespNum/NonZero,20,main = 'Number of Or an odorant inhibit',xlab = 'Percentage of OSN response',ylab = 'Frequency',col="#104E8B",cex=2,cex.main = 2,cex.lab=1.5,cex.axis=1.5)

#bar plot of data summary
OdorSummary <- data.frame(
  Odor=rownames(OdorantRespStat),
  Category=factor(rep(c("NoResponse", "Activation","Inhibition"), each= dim(OdorantRespStat)[1])),
  Counts = c(NonZero,PosiOdorRespNum,NegaOdorRespNum)
)
ggplot(OdorSummary, aes(x=Odor, y= Counts,fill=Category))+ geom_bar(stat="identity") + coord_flip()
ggsave('statisticsSelectedOrdoants.eps',width = 6,height = 12)
#overall positive response intensity
hist(log(ResponseMatrixSFR[ResponseMatrixSFR>0]),50,main = 'Histogram of Odorant-OSN response intensity',xlab = 'log(Activity)',col="#104E8B",cex=2,cex.main = 1.5,cex.lab=1.5,cex.axis=1.5)

hist(log(abs(ResponseMatrixSFR[ResponseMatrixSFR<0])),50,main = 'Histogram of Odorant-OSN response intensity',xlab = 'log(Activity)',col="#104E8B",cex=2,cex.main = 1.5,cex.lab=1.5,cex.axis=1.5)
#histgrams of all the Or and Odorant response
require(cowplot)
plots <- list()
for (i in dim(NewResp)[2]){
  #plots[[i]] <- dplot_tuningCurve(receptor = colnames(NewResp)[i], fill.receptor = "#104E8B",base_size = 12)
  plots[[i]] <- dplot_tuningCurve(receptor = colnames(NewResp)[i],response.vector = NewResp[,i], fill.receptor = "#104E8B",base_size = 12)
}
finalPlot <- plot_grid(plotlist = plots[1:20])
ggsave(finalPlot, file="OrResponse.pdf", width=16, height=20)


#histgram of data used
#first transform the statistic into a data frame
DataSummary <- data.frame(
  Or=colnames(ResponseStatistics),
  Category=factor(rep(c("NoResponse", "Activation","Inhibition"), each= dim(ResponseStatistics)[2])),
  Counts = c(ResponseStatistics[1,]-ResponseStatistics[2,]-ResponseStatistics[3,],ResponseStatistics[2,],ResponseStatistics[3,])
)
ggplot(DataSummary, aes(x=Or, y= Counts,fill=Category))+ geom_bar(stat="identity") + coord_flip()
ggsave('statisticsAllOrdoants.eps',width = 6,height = 9)

#histgram of excitatory response
qplot()
