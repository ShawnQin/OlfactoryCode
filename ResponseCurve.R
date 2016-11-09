# compare response curve
# this program plot the response curve of the multiple binding and compeitive binding model
# and also plot the curve shift due to the existene of inhibitory odors

#load liberies
library(xlsx)
library(RColorBrewer) #this package is used for color settings
library(ggplot2)

#source the functions
source("setBindAffiInhiMatx.R")
source("digitMatx.R")

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
DigitSp <- c(15,30,50,100,150,200)
alpha <- ORNmax/spontaneous - 1
OdorNum <- dim(newMatx)[1] #number of odorants

#specify the interaction time and balance method
balance_method = 1   #OSN-wise
ResType = "comp"
Beta <- 1
#set the binding affinity matrix
equilMatx <- setBindAffiInhiMatx(newMatx,DigitSp,ORNmax,alpha,Beta,balance_method,ResType)

#select one odor-OSN pair and plot their response curves
OdorantExci <- 17
OdorantInhi <- 16
SelectedOSN <- 2
KA <- equilMatx[OdorantExci,SelectedOSN]
KI <- equilMatx[OdorantInhi,SelectedOSN]
Concentration <- 10^(seq(-2,3,by=0.1))
InhiCon <- c(0.1,1,10)/KI    #backgroudn inhibitory response

#for multiple binding 
ResponseInhiBck <- matrix(0, nrow = length(Concentration),ncol = length(InhiCon))
    ResponseExci <- ORNmax/(1 + alpha/(1 + Concentration*KA))
    for (i in 1:length(InhiCon))
    ResponseInhiBck[,i] <- ORNmax/(1 + alpha*(1+InhiCon[i]*KI)/(1 + Concentration*KA))

#plot the excitatory odor response curve and with the background of inhibitory odor
myColor <- brewer.pal(9,"Blues")
dataForPlot <- data.frame(
    sys = factor(c(rep("Exci",times=length(ResponseExci)),rep("InhiBck",times=dim(ResponseInhiBck)[1]*dim(ResponseInhiBck)[2]))),
    con = rep(Concentration,4),
    ref = factor(c(rep(1,times = length(ResponseExci)),ref = sort(rep(InhiCon,length(Concentration))))),
    res = c(ResponseExci,as.vector(ResponseInhiBck))
)
p1 <- ggplot(subset(dataForPlot,dataForPlot$sys=="InhiBck"),aes(x = con, y = res, group = ref,color = ref)) + geom_line()
p1 <- p1 + geom_line(data = subset(dataForPlot,dataForPlot$sys=="Exci"),aes(x = con, y = res,color="#8B0000")) + scale_color_manual(values = c(myColor[c(3,6,9)],"#8B0000")) + scale_x_log10()
p1 <- p1 + labs(x="odor concentration",y="response")
p1 <- p1 + theme_classic()
    # + geom_line(data=subset(dataForPlot,dataForPlot$sys=="Exci"),aes(x = con, y = res,group = con))
