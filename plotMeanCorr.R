#plot mean response and correlation matrix
plotMeanCorr <- function(X,#response matrix for the original
                         Y,#response matrix for the original
                         ORNname, #names of each ORN
                         filename, #name to save the plot
                         plotOrNot=TRUE, #whether plot or not
                         baseline=30,  #the default basal activity of OSN
                         ORNmax,  #maximum activity
                         numBin=250,
                         ...
)
{
  
    #if plot the mean response of each OSN
  if(plotOrNot){
    library(ggplot2)
    newData <- data.frame(
        OrnNames = rep(factor(ORNname,as.character(ORNname)),times=2),
        system = factor(rep(c("Ori", "NoI"), each= dim(X)[2])),
        meanResp = c(apply(X,2,mean),apply(Y,2,mean))                 
    )
    
    # browser()
    gplot <- ggplot(newData,aes(x = OrnNames,y =meanResp,fill= system))+geom_bar(stat="identity", position=position_dodge())+labs(x='ORN',y="mean response")+theme(axis.text.x=element_text(size = 14,angle = 60,hjust = 1)) #geom_abline(intercept = median(range(X)),slope = 0) 
    saveFile <- paste(filename,"_hist.pdf",sep = "")
    ggsave(gplot, file=saveFile, width=10, height=8)
  }
    
    #bar plot of the entropy of each ORN
    library(entropy)  
    rangeX <- apply(X,2,range)
    rangeY <- apply(Y,2,range)
    entropyX <- rep(0,times = dim(X)[2])
    entropyY <- rep(0,times = dim(X)[2])
    cutX<- (0:numBin)*ORNmax/numBin  #with inhibiton
    cutY<- baseline-0.1 + (0:numBin)*(ORNmax-baseline+0.1)/numBin #no inhibition
    
    for(i in 1:dim(X)[2]){
        # cutX<- floor(rangeX[1,i]):ceiling(rangeX[2,i])
        # cutY <- floor(rangeY[1,i]):ceiling(rangeY[2,i])
        countX <- hist(X[,i],cutX,plot=FALSE)
        countY <- hist(Y[,i],cutY,plot=FALSE)
        entropyX[i] <- entropy.empirical(countX$counts,unit = "log2")
        entropyY[i] <- entropy.empirical(countY$counts,unit = "log2")
        #entropyX[i] <- -sum(countX$density*log2(countX$density))
        #entropyY[i] <- -sum(countY$density*log2(countY$density))
    }
    
    allEntr <- data.frame(
      entropyOri=entropyX,
      entropyNo=entropyY
    )

    
    #plot the entropy of each OSN
    if(plotOrNot==TRUE){
        allEntropy <- data.frame(
        OrnNames = rep(factor(ORNname,as.character(ORNname)),times=2),
        system = factor(rep(c("Ori", "NoI"), each= dim(X)[2])),
        entropyORN = c(entropyX,entropyY)  
      )
      gplot <- ggplot(allEntropy,aes(x = OrnNames,y =entropyORN,fill= system))+geom_bar(stat="identity", position=position_dodge())+labs(x='ORN',y="entropy(bit)")+theme(axis.text.x=element_text(size = 14,angle = 60,hjust = 1))+ylim(0, 6)
      entropyFile <- paste(filename,"_entropy.pdf",sep = "")
      ggsave(gplot, file=entropyFile, width=10, height=5)
#     browser()
    #plot the correlation matrix and correlation coefficient distribution
    #comparing of correlation coefficients
    Mx <- cor(X)
    My <- cor(Y)
    file_corr <- paste("corr_comp_",filename,".pdf",sep = "")
    plotRng<- range(c(Mx[lower.tri(Mx)],My[lower.tri(My)]))
    pdf(file_corr,width = 6,height = 6)
    plot(Mx[lower.tri(Mx)],My[lower.tri(My)],xlim = plotRng,ylim = plotRng,xlab="correlation coefficient between ORNs(O)",ylab="correlation coefficient between ORNs(N)",cex.axis=1.5,cex.lab = 1.5)
    abline(a=0,b=1,lwd=2,col="gray")
    dev.off()
    
    library(easyGgplot2) 
    file_corr_hist <-paste("corr_hist_",filename,".pdf",sep = "")
    newData <- data.frame(
        system=factor(rep(c("O", "N"), each= length(Mx[lower.tri(Mx)]))),
        corrVal= c(Mx[lower.tri(Mx)],My[lower.tri(Mx)])            
        )
    ght <- ggplot2.histogram(data=newData, xName='corrVal',
                      groupName='system',bins = 20, legendPosition="top",
                      alpha=0.5, position="identity",backgroundColor="white")
                      #removePanelGrid=TRUE,removePanelBorder=TRUE,
                      #axisLine=c(0.5, "solid", "black"))
    ght+theme(panel.grid =element_blank(),axis.line.y = element_line(color="black", size = 0.5),axis.line.x = element_line(color="black", size = 0.5),axis.text=element_text(size=18),axis.text.x = element_text(size=14))
    #ght <- ggplot(newData, aes(x=corrVal, color=system)) +
        #geom_histogram(fill="white", position="dodge")+scale_color_grey()+theme_classic()
    ggsave(ght, file=file_corr_hist, width=8, height=6)
    
    #source("rquery.cormat.R")
    fn1 <- paste(filename,"_corrMatrx_ori.pdf")
    fn2 <- paste(filename,"_corrMatrx_NoInhi.pdf")
    
    #col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
    library(corrplot)
    pdf(fn1,width = 10,height = 10)
    corrplot(Mx, method="color",#col=col(200),
             type="upper", order="hclust", 
             #addCoef.col = "black", # Add coefficient of correlation
             tl.col="black", tl.srt=45,tl.cex=1.5, #Text label color and rotation
             # Combine with significance
             #p.mat = M$p, sig.level = 0.001, insig = "blank", 
             # hide correlation coefficient on the principal diagonal
             diag=FALSE #don't plot diagnal region
    )
    dev.off()
    
    pdf(fn2,width = 10,height = 10)
    corrplot(My, method="color",#col=col(200),
             type="upper", order="hclust", 
             #addCoef.col = "black", # Add coefficient of correlation
             tl.col="black", tl.srt=45,tl.cex = 1.5, #Text label color and rotation
             # Combine with significance
             #p.mat = M$p, sig.level = 0.001, insig = "blank", 
             # hide correlation coefficient on the principal diagonal
             diag=FALSE #don't plot diagnal region
    )
    dev.off()
    }
    
    return(allEntr)
}