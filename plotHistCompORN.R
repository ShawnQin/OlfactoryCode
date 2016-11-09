#this function plot all histogram respone of the ORNs
plotHistCompORN <- function(X,
                            Y,
                            fileName="histoCompare.pdf",
                            ...){
    library(ggplot2)
    library(cowplot)
    library(easyGgplot2)
    #pdf(fileName,width = 20,height = 20) #set the size of the graphics
    #op<-par(mfrow = c(5,5))
    #dev.off()
    plots <- list()
    for(i in 1:dim(X)[2])
    {       
        xRange <- range(c(range(X[,i]),range(Y[,i])))
        newData <- data.frame(
            system = factor(rep(c("O", "I"), each= dim(X)[1])),
            spiking = c(X[,i],Y[,i])                 
            )
        #plots[[i]]<-ggplot(newData,aes(x=spiking,fill=system))+geom_histogram(bins=10,position='dodge',stat="identity")+labs(x='spiking rate',y="Frequency")+ggtitle(colnames(X)[i])
        plots[[i]] <-ggplot2.histogram(data=newData, xName='spiking',
                          groupName='system',bins = 10, legendPosition="top",
                          alpha=0.5, position="identity")
        #plots[[i]]
        #subp<-ggplot(newData, aes(x=spiking, fill=system, color=system)) +
        #    geom_histogram(position="identity",alpha=0.5)
        #print(subp,vp = viewport(0.5, 0.75, 4, 4))
        #hist(X[,i],xlim=xRange,col='skyblue',border=F,main = NULL,xlab = "spiking rate(spikes/s)",cex.lab=1.5,cex.axis=1.5)
        #hist(Y[,i],add=T,col=scales::alpha('red',.5),border=F,cex.lab=1.5,cex.axis=1.5)
        #par(pin=c(3,3),mai=c(0.6,0.6,0.2,0.6)) #set the figure size
    }
    #dev.new()
    #pdf(fileName,width = 20,height = 20)
    finalPlot <- plot_grid(plotlist = plots)
    ggsave(finalPlot, file=fileName, width=20, height=20)
    #dev.off()  #close the graphic device  
}
