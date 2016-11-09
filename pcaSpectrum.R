#priciple component analysis and comparing
#this function take two response matrix as input and calcuate their PCA spectrum
pcaSpectrum <- function(X, #original response matrix
                        Y,  #response matrix without inhibition
                        numMix, #number of mixture
                        scale=FALSE, #scale the data or not
                        fileName,  #file name to store PCA spectrum
                        plot=TRUE,
                        ...)
{
    #comparing the dimension of the data
    if(all(dim(X)==dim(Y))){
        if(scale==TRUE){
            stdThreshold <- 1e-4 #thresold of standard deviation
            allSdX<-apply(X,2,sd)
            allSdY <- apply(Y,2,sd)
            Xnew <- X[,allSdX>stdThreshold] #only retain the variable whose variance are larger than the threshold
            Ynew <- Y[,allSdY>stdThreshold]
            #for orignial data set
            pcaMix <- prcomp(Xnew,scale=T)
            #for the system without inhibition
            pcaMixNoInhi <- prcomp(Ynew, scale=T)
        }
        else{
            #for orignial data set
            pcaMix <- prcomp(X,scale=FALSE)
            #for the system without inhibition
            pcaMixNoInhi <- prcomp(Y,scale=FALSE)
        }
        
        # Eigenvalues
        eig <- (pcaMix$sdev)^2
        # Variances explained by each component in percentage
        variance.explained <- eig*100/sum(eig)
                
        # Eigenvalues
        eigNoInhi <- (pcaMixNoInhi$sdev)^2
        # Variances explained by each component in percentage
        variance.explained.NoInhi<- eigNoInhi*100/sum(eigNoInhi)
        
        #data frame for return
        pcaComp <- list(eig = eig,variance = variance.explained,eigNoInhi = eigNoInhi,varianceNoInhi=variance.explained.NoInhi)
        
        if(plot){
        #plot contributions of each components
        yRange <- range(c(variance.explained,variance.explained.NoInhi))
        pdf(paste("compare PCA spectrum",fileName,".pdf",sep=""),width = 6,height = 6) 
        par(lwd=2,cex=1.5,font.lab = 2)
        plot(1:length(variance.explained),variance.explained,type = "b",pch=20,lty = 1,col="blue",ylim = yRange,main=paste(as.character(numMix),"random odorants mixture"),xlab = "component",ylab = "variance explained (%)")
        lines(1:length(variance.explained.NoInhi),variance.explained.NoInhi,type = "b",pch=18,lty = 1,col="red")       
        legend("topright",legend = c("original","without inhibition"),lty=c(1,1),pch = c(20,18),col = c("blue","red"),box.lty=0)
        dev.off()
        
        #plot eigen values
        yRange <- range(c(eig,eigNoInhi))
        pdf(paste("PCA eigenvalues",fileName,".pdf",sep=""),width = 6,height = 6) 
        par(lwd=2,cex=1.5,font.lab = 2)
        plot(1:length(variance.explained),eig,type = "b",pch=20,lty = 1,col="blue",ylim = yRange,main=paste(as.character(numMix),"random odorants mixture"),xlab = "component",ylab = "eigenvalue")
        lines(1:length(variance.explained.NoInhi),eigNoInhi,type = "b",pch=18,lty = 1,col="red")
        legend("topright",legend = c("original","without inhibition"),lty=c(1,1),pch = c(20,18),col = c("blue","red"),box.lty=0)
        abline(a=1,b=0,lwd=2,col="gray")
        dev.off()
        }
    }
    else{
        stop("The dimension of X and Y are different, please check the data!")
    } #end if
    return(pcaComp)
}