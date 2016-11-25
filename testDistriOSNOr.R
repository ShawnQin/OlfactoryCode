testDistriOSNOr <- function(resMatr,stMethod,normalizedRate = 50,fileName,...)
{
  # this function test the distribution of odor-OSN response
  # resMatr    the repsonse matrix, all the elements should be larger than 0
  # stMethod    test the distribution, can be exponential, lognormal,or gamma
  # normalizedRate  the spiking rate set as a common mean spiking rate
  # fileName    figure name
  
  
  COL <- dim(resMatr)[2]  #number of variables
  ROW <- dim(resMatr)[1]  #number of measurement
  
  # figure name and plot setting
  saveFile = paste(fileName,"_",stMethod,".pdf",sep ="")
  pdf(saveFile,width = 8,height = 6)
  par(mfrow=c(2,2),mai=c(0.7,0.8,0.3,0.5))
  
  # response of each variables
  if(stMethod == "exp"){
    # each OSN or odor
    plot.ecdf(resMatr[,1],do.points=FALSE,xlim = c(0,300),main = "",xlab = "spiking rate (HZ)", ylab = "cdf",cex.lab = 1.5,cex.axis = 1.5)
    for (i0 in 2:dim(resMatr)[2]){
      temp <- ecdf(resMatr[,i0])
      plot(temp,verticals=TRUE, do.points=FALSE,add=TRUE)
    }
    
    # normalize
    normVal <- normalizedRate # normalized average sipking rate
    averageSpikes <- apply(resMatr, 2, mean)
    browser()
    afterScaleMatr <- resMatr*matrix(rep(normVal/averageSpikes,ROW),nrow = ROW, ncol = COL, byrow = TRUE)  #normalized to average spking rate as 100 Hz
    index <- seq(1,300,1)
    allCdf <- matrix(0,nrow = length(index),ncol = dim(afterScaleMatr)[2])
    
    fit1 <- ecdf(afterScaleMatr[,1])
    allCdf[,1] <- fit1(index)
    plot(fit1,do.points=FALSE,xlim = c(0,300),main ="",ylab = "cdf",xlab = "spiking rate (Hz)",cex.lab = 1.5,cex.axis = 1.5)
    for (i0 in 2:dim(afterScaleMatr)[2]){
      temp <- ecdf(afterScaleMatr[,i0])
      plot(temp,verticals=TRUE, do.points=FALSE,add=TRUE)
      allCdf[,i0] <- temp(index)
    }
    
    #generate same number of exponential distribution
    plot(ecdf(rexp(ROW,rate = 1/normVal)),verticals=TRUE, do.points=FALSE,main = "",xlim = c(0,300),xlab = "spiking rate (Hz)",ylab = "cdf",cex.lab = 1.5,cex.axis = 1.5)
    for (i0 in 2:dim(afterScaleMatr)[2]){
      plot(ecdf(rexp(ROW,rate = 1/normVal)),verticals=TRUE, do.points=FALSE,add=TRUE)
    }
    
    #compare the average cdf with theoretic results
    plot(index,apply(allCdf, 1,mean),main = "",xlab = "spiking rate (Hz)",ylab = "cdf",cex.lab = 1.5,cex.axis = 1.5)
    curve(pexp(x, rate = 1/normVal), col = "red", lwd =3,add = TRUE) #overlay an exponential distribution
    dev.off()
  }
  
  else if (stMethod == "lnorm"){
    #overall fitting paramter
    allFit <- fitdist(log(resMatr[resMatr>0]),"norm")
    meanStd <- allFit$estimate
    #browser()
    # plot original the cdf in log scale, only consider values lareger than 0
    plot.ecdf(log(resMatr[resMatr[,1]>0,1]),do.points=FALSE,xlim = c(0,6),main = "",xlab = "log(spiking rate)", ylab = "cdf",cex.lab = 1.5,cex.axis = 1.5)
    for (i0 in 2:dim(resMatr)[2]){
      temp <- ecdf(log(resMatr[resMatr[,i0]>0,i0]))
      plot(temp,verticals=TRUE, do.points=FALSE,add=TRUE)
    }
    
    # normalized individule cdf
    normVal <- normalizedRate # normalized average sipking rate
    averageSpikes <- apply(resMatr, 2, mean)
    afterScaleMatr <- resetPosiMatr*matrix(rep(normVal/averageSpikes,ROW),nrow = ROW, ncol = COL, byrow = TRUE)  #normalized to average spking rate as 100 Hz
    index <- seq(0,6,0.05)
    allCdf <- matrix(0,nrow = length(index),ncol = dim(afterScaleMatr)[2])
    
    fit1 <- ecdf(log(afterScaleMatr[afterScaleMatr[,1] >0,1]))
    allCdf[,1] <- fit1(index)
    plot(fit1,do.points=FALSE,xlim = c(0,6),main ="",ylab = "cdf",xlab = "log(spiking rate)",cex.lab = 1.5,cex.axis = 1.5)
    for (i0 in 2:dim(afterScaleMatr)[2]){
      temp <- ecdf(log(afterScaleMatr[afterScaleMatr[,i0]>0,i0]))
      plot(temp,verticals=TRUE, do.points=FALSE,add=TRUE)
      allCdf[,i0] <- temp(index)
    }
    
    
    
    #generate lognormal distribution
    plot(ecdf(rnorm(ROW,meanStd[1],meanStd[2])),verticals=TRUE, do.points=FALSE,main = "",xlim = c(0,6),xlab = "log(spiking rate)",ylab = "cdf",cex.lab = 1.5,cex.axis = 1.5)
    for (i0 in 2:dim(afterScaleMatr)[2]){
      plot(ecdf(rnorm(ROW,meanStd[1],meanStd[2])),verticals=TRUE, do.points=FALSE,add=TRUE)
    }
    
    #average cdf
    plot(index,apply(allCdf, 1,mean),main = "",xlab = "log(spiking rate)",ylab = "cdf",cex.lab = 1.5,cex.axis = 1.5)
    curve(pnorm(x, meanStd[1],meanStd[2]), col = "red", lwd =3,add = TRUE) #overlay an exponential distribution
    dev.off()
    
  }
  else{
    stop("currently, only exponential or lognormal distribution is supported!")
  }


}