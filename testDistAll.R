testDistAll <- function(resMatr,stMethod,fileName,...)
{
  # this function test the distribution of excitatory or inhibitory 
  # resMatr    the repsonse matrix, all the elements should be larger than 0
  # stMethod    test the distribution, can be exponential, lognormal,or gamma
  # normalizedRate  the spiking rate set as a common mean spiking rate
  # fileName    pdf file name used to store the data
  
  #set the paramters of figure
  saveFile <- paste(fileName,"_",stMethod,".pdf",sep = "")
  pdf(saveFile,width = 8,height = 6)
  par(mfrow=c(2,2),mai=c(0.7,0.8,0.3,0.5))
  
  if(stMethod == "exp"){
    staDist <- "exp"
    if(min(resMatr) < 0){
      stop("data matrix containing negative elements!")
    }
    else{
      dataList <- resMatr[resMatr>0]
    }
    xLable <- "spiking rate (Hz)"
  }
  else if(stMethod == "lnorm"){
    staDist <- "norm"  #since we will transform data into log scale before fitting
    if(min(resMatr) < 0){
      stop("data matrix containing negative elements!")
    }
    else{
      dataList <- log(resMatr[resMatr>0])
    }
    xLable <- "log(spiking rate)"
  }
  else{
    stop("current distribution has to be exp or lnorm!")
  }
  
  # fit the data
  fit <- fitdist(dataList,staDist)

  denscomp(fit,datacol = "gray",fitcol = "red",xlab = xLable,addlegend = FALSE,main = "",lwd = 1.5,cex.lab = 1.5,cex.axis = 1.5)
  cdfcomp(fit,addlegend = FALSE,lwd = 1.5,main = "",cex.lab = 1.5,cex.axis = 1.5)  #compare cdf of emperical data and theoretical
  qqcomp(fit,fitcol= "black",line01col = "red",main = "",addlegend = FALSE,lwd = 1.5,cex.lab = 1.5,cex.axis = 1.5)
  ppcomp(fit,fitcol= "black",line01col = "red",main = "",addlegend = FALSE,lwd = 1.5,cex.lab = 1.5,cex.axis = 1.5)
  dev.off()
  
}