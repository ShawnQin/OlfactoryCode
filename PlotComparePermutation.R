PlotComparePermutation = function(wMatrix,  #the orignal weight matrix
                           pMatrix,  #after permutation
                            ...){
  # this function plot the statistics of interaction matrix befor and after permutation
  # such as the number of positive reponse an OSN has
  
  #interaction matrixs have to have the same dimension
  if (dim(wMatrix)[1]== dim(pMatrix)[1] && dim(wMatrix)[2]== dim(pMatrix)[2]){
    
    ResponseStatistics <- matrix(0,4,dim(wMatrix)[2])
    ResponseStatistics[1,] <- apply(wMatrix,2,function (x) sum(x[x>0]))
    ResponseStatistics[2,] <- abs(apply(wMatrix,2,function (x) sum(x[x<0])))
    ResponseStatistics[3,] <- apply(pMatrix,2,function (x) sum(x[x>0]))
    ResponseStatistics[4,] <- abs(apply(pMatrix,2,function (x) sum(x[x<0]))) 
    
    #plot the statistics of orignal interaction matrix
    Summary1 <- data.frame(
      OSN=c(1:dim(ResponseStatistics)[2]),
      Category=factor(rep(c("Excitation","Inhibition"), each= dim(ResponseStatistics)[2])),
      Counts = c(ResponseStatistics[1,],ResponseStatistics[2,])
    )
    browser()
    p1 <- ggplot(Summary1, aes(x=OSN, y= Counts,fill=Category))+ geom_bar(stat="identity")
    p1 <- p1+ theme_classic()
    ggsave('statisticsOSNresponse.eps',width = 6,height = 2)
    #dev.off()
    
    Summary2 <- data.frame(
      OSN=c(1:dim(ResponseStatistics)[2]),
      Category=factor(rep(c("Excitation","Inhibition"), each= dim(ResponseStatistics)[2])),
      Counts = c(ResponseStatistics[3,],ResponseStatistics[4,])
    )
    ggplot(Summary2, aes(x=OSN, y= Counts,fill=Category))+ geom_bar(stat="identity")+ theme_classic()
#    + theme(axis.line = element_line(colour = "black"))
    
    ggsave('statisticsOSNresponsePermution.eps',width = 6,height = 2)
    #dev.off()
  }
  else{
    stop('the two interaction matrixs have to have the same dimensions!')
  }
}
