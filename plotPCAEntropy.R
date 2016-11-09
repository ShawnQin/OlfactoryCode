#this function plot the result of PCA and total entropy comparing
plotPCAEntropy<-function(compNumMix, #vetor of number of odors
               sumEntropy, #dada frame of sum of entropy
               pcacmp,     #data frame of pca results
               balanceMethod, #balance method used
               odorCon,    #odor concentration
               Rtype,   #interaction type, competitive or multiple binding
               
               ...)
{
  #set the color
  myColor <- brewer.pal(11,"Spectral")
  
  if(!missing(sumEntropy)){
      #plot the sum of entropy
    yRange <-c(40,120)
    par(bty="l")
    plot(compNumMix,sumEntropy[,1],type = "b",pch=16,lty=1,lwd=2,cex=2,col=myColor[1],ylim = yRange,xlab = "Number of odor",ylab = "Total Entropy (bit)")
    lines(compNumMix,sumEntropy[,2],type = "b",pch=16,lty=1,lwd=2,cex=2,col=myColor[10])
    # dev.off()
  }
  #use ggplot2
  # sumEntro<data.frame(sys=factor(rep(c("Ori","No"),each=dim())))
  # peig <- ggplot(subset(pcacmp,pcacmp$system=="Ori"),aes(x=component, y=eigenvalus, group=numMix)) + geom_line(data=subset(pcacmp,pcacmp$system=="Ori"),linetype="solid",color=myColor[1])+geom_point(data=subset(pcacmp,pcacmp$system=="Ori"),aes(shape=numMix),color=myColor[1],size=4)+theme(legend.position="top")
  
  
  if(!missing(pcacmp)){
  #plot pca result, variance explained and eigen value
    pvar <- ggplot(subset(pcacmp,pcacmp$system=="Ori"),aes(x=component, y=variance_expl, group=numMix)) + geom_line(data=subset(pcacmp,pcacmp$system=="Ori"),linetype="solid",color=myColor[1])+geom_point(data=subset(pcacmp,pcacmp$system=="Ori"),aes(shape=numMix),color=myColor[1],size=4)+theme(legend.position="top") 
  
    pvar+ geom_line(data=subset(pcacmp,pcacmp$system=="Noi"),aes(x=component, y=variance_expl, group=numMix),linetype="dashed",color=myColor[10])+geom_point(data=subset(pcacmp,pcacmp$system=="Noi"),aes(shape=numMix),color=myColor[10],size=4)#+geom_hline(aes(yintercept = 1),linetype=3)
  
    # browser()
    varFile <-paste("comp_var_method",as.character(balanceMethod),"con1_",as.character(odorCon),"_",Rtype,".pdf",sep = "")
    ggsave(varFile,width = 8,height = 6)
  

    peig <- ggplot(subset(pcacmp,pcacmp$system=="Ori"),aes(x=component, y=eigenvalue, group=numMix)) + geom_line(data=subset(pcacmp,pcacmp$system=="Ori"),linetype="solid",color=myColor[1])+geom_point(data=subset(pcacmp,pcacmp$system=="Ori"),aes(shape=numMix),color=myColor[1],size=4)+theme(legend.position="top")
  
    peig+ geom_line(data=subset(pcacmp,pcacmp$system=="Noi"),aes(x=component, y=eigenvalue, group=numMix),linetype="dashed",color=myColor[10])+geom_point(data=subset(pcacmp,pcacmp$system=="Noi"),aes(shape=numMix),color=myColor[10],size=4)#+geom_hline(aes(yintercept = 1),linetype=3)
  
    eigFile <-paste("comp_eig_method",as.character(balanceMethod),"con1_",as.character(odorCon),"_",Rtype,".pdf",sep = "")
  ggsave(eigFile,width = 8,height = 6)
  
  
  # #effective dimension
  effectivDim <- aggregate(pcacmp$eigenvalue,by=list(system = pcacmp$system,numMix=pcacmp$numMix),function(x) sum(x>=1))
  #   effectivDim <-data.frame(
  #     system=factor(c("Ori","Ori","Ori", "Noi","Noi","Noi")),
  #     numMix=c(c(5,20,35),c(5,20,35)),
  #     dimen=c(c(10,7,3),c(8,5,3))
  # )
  # 
    gEffDim<-ggplot(effectivDim,aes(x=numMix,y=x,group=system,color=system)) + geom_point(aes(x=numMix,y=x,shape=numMix),size=4) + geom_line() + scale_y_continuous(limits = c(0,10)) + scale_color_manual(values=c(myColor[10],myColor[1]))
    gEffDim <- gEffDim + labs(x="number of odor", y = "effective coding dimension")
    effDimFile <-paste("comp_effDim_method",as.character(balanceMethod),"con1_",as.character(odorCon),"_",Rtype,".pdf",sep = "")
    ggsave(effDimFile,width = 4,height = 3)
  }
}