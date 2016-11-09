#this small pice of script extract the hist count of the odor Response of each OSN
#last revised on June 16th,2016
histCount<-function(
  numBin=10,
  ORNmax=250,
  baseline=30,
  X,
  Y,
  ...
){
numBin<-10
ORNmax <-250
baseline <-30
# cutX<- (0:numBin)*ORNmax/numBin  #with inhibitonn
cutY<- baseline-0.1 + (0:numBin)*(ORNmax-baseline+0.1)/numBin #no inhibitio
# browser()
countX <-matrix(data=NA,nrow=42,ncol=numBin)
countY <-matrix(data=NA,nrow=42,ncol=numBin)
for(i in 1:dim(X)[2]){
  cutX<-(0:numBin)*max(cbind(X[,i],Y[,i]))/numBin
  temp1 <- hist(X[,i],cutX,plot=FALSE)
  countX[2*i-1,]<-temp1$breaks[-1]
  countX[2*i,]<-temp1$counts
  temp2<- hist(Y[,i],cutX,plot=FALSE)
  countY[2*i-1,]<-temp2$breaks[-1]
  countY[2*i,]<-temp2$counts
}
write.xlsx(countX,'histSpikEachOSN-multi-N20-Sample.xlsx')
write.xlsx(countY,'histSpikEachOSN-noInhi-multi-N20-Sample.xlsx')
}