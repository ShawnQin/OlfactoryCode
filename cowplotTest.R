#plot hist test cowplot
ini<-read.table('ini')
fina<-read.table('finalTOP10')
t1<-ini[1,c(1:8)]
t2<-fina[1,c(1:8)]
ti<-data.frame(table(ini[,1]))
tf<-data.frame(table(fina[,1]))
ti$DynaOrNot<-0
tf$DynaOrNot<-1
temp<-rbind(ti,tf)
temp$Name<-paste('Rank',1,seq="")
for(i in c(2:10)) {
    ti<-data.frame(table(ini[,i]))
    tf<-data.frame(table(fina[,i]))
    ti$DynaOrNot<-0
    tf$DynaOrNot<-1
    t<-rbind(ti,tf)
    t$Name<-paste('Rank',i,seq="")
    temp<-rbind(temp,t)
}