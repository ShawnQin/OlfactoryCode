#this fucntion maps the response at ORN to PNs
PNRatesFun<-function(ORNrates,Rmax,n,sigam,m){
    #ORNrates[ORNrates<0] <- 0 
    sumOfRats<-sum(ORNrates)
    result<-Rmax*ORNrates^n/(ORNrates^n + sigma^n + (m*sumOfRats/190)^n)
    return(result)
}