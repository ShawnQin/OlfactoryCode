#sampling mixture from the 107 odorants with same concentration
#and return the response matrix
mixResp<- function(MixNum=5,#number of odorants mixture
                   SampNum=10000,#number of mixtures sampled
                   ORNmax,
                   alpha,
                   weightMatx,
                   equilMatx,
                   conc=1, #concentration of each odorant
                   ResType="multi",#respononse to mixture, mutliple binding sites, or competitive binding
                   ...
)
{   
    OdorNum <- dim(weightMatx)[1]   #total number of odorants
    allResp <- matrix(0,SampNum,dim(weightMatx)[2])  #initialize the response matrix of mixture
    #browser()
    allInx <- matrix(0,SampNum,MixNum) #initialize the sampled odorants index
    #set.seed(1)  #to test whether two repeat calculation the same or not?
    for(i in 1:SampNum){
        allInx[i,]<- sample(1:OdorNum,MixNum)
        allResp[i,]<-mixRespFunc(allInx[i,],ORNmax,alpha,weightMatx,equilMatx,conc,ResType)
    }
    # browser()
    colnames(allResp) <- colnames(weightMatx) #set the name of each column
    return(allResp)
}