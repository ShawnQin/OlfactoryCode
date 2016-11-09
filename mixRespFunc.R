#define the respone fucntion of mixture, provided by Yuhai Tu
mixRespFunc <- function(odorMix, #mix of odorant, can be indexed through 1 to 110
                        ORNmax,
                        alpha,
                        weightMatx,
                        equilMatx,
                        conc = 1,  #concentration of odorant
                        ResType="multi", #respononse to mixture, mutliple binding sites, or competitive binding
                        ...){
    
    K<- rbind(equilMatx[odorMix,])
    W <- rbind(weightMatx[odorMix,])
    if(ResType=="multi"){
        interac <- rep(1,times = dim(weightMatx)[2])  #initialize the response
        for(i in 1:length(odorMix)){
            interac <- interac*(1 + conc*K[i,])^(-W[i,])
        }
        response<- ORNmax/(1+alpha*interac)
        return(response)
    }
    else if(ResType=="comp"){
        reac_ini<- rep(0,times = dim(weightMatx)[2])  #initialize the response
        reac_exc<- rep(0,times = dim(weightMatx)[2])  #initialize the response
        for(i in 1:dim(weightMatx)[2]){
            for(j in 1:length(odorMix)){
                if(W[j,i] > 0){
                    reac_exc[i] = reac_exc[i] + conc*K[j,i]
                }
                else if(W[j,i] < 0){
                    reac_ini[i] = reac_ini[i] + conc*K[j,i]
                }
            }
        }
        response <- ORNmax/(1+alpha*(1+reac_ini)/(1+reac_exc))
        return(response)
    }
    else{
        stop("method has to be multi or comp!")
    }
}