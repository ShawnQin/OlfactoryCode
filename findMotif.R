findMotif <- function(newMatx,plotfigure = FALSE)
{
    #=========================================================
    # this function find the defined 7 two nodes interaction motif both in the orignal 
    # and randomly shuffled matrix
    # adjMatix          an interaction matrix
    #plotfigure         whether plot the figure or not
    #==========================================================
    
    #check if the input is matrix
    if(!is.matrix(newMatx)){
        stop('input has to be numeric matrix')
    }
    # adjacent matrix, wiht elements 0, -1 and 1
    adjMatx <- newMatx
    adjMatx[adjMatx > 0] <- 1
    
    #first find the "bi-fan" motif
    countBiFan <- 0  #used to count all the combination of observed bi-fan motif
    negBiFanCount <-  0  # we are extramely interested in cross negative response bi-fan structure
    recordMatx <- matrix(NA,nrow = 0, ncol = 8) #i, j, OSN1, OSN2, resp-i1,respon-i2, respon-j1,respon-j2
    allBifanMatx <- matrix(NA,nrow = 0, ncol = 8) # all combination of bi-fan motif
    allCategoryCount <- rep(0, 7) # 7 different two-node motif
    for (i in 1:(dim(adjMatx)[1]-1)){
        for (j in (i+1):dim(adjMatx)[1]){
            inx1 <- which(adjMatx[i,] != 0) #nonzero elements for odor i
            inx2 <- which(adjMatx[j,] != 0) #nonzero elements for odor j
            sharedOSN <- intersect(inx1, inx2) #shared OSN
            if(length(sharedOSN) >= 2){
                combinationBiFan <- combn(sharedOSN,2)  # each column is a combination
                countBiFan <- countBiFan + dim(combinationBiFan)[2]
                
                for (k in 1:dim(combinationBiFan)[2]){
                    # <- rbind(allBifanMatx,c(i,j,newMatx[i,combinationBiFan[,k]],newMatx[j,combinationBiFan[,k]]))
                    R1 <- adjMatx[i,combinationBiFan[,k]]
                    R2 <- adjMatx[j,combinationBiFan[,k]]
                    if (R1[1]*R2[1] < 0 && R1[2]*R2[2] < 0 && R1[1]*R1[2] < 0){
                        negBiFanCount <- negBiFanCount + 1
                        allCategoryCount[4] <- allCategoryCount[4] + 1
                        recordMatx <- rbind(recordMatx,c(i,j,combinationBiFan[,k],newMatx[i,combinationBiFan[,k]],newMatx[j,combinationBiFan[,k]]))
                    }
                    else if (all(c(R1,R2) > 0)){
                        allCategoryCount[1] <- allCategoryCount[1] + 1
                    }
                    else if (sum(c(R1, R2) > 0) == 3){
                        allCategoryCount[2] <- allCategoryCount[2] + 1
                    }
                    else if (sum(c(R1, R2) > 0) == 2 && (all(R1>0) | all(R1 < 0))){
                        allCategoryCount[3] <- allCategoryCount[3] + 1
                    }
                    else if (sum(c(R1, R2) > 0) == 2 && R1[1]*R2[1]>0){
                        allCategoryCount[5] <- allCategoryCount[5] + 1
                    }
                    else if (sum(c(R1, R2) > 0) == 1){
                        allCategoryCount[6] <- allCategoryCount[6] + 1
                    }
                    else if (sum(c(R1, R2) > 0) == 0){
                        allCategoryCount[7] <- allCategoryCount[7] + 1
                    }
                }
            }
        }
    }
    
    
    numShuffle <- 100  #try 100 times
    allNumberExist <- matrix(0, nrow = numShuffle, ncol = 2)
    #shuffling matrix
    allCategoryCountShuffled <- matrix(0, nrow = numShuffle, ncol = 7) # 7 different two-node motif
    for (l0 in 1:numShuffle){
        # shuffledMatx <- ShuffledData(newMatx,WhichWay = "both")
        shuffledMatx <- RandomizedMatrx(newMatx,method="switch")
        adjMatxShuffled <- shuffledMatx
        adjMatxShuffled[shuffledMatx > 0] <-1
        countBiFanShuffled <- 0
        negBiFanCountShuffled <- 0
        for (i in 1:(dim(adjMatxShuffled)[1]-1)){
            for (j in (i+1):dim(adjMatxShuffled)[1]){
                inx1 <- which(adjMatxShuffled[i,] != 0) #nonzero elements for odor i
                inx2 <- which(adjMatxShuffled[j,] != 0) #nonzero elements for odor j
                sharedOSN <- intersect(inx1, inx2) #shared OSN
                if(length(sharedOSN) >= 2){
                    combinationBiFan <- combn(sharedOSN,2)  #each column is a combination
                    countBiFanShuffled <- countBiFanShuffled + dim(combinationBiFan)[2]
                    for (k in 1:dim(combinationBiFan)[2]){
                        R1 <- adjMatxShuffled[i,combinationBiFan[,k]]
                        R2 <- adjMatxShuffled[j,combinationBiFan[,k]]
                        if (R1[1]*R2[1] < 0 && R1[2]*R2[2] < 0 && R1[1]*R1[2] < 0){
                            negBiFanCountShuffled <- negBiFanCountShuffled + 1
                            allCategoryCountShuffled[l0,4] <- allCategoryCountShuffled[l0,4] + 1
                            # recordMatx <- rbind(recordMatx,c(i,j,combinationBiFan[,k]))
                        }
                        else if (all(c(R1,R2) > 0)){
                            allCategoryCountShuffled[l0,1] <- allCategoryCountShuffled[l0,1] + 1
                        }
                        else if (sum(c(R1, R2) > 0) == 3){
                            allCategoryCountShuffled[l0,2] <- allCategoryCountShuffled[l0,2] + 1
                        }
                        else if (sum(c(R1, R2) > 0) == 2 && (all(R1>0) | all(R1 < 0))){
                            allCategoryCountShuffled[l0,3] <- allCategoryCountShuffled[l0,3] + 1
                        }
                        else if (sum(c(R1, R2) > 0) == 2 && R1[1]*R2[1]>0){
                            allCategoryCountShuffled[l0,5] <- allCategoryCountShuffled[l0,5] + 1
                        }
                        else if (sum(c(R1, R2) > 0) == 1){
                            allCategoryCountShuffled[l0,6] <- allCategoryCountShuffled[l0,6] + 1
                        }
                        else if (sum(c(R1, R2) > 0) == 0){
                            allCategoryCountShuffled[l0,7] <- allCategoryCountShuffled[l0,7] + 1
                        }
                    }
                }
                
            }
        }
        allNumberExist[l0,] <- c(countBiFanShuffled,negBiFanCountShuffled)
    }
    
    # Zscore of different moitif
    Zscore <- abs((colMeans(allCategoryCountShuffled)-allCategoryCount))/apply(allCategoryCountShuffled, 2, sd)
    # statistics of different motif in ranomized network
    cateRandStat <- data.frame(
        sys = c(rep('random',7),rep('original',7)),
        motif = rep(c('type1', 'type2','type3','type4','type5','type6','type7'),2),
        meanNum = c(colMeans(allCategoryCountShuffled),allCategoryCount),
        std = c(apply(allCategoryCountShuffled, 2, sd),rep(NA,7))
    )
    
    if(plotfigure){
    #bar plot of total bi-fan motif and negative motif
        mycolor <- brewer.pal(11,"Spectral")
        motifStatistics <- data.frame(
        category = c("total", "cross negative"),
        meanNumber = colMeans(allNumberExist),
        std = c(sd(allNumberExist[,1]),sd(allNumberExist[,2]))
        )
        gmtif <- ggplot(motifStatistics,aes(x=category,y=meanNumber,color=category)) + geom_bar(stat="identity", color="black", position=position_dodge()) + geom_errorbar(aes(ymin=meanNumber, ymax=meanNumber + std), width=.2,position=position_dodge(.9)) 
    
        ## bar plot of different motif in randomized network
        gmotifrand <- ggplot(cateRandStat,aes(x = motif,y = meanNum,fill = sys,width = 0.9)) + geom_bar(stat = "identity",color = "black", position = position_dodge()) + geom_errorbar(aes(ymin=meanNum,ymax = meanNum + std), width = 0.2, position = position_dodge(.9))
        gmotifrand <- gmotifrand + theme_classic() + scale_fill_manual(values = mycolor[c(1,10)]) + theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5)) + labs(x = "motif",y="counts")+ theme(axis.text = element_text(size=14), axis.title = element_text(size=18,face="bold"))
        }
    
    #return all the data in a list
    motifStsts <- list(bifanRecord = recordMatx,bifanMatri = allCategoryCount, allMotif = allNumberExist, allMotifShuff = allCategoryCountShuffled)
    return(motifStsts)
}