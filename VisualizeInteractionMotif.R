# this program visulize the odor-OSN interaction
# and try to find motif
# this is a bipartite network
# we use the interaction matrix from Hallem 2006 paper
# last revised on 10/27/2016
##===================================================##

#load the packages
#load liberies
library(xlsx)
library(RColorBrewer)  #this package is used for color settings
library(ggplot2)
library(igraph)
library(network)
library(sna)
library(visNetwork)

#source some functions
source("ShuffledData.R")
source("RandomizedMatrx.R")

#load data file
digiFile<-"/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/CarslonORNdigit.xlsx"
rawData1<-read.xlsx(digiFile,1)
rawMatx <- as.matrix(rawData1[,3:length(rawData1)])
chemType <- as.character(rawData1[,2])
rownames(rawMatx)<- rawData1$odorantName[1:110]

#delete odorants that don't elicite any response
adjustSpikingMatx<-rawMatx     #estimated absolute ORNs spiking rate
temp<-apply(abs(adjustSpikingMatx),1,sum)
newMatx <- adjustSpikingMatx[which(temp >0),]
allType <- chemType[which(temp >0)]
##====================== finding motif ===============
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

#================== randomly shuffling ================
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

#bar plot of total bi-fan motif and negative motif
mycolor <- brewer.pal(11,"Spectral")
motifStatistics <- data.frame(
    category = c("total", "cross negative"),
    meanNumber = colMeans(allNumberExist),
    std = c(sd(allNumberExist[,1]),sd(allNumberExist[,2]))
)
gmtif <- ggplot(motifStatistics,aes(x=category,y=meanNumber,color=category)) + geom_bar(stat="identity", color="black", position=position_dodge()) + geom_errorbar(aes(ymin=meanNumber, ymax=meanNumber + std), width=.2,
                  position=position_dodge(.9)) 

## bar plot of different motif in randomized network
gmotifrand <- ggplot(cateRandStat,aes(x = motif,y = meanNum,fill = sys,width = 0.9)) + geom_bar(stat = "identity",color = "black", position = position_dodge()) + geom_errorbar(aes(ymin=meanNum,ymax = meanNum + std), width = 0.2, position = position_dodge(.9))
gmotifrand <- gmotifrand + theme_classic() + scale_fill_manual(values = mycolor[c(1,10)]) + theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5)) + labs(x = "motif",y="counts")+ theme(axis.text = element_text(size=14), axis.title = element_text(size=18,face="bold"))



##====================== visualization the interaction network===============
## trime the data can be supported in igraph
## dimension of bipartite networks
NumOdors <- length(unique(as.vector(recordMatx[,c(1,2)])))
NumOr <- length(unique(as.vector(recordMatx[,c(3,4)])))

allOdorName <- rownames(newMatx)
allOrName <- colnames(newMatx)
# allType, the chemical type of these 107 odorants
## initialize the adjecent matrix
links <- matrix(NA,nrow = 107,ncol = 24)
colnames(links) <- colnames(newMatx)
rownames(links) <- rownames(newMatx)
for(i in 1:dim(recordMatx)[1]){
   oneRow <- recordMatx[i,]
   links[oneRow[1],oneRow[3]] <- oneRow[5]
   links[oneRow[1],oneRow[4]] <- oneRow[6]
   links[oneRow[2],oneRow[3]] <- oneRow[7]
   links[oneRow[2],oneRow[4]] <- oneRow[8]
}

#get rid of irrelavent rows and columns
irRow <- apply(links,1,function(x) all(is.na(x)))
irCol <- apply(links,2,function(x) all(is.na(x)))
links2 <- links[-which(irRow),-which(irCol)]
relType <- allType[-which(irRow)]  #chemical types

#final chemical types 
fType <- allType[-which(irRow)]

nodes <- data.frame(
      id = c(rownames(links2),colnames(links2)),
      cType = c(as.character(fType),rep(NA,times = dim(links2)[2])),
      category = c(rep(1,times=dim(links2)[1]),rep(2,times=dim(links2)[2]))
)

#use igraph to plot the bipartite network
net <- graph_from_incidence_matrix(links2)
table(V(net)$type)
#mycolor <- brewer.pal(11,"Spectral")
mycolor <- c(brewer.pal(9,"Set1"),"black")
colrs <- mycolor[c(1,10)]
V(net)$color <- mycolor[nodes$cType]
plot(net,vertex.label=NA,vertex.size = 4,layout=layout.bipartite)


tx <- links2
tx[is.na(tx)] <- 0
link3 <- as.data.frame(which(tx!=0,arr.ind = T))
odorID <- rownames(tx)[link3[,1]]
orID <- colnames(tx)[link3[,2]]
wt <- abs(tx[as.matrix(link3)])
edgeType <- tx[as.matrix(link3)]
edgeType[edgeType > 0] <- 1  #excitatory
edgeType[edgeType < 0] <- 2  #inhibitory
link3[] <- cbind(odorID,orID)
link3$width <- wt
link3$color <- c("grey","red")[edgeType]
colnames(link3) <- c("from","to","width","color")
#nodes3 <- c(unique(link3[,1]),unique(link3[,2]))
#nodes$id <- c(unique(link3[,1]),unique(link3[,2]))


nodes$label <- c(rownames(tx),colnames(tx))
nodes$shape <- c("square","circle")[nodes$category]
nodes$color.background <- mycolor[nodes$cType]
nodes$color.border <- "black"
nodes$color.highlight.background <-"orange"
nodes$color.highlight.border <- "darked"
nodes$size <- rep(20,times=length(nodes$shape))
#nodes$label <- rep(NA,times=length(nodes$shape))

visNetwork(nodes,link3)

visNetwork(nodes, link3) %>% visOptions(highlightNearest = TRUE,selectedBy = "cType")


#sub networks of odors elicit negative and positive response on at least three OSNs
# based on recordMatrx
nOdor <- table(recordMatx[,1])
bifanInd <- as.data.frame(recordMatx[,1:2])
repTable <- aggregate(list(numdup=rep(1,nrow(bifanInd))),bifanInd,length)
selecTable <- subset(repTable,numdup>3)
allSubNet <- c()
par(mfrow=c(4,3), mar=c(1,1,1,1))
for (i in 1:dim(selecTable)[1]){
    subTab <- recordMatx[recordMatx[,1] == selecTable[i,1] & recordMatx[,2] == selecTable[i,2],]
    links4 <- as.data.frame(unique(rbind(subTab[,c(1,3,5)],subTab[,c(1,4,6)],subTab[,c(2,3,7)],subTab[,c(2,4,8)])))
    links4[,1:2] <- cbind(rownames(adjMatx)[links4[,1]],colnames(adjMatx)[links4[,2]])
    colnames(links4) <- c("from","to","weight")
    links4$color <- c("red","grey")[as.factor(links4$weight>0)]
    links4$width <- abs(links4$weight)*1.5
    
    nodes4 <- data.frame(id=c(unique(links4[,1]),unique(links4[,2])),color=c("blue","red")[as.factor(c(1,1,rep(2,times = length(unique(as.vector(subTab[,3:4]))))))])
    nodes4$shape <- c("square","circle")[as.factor(nodes4$color)]
    nodes4$color <- c(mycolor[allType[subTab[1,1:2]]],rep("#E0E0E0",length(unique(as.vector(subTab[,3:4])))))
    nodes4$label <- nodes4$id
    nodes4$size <- rep(20,times=length(nodes4$id))
    allSubNet[[i]] <- visNetwork(nodes4,links4)
    #browser()
    rm(nodes4)
    rm(links4)
    rm(subTab)
    #
    
    #nodeName <- c(rownames(adjMatx)[subTab[1,1:2]],unique(links4[,2]))
    
    # bipartite
    # bmx <-acast(as.data.frame(links4),from ~ to)
    # bmx <- matrix(as.numeric(unlist(bmx)),nrow=nrow(bmx))
    # 
    # net5 <- graph_from_incidence_matrix(bmx)
    # V(net5)$label <- nodeName
    # plot(net5,layout=layout.bipartite)
}


#calculate the degree distribution based on links2 information 
dm <- links2
dm[is.na(dm)] <- 0
posiDegree<- apply(dm>0, 1, sum)
negaDegree <- apply(dm<0, 1, sum)
OrPosiDegr <- apply(dm>0, 2, sum)
OrNegaDegr <- apply(dm<0, 2, sum)
OdorDegreeFrame <- data.frame(
    odor= rep(factor(rownames(links2),levels = rownames(links2)),times = 2),
    type = factor(rep(c("positive","negative"),each=length(posiDegree))),
    deg = c(posiDegree,negaDegree),
    category = c(relType,relType)
)
degPlot <- ggplot(OdorDegreeFrame, aes(x=odor, y= deg,fill=type))+ geom_bar(stat="identity") + coord_flip()  #here use reorder to oder the bar "reorder(odor, -deg)"
degPlot <- degPlot + theme_classic() + theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5)) + labs(x = "odorant",y="degree")+ theme(axis.text.y = element_text(size=10,color = mycolor[factor(relType)]),axis.text.x = element_text(size=14), axis.title = element_text(size=18,face="bold"))
ggsave(degPlot,file = "degreeDistBifan.pdf",width = 6,height = 7)

# the degree of OSN
OSNDegreeFrame <- data.frame(
    OSN= rep(colnames(links2),times = 2),
    type = factor(rep(c("positive","negative"),each=length(OrPosiDegr))),
    deg = c(OrPosiDegr,OrNegaDegr)
)
degOrPlot <- ggplot(OSNDegreeFrame, aes(x= reorder(OSN,-deg), y= deg,fill=type))+ geom_bar(stat="identity") + coord_flip()  #here use reorder to oder the bar "reorder(odor, -deg)"
degOrPlot <- degOrPlot + theme_classic() + theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5)) + labs(x = "odorant",y="degree")+ theme(axis.text = element_text(size=14),axis.title = element_text(size=18,face="bold"))
ggsave(degOrPlot,file = "degreeOSNDistBifan.pdf",width = 4,height = 5)


#=============================================================
# some odorants elicits much strong responses, those at the tail of distribution
# sort the odorants that can elicits strong OSN responses
odorRespSumm <- apply(abs(adjMatx)>0, 1, sum)
odorSum <- data.frame(
        odor = factor(rownames(adjMatx)[-which(odorRespSumm<=9)],levels = rownames(adjMatx)[-which(odorRespSumm<=9)]),
        repNum = odorRespSumm[odorRespSumm>9],
        category = allType[-which(odorRespSumm<=9)]
)

mycolor <- c(brewer.pal(8,"Dark2"),"#A6CEE3","black")
odorSumPlot <- ggplot(odorSum, aes(x=odor, y= repNum))+ geom_bar(stat="identity") + coord_flip()  #here use reorder to oder the bar "reorder(odor, -deg)"
odorSumPlot <- odorSumPlot + theme_classic() + theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5)) + labs(x = "odorant",y="number OSN response")+ theme(axis.text.y = element_text(size=10,color = mycolor[odorSum$category]),axis.text.x = element_text(size=14), axis.title = element_text(size=18,face="bold"))
ggsave(odorSumPlot,file = "strongOdor.pdf",width = 5,height = 5)
# net <- graph_from_data_frame(d=link3, vertices=nodes, directed=T)
# plot(net,edge.arrow.size=.4,vertex.label=NA,edge.color = link3$color,vertex.color = nodes$color.background,layout=layout_randomly)
# #second method to generate a igraph
# links3 <- data.frame(from = integer(),to = integer(),type = charsacter(),weight = integer())
# for(i in 1:dim(recordMatx)[1]){
#     #oneRow <- recordMatx[i,]
#     links3$from[i] = recordMatx[i,1]
#     links3$to[i] = recordMatx[i,2]
#     if(recordMatx[i,3] >0){
#       links3$type[i] = "activation"
#     }
#     else{
#       links3$type[i] = "inhibition"
#     }
#     links3$weight[i] = recordMatx[i,3]
# }
# mycolor <- brewer.pal(11,"Spectral")
# posiInter <- newMatx
# posiInter[posiInter < 0] <- 0
# posiInter[posiInter > 0] <-1
# netP <- graph_from_incidence_matrix(posiInter)
# vcolors <- c(rep("gray",dim(posiInter)[1]),rep("red",dim(posiInter)[2]))
# V(netP)$color <- vcolors
# V(netP)$label <- NA
# V(netP)$size <- 3
# 
# plot(netP)
# plot(netP, vertex.label=NA, vertex.size=4, layout=layout.bipartite) 
# plot(netP, edge.arrow.mode=0, layout="layout_as_bipartite", main=layout)
# 
# allInter <- newMatx
# allInter[allInter!=0] <- 1
# plot(netA,layout=layout.bipartite,vertex.label=NA,vertex.size=4,vertex.color=c(mycolor[1],mycolor[10])[V(netA)$type+1],vertex.frame.color = NA)

# netDataFrame <- aggregate(allInter,)


#delete the ORNs that only have inhibitory response or no repoonse
# tempORN <- apply(apply(newMatx,2,range),2,sum)
# newMatx <-newMatx[,-which(tempORN==-1)]

# ##=======================example ===========================
# filePath <- "/Users/shan/Desktop/polnet2016/Data\ files/"
# fileName1 <- paste(filePath,"Dataset1-Media-Example-NODES.csv",sep = "")
# fileName2 <- paste(filePath,"Dataset1-Media-Example-EDGES.csv",sep = "")
# nodes <- read.csv(fileName1, header=TRUE, as.is=TRUE)
# links <- read.csv(fileName2, header=TRUE, as.is=TRUE)
# links <- aggregate(links[,3], links[,-3], sum)
# links <- links[order(links$from, links$to),]
# colnames(links)[4] <- "weight"
# rownames(links) <- NULL
# 
# ## another data set
# fileName3 <- paste(filePath,"Dataset2-Media-User-Example-NODES.csv",sep = "")
# fileName4 <- paste(filePath,"Dataset2-Media-User-Example-EDGES.csv",sep = "")
# nodes2 <- read.csv(fileName3, header=T, as.is=T)
# links2 <- read.csv(fileName4, header=T, row.names=1)
# links2 <- as.matrix(links2)
# 
# net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
# plot(net)  #plot, but not a pretty network
# net <- simplify(net, remove.multiple = F, remove.loops = T)
# plot(net, edge.arrow.size=.4,vertex.label=NA)
# 
# net2 <- graph_from_incidence_matrix(links2)
# 
# colrs <- c("gray50", "tomato", "gold")
# V(net)$color <- colrs[V(net)$media.type]
# 
# # Compute node degrees (#links) and use that to set node size:
# deg <- degree(net, mode="all")
# V(net)$size <- deg*3
# # We could also use the audience size value:
# V(net)$size <- V(net)$audience.size*0.6
# 
# # The labels are currently node IDs.
# # Setting them to NA will render no labels:
# V(net)$label <- NA
# 
# # Set edge width based on weight:
# E(net)$width <- E(net)$weight/6
# 
# #change arrow size and edge color:
# E(net)$arrow.size <- .2
# E(net)$edge.color <- "gray80"
# E(net)$width <- 1+E(net)$weight/12
# plot(net)
# #
# # layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1]
# # # Remove layouts that do not apply to our graph.
# # layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
# # par(mfrow=c(3,4), mar=c(1,1,1,1))
# # for (layout in layouts) {
# #     print(layout)
# #     l <- do.call(layout, list(net))
# #     plot(net, edge.arrow.mode=0, layout=l, main=layout) }
# #
# # plot(net2, vertex.label=NA, vertex.size=7, layout=layout.bipartite)
# 
# 
# #====== animation
# nodes$shape <- "dot"  
# nodes$shadow <- TRUE # Nodes will drop shadow
# nodes$title <- nodes$media # Text on click
# nodes$label <- nodes$type.label # Node label
# nodes$size <- nodes$audience.size # Node size
# nodes$borderWidth <- 2 # Node border width
# 
# nodes$color.background <- c("slategrey", "tomato", "gold")[nodes$media.type]
# nodes$color.border <- "black"
# nodes$color.highlight.background <- "orange"
# nodes$color.highlight.border <- "darkred"
# 
# visNetwork(nodes, links)
# 
# ##===============
# 
# # Random bipartite graph
# inc <- matrix(sample(0:1, 50, replace=TRUE, prob=c(2,1)), 10, 5)
# g <- graph.incidence(inc)
# plot(g, layout=layout.bipartite,
#      vertex.color=c(mycolor[1],mycolor[10])[V(g)$type+1])
# 
# # Two columns
# lay <- layout.bipartite(g)
# plot(g, layout=lay[,2:1])
