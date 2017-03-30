visualBifan <-function(recordMatx,allOdorName,allOrName,allType,visualNet = FALSE)
{
    # ================================================
    # this fucntion visualize the several strong bifan motifs
    # recordMatx      a two column matrix indicates the interaction of bifan motif
    # allOdorName     all the names of odorants
    # allOrName       all the OSN names
    # allType         all the chemical categories of the odorants
    # visualNet       whethe plot the visualized graph
    # =================================================
    
    # dependent libraries
    
    ## trime the data can be supported in igraph
    ## dimension of bipartite networks
    NumOdors <- length(unique(as.vector(recordMatx[,c(1,2)])))
    NumOr <- length(unique(as.vector(recordMatx[,c(3,4)])))
    #browser()
    
    #allOdorName <- rownames(newMatx)
    #allOrName <- colnames(newMatx)
        
    # allType, the chemical type of these 107 odorants
    ## initialize the adjecent matrix
    links <- matrix(NA,nrow = 107,ncol = 24)
    # colnames(links) <- colnames(newMatx)
    # rownames(links) <- rownames(newMatx)
    colnames(links) <- allOrName
    rownames(links) <- allOdorName
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
    #plot(net,vertex.label=NA,vertex.size = 4,layout=layout.bipartite)
    
    
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
    
    if(visualNet){
        visNetwork(nodes,link3)
        # dynamical interactive graph
        #visNetwork(nodes, link3) %>% visOptions(highlightNearest = TRUE,selectedBy = "cType")
    }
    
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
    }
    
    #calculate the degree distribution based on links2 information 
    colorOrder <- sort(unique(allType),index.return = T)
    #mycolor <- c(brewer.pal(8,"Dark2"),"#A6CEE3","black")
    mycolor2 <- colors()[c(567,566,456,24,88,225,142,121,136,614)]
    mycolor2 <- mycolor2[colorOrder$ix]
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
    
    #figure folder
    figFolder <- "/Users/shan/Documents/GoogleDrive/olfactoryCoding/figures/motif/"
    degPlot <- ggplot(OdorDegreeFrame, aes(x=odor, y= deg,fill=type))+ geom_bar(stat="identity") + coord_flip()  #here use reorder to oder the bar "reorder(odor, -deg)"
    degPlot <- degPlot + theme_classic() + theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5)) + labs(x = "odorant",y="degree")+ theme(axis.text.y = element_text(size=10,color = mycolor2[relType]),axis.text.x = element_text(size=14), axis.title = element_text(size=18,face="bold"))
    ggsave(degPlot,file = paste(figFolder,"degreeDistBifan.pdf",sep = ""),width = 6,height = 7)
    
    # the degree of OSN
    OSNDegreeFrame <- data.frame(
        OSN= rep(colnames(links2),times = 2),
        type = factor(rep(c("positive","negative"),each=length(OrPosiDegr))),
        deg = c(OrPosiDegr,OrNegaDegr)
    )
    degOrPlot <- ggplot(OSNDegreeFrame, aes(x= reorder(OSN,-deg), y= deg,fill=type))+ geom_bar(stat="identity") + coord_flip()  #here use reorder to oder the bar "reorder(odor, -deg)"
    degOrPlot <- degOrPlot + theme_classic() + theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5)) + labs(x = "odorant",y="degree")+ theme(axis.text = element_text(size=14),axis.title = element_text(size=18,face="bold"))
    ggsave(degOrPlot,file = paste(figFolder,"degreeOSNDistBifan.pdf"),width = 4,height = 5)
    
    

}