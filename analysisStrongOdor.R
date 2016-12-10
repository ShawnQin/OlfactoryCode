analysisStrongOdor <- function(adjMatx,allTypE)
{  
    # =================================================
    # this fucntion compare and plot the most strongest odorants
    # use a chemical metric to compare their 
    # return a list of three types of strong odorants
    
    # adjMatx       an interaction matrix
    # allType       all the chemical categories
    # =================================================
    
    # set colors
    colorOrder <- sort(unique(allType),index.return = T)
    #mycolor <- c(brewer.pal(8,"Dark2"),"#A6CEE3","black")
    mycolor2 <- colors()[c(567,566,456,24,88,225,142,121,136,614)]
    mycolor2 <- mycolor2[colorOrder$ix]
    
    # some odorants elicits much stronger responses, those at the tail of distribution
    # sort the odorants that can elicits strong OSN responses
    
    figFolder <- "/Users/shan/Documents/GoogleDrive/olfactoryCoding/figures/motif/"
    odorRespSumm <- apply(abs(adjMatx)>0, 1, sum)
    odorSum <- data.frame(
        odor = factor(rownames(adjMatx)[-which(odorRespSumm<=9)],levels = rownames(adjMatx)[-which(odorRespSumm<=9)]),
        repNum = odorRespSumm[odorRespSumm>9],
        category = allType[-which(odorRespSumm<=9)]
    )

    odorSumPlot <- ggplot(odorSum, aes(x=odor, y= repNum))+ geom_bar(stat="identity") + coord_flip()  #here use reorder to oder the bar "reorder(odor, -deg)"
    odorSumPlot <- odorSumPlot + theme_classic() + theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5)) + labs(x = "odorant",y="number OSN response")+ theme(axis.text.y = element_text(size=10,color = mycolor2[odorSum$category]),axis.text.x = element_text(size=14), axis.title = element_text(size=18,face="bold"))
    ggsave(odorSumPlot,file = paste(figFolder,"strongOdor.pdf",sep = ""),width = 5,height = 5)
    
    ## strong positive responses odorants
    odorPosiSumm <- apply(adjMatx>0, 1, sum)
    posiThres <- 9     #threshold of positive odorant
    posiSum <- data.frame(
        odor = factor(rownames(adjMatx)[-which(odorPosiSumm<=posiThres)],levels = rownames(adjMatx)[-which(odorPosiSumm<=posiThres)]),
        repNum = odorPosiSumm[odorPosiSumm>posiThres],
        category = allType[-which(odorPosiSumm<=posiThres)]
    )
    odorPosiPlot <- ggplot(posiSum, aes(x=odor, y= repNum))+ geom_bar(stat="identity") + coord_flip()  #here use reorder to oder the bar "reorder(odor, -deg)"
    odorPosiPlot <- odorPosiPlot + theme_classic() + theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5)) + labs(x = "odorant",y="number OSN response")+ theme(axis.text.y = element_text(size=12,color = mycolor2[posiSum$category]),axis.text.x = element_text(size=14), axis.title = element_text(size=18,face="bold"))
    ggsave(odorPosiPlot,file = paste(figFolder,"strongPosiOdor.pdf"),width = 5,height = 4)
    
    
    ## inhibitory odorants, some odorants elicit more inhibiotry responses
    odorInhiSumm <- apply(adjMatx<0, 1, sum)
    threshold <- 4   #threhold number of inhibitory responses 
    inhiSum <- data.frame(
        odor = factor(rownames(adjMatx)[-which(odorInhiSumm<=threshold)],levels = rownames(adjMatx)[-which(odorInhiSumm<=threshold)]),
        repNum = odorInhiSumm[odorInhiSumm>threshold],
        category = allType[-which(odorInhiSumm<=threshold)]
    )
    odorInhiPlot <- ggplot(inhiSum, aes(x=odor, y= repNum))+ geom_bar(stat="identity") + coord_flip()  #here use reorder to oder the bar "reorder(odor, -deg)"
    odorInhiPlot <- odorInhiPlot + theme_classic() + theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5)) + labs(x = "odorant",y="number OSN response")+ theme(axis.text.y = element_text(size=12,color = mycolor2[inhiSum$category]),axis.text.x = element_text(size=14), axis.title = element_text(size=18,face="bold"))
    ggsave(odorInhiPlot,file = paste(figFolder,"strongInhiOdor.pdf"),width = 5,height = 4)
    
    strongOdor <- list(odorSum = odorSum,posiSum = posiSum,inhiSum = inhiSum)
    return(strongOdor)
    
}