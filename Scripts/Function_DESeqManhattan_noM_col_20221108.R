suppressPackageStartupMessages(library(qqman))
DESeqManhattan_noM_col <- function(myAnnoData, mySigGenes=NULL, myTitle="Manhattan Plot", fileName="manhattan.pdf",myYlim=NULL,pngSize = c(5,5), pointsize=1,axissize = 1, labelsize=1,ylabel="Slope of regression (log2)", isFemale = FALSE,pdfWidth=7, pdfHeight=5){
  #replace PAR chromosomes
  my.chrom <- as.character(myAnnoData$chr)
  for(i in 1:(dim(myAnnoData)[1])){
    if(rownames(myAnnoData[i,]) %in% PAR_genes_all){
      my.chrom[i] <- "PAR"
    }
  }
  
  if(isFemale == TRUE){
    #convert chromosomes to numbers, remove M, set X as 23, PAR as 24
    chromosomes = c(1:22)
    for(myChr in chromosomes){
      my.chrom <- replace(my.chrom, my.chrom==paste("chr",as.character(myChr),sep=""), myChr)
    }
    my.chrom <- replace(my.chrom, my.chrom=="chrX", "23")
    my.chrom <- replace(my.chrom, my.chrom=="PAR", "24")
    #get the dataframe ready
    my.df <- data.frame(rownames(myAnnoData),as.numeric(my.chrom),myAnnoData$start, myAnnoData$padj, myAnnoData$log2FoldChange)
    colnames(my.df) <- c("SNP","CHR", "BP","P","log2FC")
    my.df <- na.omit(my.df)
    my.df <- my.df[my.df$CHR != "chrM",]
    
    #make the manhattan plot
    pdf(file=fileName, width=pdfWidth, height=pdfHeight)
    
    if(is.null(mySigGenes)){
      ASRQQmanhattan(my.df, p="log2FC", logp=FALSE, ylab=ylabel,main=myTitle, cex=pointsize,
                     cex.axis=axissize, col=c("#99999920", "#99999920"), chrlabs=c(1:22,"X","PAR"), cex.lab=labelsize, 
                     myLim = myYlim)
    }
    else{
      
      #highlight significant genes
      sig_genes <- rownames(mySigGenes)
      ASRQQmanhattan(my.df, p="log2FC", logp=FALSE, ylab=ylabel,main=myTitle, cex=pointsize,
                     cex.axis=axissize, col=c("#99999920", "#99999920"), colhighlight=c("#0072B2","#D55E00"),
                     chrlabs=c(1:22,"X","PAR"), highlight = sig_genes, cex.lab=labelsize, myLim=myYlim)
    }
    dev.off()
  }
  
  if(isFemale == FALSE){
    #convert chromosomes to numbers, remove M, set X as 23, Y as 24, PAR as 25
    chromosomes = c(1:22)
    for(myChr in chromosomes){
      my.chrom <- replace(my.chrom, my.chrom==paste("chr",as.character(myChr),sep=""), myChr)
    }
    my.chrom <- replace(my.chrom, my.chrom=="chrX", "23")
    my.chrom <- replace(my.chrom, my.chrom=="chrY", "24")
    my.chrom <- replace(my.chrom, my.chrom=="PAR", "25")
    #get the dataframe ready
    my.df <- data.frame(rownames(myAnnoData),as.numeric(my.chrom),myAnnoData$start, myAnnoData$padj, myAnnoData$log2FoldChange)
    colnames(my.df) <- c("SNP","CHR", "BP","P","log2FC")
    my.df <- na.omit(my.df)
    my.df <- my.df[my.df$CHR != "chrM",]
    
    #make the manhattan plot
    pdf(file=fileName, width=pdfWidth, height=pdfHeight)
    
    if(is.null(mySigGenes)){
      ASRQQmanhattan(my.df, p="log2FC", logp=FALSE, ylab=ylabel,main=myTitle, cex=pointsize,
                     cex.axis=axissize, col=c("#99999920", "#99999920"), chrlabs=c(1:22,"X","Y","PAR"), cex.lab=labelsize, 
                     myLim = myYlim)
    }
    else{
      
      #highlight significant genes
      sig_genes <- rownames(mySigGenes)
      ASRQQmanhattan(my.df, p="log2FC", logp=FALSE, ylab=ylabel,main=myTitle, cex=pointsize,
                     cex.axis=axissize, col=c("#99999920", "#99999920"), colhighlight=c("#0072B2","#D55E00"),
                     chrlabs=c(1:22,"X","Y","PAR"), highlight = sig_genes, cex.lab=labelsize, myLim=myYlim)
    }
    dev.off()
  }
}