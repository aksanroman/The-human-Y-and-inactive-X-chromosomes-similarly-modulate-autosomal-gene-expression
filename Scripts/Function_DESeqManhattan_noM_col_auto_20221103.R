library(qqman)
DESeqManhattan_noM_col_auto <- function(myAnnoData, mySigGenes=NULL, myTitle="Manhattan Plot", fileName="manhattan.pdf",myYlim=NULL,pngSize = c(5,5), pointsize=1,axissize = 1, labelsize=1,ylabel="Slope of regression (log2)", pdfWidth=7, 
                                        pdfHeight=5, highlight_col = c("#0072B290","#D55E0090"), reg_col = c("#99999920", "#99999920")){
  
  #replace PAR chromosomes
  my.chrom <- as.character(myAnnoData$chr)
  
  #convert chromosomes to numbers, remove M, remove X, remove Y
  chromosomes = c(1:22)
  for(myChr in chromosomes){
    my.chrom <- replace(my.chrom, my.chrom==paste0("chr",as.character(myChr)), myChr)
  }
  
  #get the dataframe ready
  my.df <- data.frame(rownames(myAnnoData),as.numeric(my.chrom),myAnnoData$start, myAnnoData$padj, myAnnoData$log2FoldChange)
  colnames(my.df) <- c("SNP","CHR", "BP","P","log2FC")
  my.df <- na.omit(my.df)
  my.df <- my.df[! my.df$CHR %in% c("chrM","chrX","chrY"),]
  
  #make the manhattan plot
  pdf(file=fileName, width=pdfWidth, height=pdfHeight)
  
  if(is.null(mySigGenes)){
    ASRQQmanhattan(my.df, p="log2FC", logp=FALSE, ylab=ylabel,main=myTitle, cex=pointsize,
                   cex.axis=axissize, col=reg_col, chrlabs=as.character(c(1:22)), cex.lab=labelsize, 
                   myLim = myYlim)
  } else{
    
    #highlight significant genes
    sig_genes <- rownames(mySigGenes)
    ASRQQmanhattan(my.df, p="log2FC", logp=FALSE, ylab=ylabel,main=myTitle, cex=pointsize,
                   cex.axis=axissize, col=reg_col, colhighlight= highlight_col,
                   chrlabs=as.character(c(1:22)), highlight = sig_genes, cex.lab=labelsize, myLim=myYlim)
  }
  dev.off()
}