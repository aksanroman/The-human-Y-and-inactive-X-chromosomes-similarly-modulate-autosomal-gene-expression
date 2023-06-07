suppressPackageStartupMessages(library(dplyr))
associateWrite <- function(myOrderedFile){
  myOrderedFile_df <- as.data.frame(myOrderedFile)
  myOrderedFile_anno <- cbind(myOrderedFile_df, geneAnno[match(rownames(myOrderedFile_df), geneAnno$Gene),])
}

deseqFactorAnalysis <- function(deseqObject, myVariable, comparisonList, expressedGeneList, p_value=0.05, myTitle="title"){
  for (i in 1:dim(comparisonList)[1]){
    print(paste("Data for the ",comparisonList[i,1]," vs ",comparisonList[i,2]," comparison."))
    #get results for the comparison we want to make
    #register(MulticoreParam(16))
    #comparison_result <- results(deseqObject, contrast=c(myVariable,comparisonList[i,1], comparisonList[i,2]), parallel = TRUE)
    comparison_result <- results(deseqObject, contrast=c(myVariable,comparisonList[i,1], comparisonList[i,2]))
    comparison_result <- comparison_result[rownames(comparison_result) %in% geneAnno$Gene,]
    #get expressed results
    comparison_result_exp <- comparison_result[rownames(comparison_result) %in% expressedGeneList,]
    
    #order the results
    comparison_exp_ordered <- comparison_result_exp[order(comparison_result_exp$padj),]
    head(comparison_exp_ordered)
    summary(comparison_exp_ordered)
    
    # print("Making MA plot")
    # #make MA plot
    # pdf(file=paste(myTitle,"_",comparisonList[i,1],"vs",comparisonList[i,2],".pdf",sep=""))
    # plotMA(comparison_exp_ordered, main=paste(comparisonList[i,1],"vs",comparisonList[i,2],sep=" "), ylim=c(-10,4))
    # dev.off()
    # 
    # print("Making volcano plot")
    # #make volcano plots
    # DESeqVolcano(comparison_result_exp,sig_level=p_value,fileName=paste(myTitle,"_",comparisonList[i,1],"vs",comparisonList[i,2],"_",as.character(p_value),"_volcano.pdf",sep=""))
    # dev.off()
    # DESeqVolcano(comparison_result_exp,sig_level=0.0001, myLabel = TRUE, fileName=paste(myTitle,"_",comparisonList[i,1],"vs",comparisonList[i,2],"_0.0001_volcano.pdf",sep=""))
    # dev.off()
    # 
    print("Annotate files")
    
    
    #annotate and write out files
    comparison_exp_anno <- associateWrite(comparison_exp_ordered)
    write.table(comparison_exp_anno, file=paste(myTitle,"_",comparisonList[i,1],"vs",comparisonList[i,2],"_DESeq_expressed.txt",sep=""),quote=FALSE,sep="\t")
    comparison_anno_exp_sig <- comparison_exp_anno[comparison_exp_anno$padj<p_value,]
    comparison_anno_exp_sig <- na.omit(comparison_anno_exp_sig)
    write.table(comparison_anno_exp_sig, file=paste(myTitle,"_",comparisonList[i,1],"vs",comparisonList[i,2],"_DESeq_exp_p",as.character(p_value),".txt",sep=""), quote=FALSE, sep="\t")
    
    #make manhattan plots 
    DESeqManhattan_noM_col(comparison_exp_anno, mySigGenes = comparison_anno_exp_sig, myTitle=paste(comparisonList[i,1],"vs",comparisonList[i,2],"padj<",as.character(p_value),sep=" "),fileName=paste(myTitle,"_",comparisonList[i,1],"vs",comparisonList[i,2],"_",as.character(p_value),"_manhattan.pdf",sep=""))
  }
}