allResults <- function(mydds,myVar,myAnno,expGenes,myalpha = 0.05, myTitle, 
                       myYlim = c(-1,1), pdfWidth = 10, pdfHeight=7, highlight_col = c("#00000080","#A9A9A980"), bg_col = c("#00000000","#00000000")) {
  #A function to get all of the results from DESeq2 output.
  #mydds <- deseq object (name of dds object)
  #myVar <- which variable you want to get the results for (character)
  #myAnno <- data file for the annotated genes (geneAnno)
  #myalpha <- adjusted p-value cutoff (numeric)
  #myTitle <- what you want the filenames to start with for the output results (Character)
  
  # mydds <- dds_x_y_lcl
  # myVar <- "x_count"
  # myAnno <- geneAnno
  # expGenes <- lcl_expressed
  # myalpha = 0.05
  # myTitle = "lcl_x_count"
  # pdfWidth = 10
  # pdfHeight = 5
  # myYlim = c(-1.5,1.5)
  
  #get results in my annotation
  res <- results(mydds, name = myVar, alpha = myalpha)
  res <- res[rownames(res) %in% myAnno$Gene, ]
  summary(res) 
  res_exp <- res[rownames(res) %in% expGenes$expressedGenes, ]
  res_exp_ordered <- res_exp[order(res_exp$padj), ]
  
  #associate with annotation and write out the expressed genes
  res_anno_exp <- merge(x=as.data.frame(res_exp_ordered), y=geneAnno, by.x=0, by.y="Gene") #Gene v84
  res_anno_exp <- res_anno_exp[! duplicated(res_anno_exp$Row.names),]
  rownames(res_anno_exp) <- res_anno_exp$Row.names
  colnames(res_anno_exp)[1] <- "Gene"
  write.table(res_anno_exp, file = paste0(myTitle, "_expressed.txt"), quote = FALSE, sep = "\t")
  
  res_anno_exp_auto <- res_anno_exp[rownames(res_anno_exp) %in% expGenes$expAutoGenes,]
  write.table(res_anno_exp_auto, file = paste0(myTitle, "_auto_expressed.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

  #write out significant list
  res_anno_exp_sig <- na.omit(res_anno_exp[res_anno_exp$padj < myalpha, ])
  
  res_anno_exp_sig_auto <- res_anno_exp_sig[rownames(res_anno_exp_sig) %in% expGenes$expAutoGenes,]
  # write.table(res_anno_exp_sig_auto, file = paste0(myTitle, "_p_", as.character(myalpha), "_auto.txt"),quote = FALSE,sep = "\t")
  
  #write out autosomal significant genes list
  res_anno_exp_sig_auto <-na.omit(res_anno_exp_sig[!res_anno_exp_sig$chr %in% c("chrX", "chrY"), ])
  write.table(rownames(res_anno_exp_sig_auto),file = paste0(myTitle,"_p_", as.character(myalpha),"_auto_genes.txt"),
              quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
  
  myRowNames <- c("All_genes", "Auto_genes")
  sigGenes <- dim(res_anno_exp_sig)[1]
  sigGenes_up <- dim(res_anno_exp_sig[res_anno_exp_sig$log2FoldChange > 0,])[1]
  sigGenes_down <- dim(res_anno_exp_sig[res_anno_exp_sig$log2FoldChange < 0,])[1]
  sig_autoGenes <- dim(res_anno_exp_sig_auto)[1]
  sig_autoGenes_up <- dim(res_anno_exp_sig_auto[res_anno_exp_sig_auto$log2FoldChange > 0,])[1]
  sig_autoGenes_down <- dim(res_anno_exp_sig_auto[res_anno_exp_sig_auto$log2FoldChange < 0,])[1]
  
  allGenes <- c(sigGenes, sig_autoGenes)
  upGenes <- c(sigGenes_up, sig_autoGenes_up)
  downGenes <- c(sigGenes_down, sig_autoGenes_down)
  myExpGenes <- c(length(expGenes$expressedGenes), length(expGenes$expAutoGenes))
  
  results_table <- data.frame("Sig_genes" = allGenes, "Up" = upGenes, "Down" = downGenes, "Expressed" = myExpGenes)
  row.names(results_table) <- myRowNames
  
  print(results_table)
  percentGenes <- 100 * ( sigGenes / length(expGenes$expressedGenes ))
  print(paste0("% of responsive expressed genes = ", format( round(percentGenes,2), nsmall = 2)))
  
  print("Top 50 significant autosome genes (by p-val):")
  sig_auto <- res_anno_exp_sig_auto[order(res_anno_exp_sig_auto$padj, decreasing = FALSE),]
  sig_auto50 <- sig_auto[1:50,]
  sig_auto50$padj <- format(sig_auto50$padj, scientific=TRUE)
  sig_auto50$log2FoldChange <- round(sig_auto50$log2FoldChange, 3)
  print(sig_auto50[,c("log2FoldChange", "padj","chr")])
  
  DESeqManhattan_noM_col_auto(res_anno_exp_auto, mySigGenes = res_anno_exp_sig_auto,
                              myTitle=myTitle,fileName=paste0(myTitle,"_p_",as.character(myalpha),"_manhattan_auto_noYlim.pdf"),
                              pointsize = 0.35,axissize = 0.5, labelsize = 0.5, ylabel = paste0("Log2 fold-change per ",myVar), pdfWidth = pdfWidth,
                              pdfHeight = pdfHeight, highlight_col = highlight_col, reg_col = bg_col)
  
  DESeqManhattan_noM_col_auto(res_anno_exp_auto, mySigGenes = res_anno_exp_sig_auto,
                              myTitle=myTitle,fileName=paste0(myTitle,"_p_",as.character(myalpha),"_manhattan_auto_Ylim.pdf"),
                              pointsize = 0.35,axissize = 0.5, labelsize = 0.5, ylabel = paste0("Log2 fold-change per ",myVar), myYlim = myYlim,pdfWidth = pdfWidth,
                              pdfHeight = pdfHeight, highlight_col = highlight_col, reg_col = bg_col)
   
  normalized_counts <- counts(mydds,normalized=TRUE)
  write.table(normalized_counts,file = paste0(myTitle,"_normcounts.txt"),
              quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
  
  log_norm_counts <- log2(counts(mydds,normalized=TRUE)+1)
  
    myList <- list("res_anno_exp" = res_anno_exp, "res_anno_exp_auto" = res_anno_exp_auto, "res_anno_exp_sig" = res_anno_exp_sig,
                   "res_anno_exp_sig_auto" = res_anno_exp_sig_auto, "norm_counts" = 
                     normalized_counts, "log_norm_counts" = log_norm_counts)

  return(myList)
}


