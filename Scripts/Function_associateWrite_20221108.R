suppressPackageStartupMessages(library(dplyr))
associateWrite <- function(myOrderedFile){
  myOrderedFile_df <- as.data.frame(myOrderedFile)
  myOrderedFile_anno <- cbind(myOrderedFile_df, geneAnno[match(rownames(myOrderedFile_df), geneAnno$Gene),])
}