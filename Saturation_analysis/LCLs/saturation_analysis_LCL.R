#This is a script to determine whether I have saturated the analysis.
#Start with LCLs X-responsive genes
inputValue <- commandArgs(trailingOnly = TRUE)

# inputValue <- c(10,100)
print(as.character(inputValue[1]))
print(as.character(inputValue[2]))

set.seed(1)

#### 0. Setup ####

#PATH TO GITHUB FILES 
myPath <- #PATH TO GITHUB FILES 

suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("BiocParallel"))
register(MulticoreParam(16))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("rlist"))

annofile <-read.delim(file=paste0(myPath,"Annotations/annotation_with_ERCC.txt"),sep = " ")
tx2gene <-data.frame("TXNAME" = annofile$transcript_id,"GENEID" = annofile$gene_name)

#bring in expressed genes
load(file=paste0(myPath,"Linear_regressions/LCLs/lcl_expressed.rda"))

#### 1. Bring in metadata tables for aneuploidy samples ####

#read in the metadata table
metadata <-read.delim(paste0(myPath,"Linear_regressions/metadata.txt"), stringsAsFactors = FALSE)

#restrict to cell type and karyotype
metadata_lcl <- metadata[metadata$karyotype %in% c("X","XX","XXX","XXXX","XXXXY","XXXY","XXY","XXYY","XY","XYY","XYYYY") & metadata$cell_type == "LCL" & metadata$Technical_replicate == "N",]
rownames(metadata_lcl) <- metadata_lcl$sample


#### 2. Bring in tximport ####
load(file=paste0(myPath,"/Linear_regressions/LCLs/txi_lcl.rda"))

#### 3. Choose subset of the samples and run deseq #### 

#create results table
saturationResults <- NULL
resultSampleTables <- vector("list")

#choose this many samples 100 times
i <- as.numeric(inputValue[1])
k <- as.numeric(inputValue[2])

print(Sys.time())

j <- 1
repeat {
  #randomly select samples
  mySamples <-sample(1:dim(metadata_lcl)[1], i, replace = FALSE)
  sample_table <- metadata_lcl[mySamples,c("sample","x_count","y_count","batch_libprep")]
  
  #Test whether matrix is full rank
  mod_npx <- model.matrix( ~ batch_libprep + y_count + x_count, sample_table)
  if (is.fullrank(mod_npx)) {
    
    # #Test whether there are replicates
    if(any(duplicated(sample_table[,2:4]))){
      
      ##Restrict to chosen samples
      txi_samp <- txi_lcl
      txi_samp$abundance <- txi_samp$abundance[,sample_table$sample]
      txi_samp$counts <- txi_samp$counts[,sample_table$sample]
      txi_samp$length <- txi_samp$length[,sample_table$sample]

      ## Perform deseq
      dds <- DESeqDataSetFromTximport(txi = txi_samp, colData = sample_table, design = ~ batch_libprep + y_count + x_count)
      dds <- DESeq(dds, parallel = TRUE)
      
      #Get autosomal results
      x_results <- data.frame(results(dds, name="x_count", alpha=0.05))
      x_results_auto <- na.omit(x_results[rownames(x_results) %in% lcl_expressed$expAutoGenes,])
      x_results_auto_sig <- data.frame(x_results_auto[x_results_auto$padj < 0.05,])
      x_results_auto_sig_up <- x_results_auto_sig[x_results_auto_sig$log2FoldChange > 0,]
      x_results_auto_sig_down <- x_results_auto_sig[x_results_auto_sig$log2FoldChange < 0,]
      
      #Do analysis of Y genes
      y_results <- data.frame(results(dds, name="y_count", alpha=0.05))
      y_results_auto <- na.omit(y_results[rownames(y_results) %in% lcl_expressed$expAutoGenes,])
      y_results_auto_sig <- y_results_auto[y_results_auto$padj < 0.05,]
      y_results_auto_sig_up <- y_results_auto_sig[y_results_auto_sig$log2FoldChange > 0,]
      y_results_auto_sig_down <- y_results_auto_sig[y_results_auto_sig$log2FoldChange < 0,]
      
      #Record the number of significant X or Y responsive genes in a table
      sigGenes <- c(dim(sample_table)[1], 
                    dim(x_results_auto_sig)[1], dim(x_results_auto_sig_up)[1], dim(x_results_auto_sig_down)[1],
                    dim(y_results_auto_sig)[1], dim(y_results_auto_sig_up)[1], dim(y_results_auto_sig_down)[1])
      names(sigGenes) <- c("Sample_num", 
                           "x_auto_sig","x_auto_up","x_auto_down",
                           "y_auto_sig","y_auto_up","y_auto_down")
      saturationResults <- rbind(saturationResults, sigGenes)
      
      resultSampleTables[paste0("round_",as.character(i),"_",as.character(j))] <- list(sample_table)
      
      #print the results out
      print(paste0("Finished round ",as.character(i),",",as.character(j),". # of X responsive auto genes: ", as.character(dim(x_results_auto_sig)[1]),"; # of Y-responsive auto genes: ",as.character(dim(y_results_auto_sig)[1])))
      
      #advance j
      j <- j + 1
      
      print(j)
      
    } else{
      print("No replicates, drawing samples again.")
      # print(sample_table)
    }
    
  } else{
    print("Model not full rank, drawing samples again.")
  }
  
  if(j == (k + 1)){
    break
  }
  
}


##outside for loop
saturationResults <- as.data.frame(saturationResults)
colnames(saturationResults) <- c("Sample_num", 
                                 "x_auto_sig","x_auto_up","x_auto_down",
                                 "y_auto_sig","y_auto_up","y_auto_down")
write.table(saturationResults, file=paste0(myPath,"Saturation_analysis/LCLs/lcl_saturationResults_100iterations_", as.character(i),".txt"), quote=FALSE, col.names = TRUE, row.names = FALSE, sep="\t")

list.save(resultSampleTables,file = paste0(myPath,"Saturation_analysis/LCLs/lcl_resultSampleTables_100iterations_",as.character(i),".rdata"))

Sys.time()

sessionInfo()

