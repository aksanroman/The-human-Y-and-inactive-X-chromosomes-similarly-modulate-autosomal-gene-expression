tpm1_genes_xx_xy_median <- function(tximport_file, XX_samples, XY_samples) {
    kallisto_output_TPM_XY <- tximport_file$abundance[,c(XY_samples)]
    kallisto_output_TPM_XX <- tximport_file$abundance[,c(XX_samples)]
    median_TPM_XY <- apply(kallisto_output_TPM_XY,1,median)
    median_TPM_XX <- apply(kallisto_output_TPM_XX,1,median)
    #Put all of the medians together in a dataframe of expression per karyotype - numeric
    expression_perKaryotype <- data.frame(median_TPM_XY, median_TPM_XX)
    #Isolate the rows of the expression_TPM dataset that have at least one karyotype passing the threshold
    expressedGenes <- expression_perKaryotype[expression_perKaryotype$median_TPM_XY >= 1 | expression_perKaryotype$median_TPM_XX >= 1,]
    #Restrict by gene annotation file
    expressedGenes <- expressedGenes[rownames(expressedGenes) %in% geneAnno$Gene,]

    #See how many genes are expressed in at least one karyotype
    num_expressedGenes <- dim(expressedGenes)[1]
    print(paste("Number of 'expressed' genes with median TPM >= 1 in either XX or XY samples: ", as.character(num_expressedGenes), sep=""))

    #Find expressed sex chromosome genes
    expSexChromGenes <- expressedGenes[rownames(expressedGenes) %in% sexChrom_genes,]
    print(paste("Number of 'expressed' sex chromosome genes: ", as.character(dim(expSexChromGenes)[1]), sep=""))

    #Find expressed Y chromosome genes
    expYGenes <- expressedGenes[rownames(expressedGenes) %in% Y_genes_all,]
    print(paste("Number of 'expressed' Y chromosome genes: ", as.character(dim(expYGenes)[1]), sep=""))

    #Find expressed X chromosome genes
    expXGenes <- expressedGenes[rownames(expressedGenes) %in% X_genes_all,]
    print(paste("Number of 'expressed' X chromosome genes: ", as.character(dim(expXGenes)[1]), sep=""))

    #Find expressed PAR genes
    expPARGenes <- expressedGenes[rownames(expressedGenes) %in% PAR_genes_all,]
    print(paste("Number of 'expressed' PAR genes (small set): ", as.character(dim(expPARGenes)[1]), sep=""))

    #Find expressed autosome genes
    expAutoGenes <- expressedGenes[rownames(expressedGenes) %in% autosome_genes,]
    print(paste("Number of 'expressed' autosome genes: ", as.character(dim(expAutoGenes)[1]), sep=""))

    myList <- list("expressedGenes" = rownames(expressedGenes), "expSexChromGenes" = rownames(expSexChromGenes),
                   "expYGenes" = rownames(expYGenes), "expXGenes" = rownames(expXGenes), "expPARGenes" = rownames(expPARGenes),
                   "expAutoGenes" = rownames(expAutoGenes))

    return(myList)
  }
