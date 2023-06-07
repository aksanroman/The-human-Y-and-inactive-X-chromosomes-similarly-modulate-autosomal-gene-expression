# Load required packages
# library(data.table)

#Define functions
get_subjectID <- function(myFile){
  myResult <- strsplit(myFile,"/")[[1]][9]
  return(myResult)
}

get_haplotype <- function(myFile){
  myResult <- strsplit(myFile,"/")[[1]][8]
  return(myResult)
}

get_sampID <- function(myFile){
  myResult <- strsplit(myFile,"/")[[1]][10]
  return(myResult)
}

#Bring in list of file paths:

bedgraph_paths <- read.delim(file="#PATH TO GC-CORRECTED CHR Y BED FILES FROM 1000 GENOMES", stringsAsFactors = FALSE, header = FALSE)$V1
#GC-corrected bed files obtained from Teitz et al 2018 American Journal of Human Genetics "Selection Has Countered High Mutability to Preserve the Ancestral Copy Number of Y Chromosome Amplicons in Diverse Human Lineages."

#Make dataframe to add normalized DYZ1 coverage to:
results_table <- data.frame("path" = bedgraph_paths, "subjectID" = sapply(X = bedgraph_paths, FUN = get_subjectID), "sampleID" = sapply(X = bedgraph_paths, FUN = get_sampID), 
                            "Yhaplotype" = sapply(X = bedgraph_paths, FUN = get_haplotype),"Single_copy_depth" = numeric(length(bedgraph_paths)), 
                            "DYZ1_depth" = numeric(length(bedgraph_paths)),"DYZ1_norm_coverage" = numeric(length(bedgraph_paths)), 
                            "Num_samples" = numeric(length(bedgraph_paths)), stringsAsFactors = FALSE)

# Define regions of interest
dyz1_region <- c(56676260, 56771430)  # start position, end position
cntrl_region <- c(14500000, 15500000)

#List of subjects with duplicate samples
duplicate_samples <- results_table[duplicated(results_table$subjectID),"subjectID"] #215 of them

#Make new results table for the combined results for duplicated samples
results_table_unique <- results_table[! duplicated(results_table$subjectID),c(-1,-3)]
row.names(results_table_unique) <- results_table_unique$subjectID
results_table_dups <- results_table_unique[results_table_unique$subjectID %in% duplicate_samples,] #189 subjects with >1 sample!

#Go through all of the 1225 unique subjects!
for(mySubject in results_table_unique$subjectID){
  cat(paste0("Processing subject: ", mySubject,"\n"))
  
  #Get list of sample IDs corresponding to each subject
  results_table_mySubject <- results_table[results_table$subjectID == mySubject,]
  
  #Initialize results list for each subject
  results_list_control <- list()
  results_list_dyz1 <- list()
  
  #Go through all samples for that subject
  for(mySample in results_table_mySubject$sampleID){
    cat(paste0("Processing sample: ", mySample,"\n"))
    
    #Get list of files corresponding to the sample
    bedgraph_file <- results_table[results_table$sampleID == mySample,"path"]
    
    # Read bedgraph file into data table
    bedgraph <- scan(bedgraph_file, what = list("", 0, 0, 0), sep="\t")
    bedgraph <- data.frame(chrom = bedgraph[[1]], start = bedgraph[[2]], end = bedgraph[[3]], coverage = bedgraph[[4]])
    
    #get part of bedgraph in regions of interest
    dyz1_bed <- bedgraph[bedgraph$start >= dyz1_region[1] & bedgraph$end < dyz1_region[2], ]
    cntrl_bed <- bedgraph[bedgraph$start >= cntrl_region[1] & bedgraph$end < cntrl_region[2], ]
    
    #get part of bedgraph that represents single or multiple base pairs
    dyz1_bed_single <- dyz1_bed[(dyz1_bed$end - dyz1_bed$start) == 1, ]
    dyz1_bed_multiple <- dyz1_bed[(dyz1_bed$end - dyz1_bed$start) > 1, ]
    #how many new rows to make for dyz1?
    dyz1_bed_multiple_calc <- sum(dyz1_bed_multiple$end - dyz1_bed_multiple$start )
    dyz1.cov <- numeric(length = dyz1_bed_multiple_calc)
    dyz1.start  <- numeric(length = dyz1_bed_multiple_calc)
    
    cntrl_bed_single <- cntrl_bed[(cntrl_bed$end - cntrl_bed$start) == 1, ]
    cntrl_bed_multiple <- cntrl_bed[(cntrl_bed$end - cntrl_bed$start) > 1, ]
    #how many new rows to make for cntrl?
    cntrl_bed_multiple_calc <- sum(cntrl_bed_multiple$end - cntrl_bed_multiple$start)
    cntrl.cov <- numeric(length = cntrl_bed_multiple_calc)
    cntrl.start <- numeric(length = cntrl_bed_multiple_calc)
    
    #For each row, split into single bp regions with the same coverage; append to dyz1_bed_single
    i <- 1
    for(myRow in c(1:dim(dyz1_bed_multiple)[1])){
      multipleEntry <- dyz1_bed_multiple[myRow,]
      useRegions <- c(multipleEntry$start: (multipleEntry$end - 1))
      num_regions <- length(useRegions)
      dyz1.cov[i:(i + num_regions - 1)] <- multipleEntry$coverage
      dyz1.start[i:(i + num_regions - 1)] <- useRegions
      i <- i + num_regions
    }
    
    cat("Done with splitting DYZ1 coverage\n")
    
    j <- 1
    for(myRow in c(1:dim(cntrl_bed_multiple)[1])){
      multipleEntry <- cntrl_bed_multiple[myRow,]
      useRegions <- c(multipleEntry$start: (multipleEntry$end - 1))
      num_regions <- length(useRegions)
      cntrl.cov[j:(j + num_regions - 1)] <- multipleEntry$coverage
      cntrl.start[j:(j + num_regions - 1)] <- useRegions
      j <- j + num_regions
    }
    
    cat("Done with splitting cntrl region coverage \n\n")
    
    #Generate the new single bed file
    dyz1_split <- data.frame( "chrom" = "chrY", "start" = dyz1.start, "end" = dyz1.start + 1, "coverage" = dyz1.cov, stringsAsFactors = FALSE)
    dyz1_revised_bed <- rbind(dyz1_split,dyz1_bed_single)
    
    cntrl_split <- data.frame( "chrom" = "chrY", "start" = cntrl.start, "end" = cntrl.start + 1, "coverage" = cntrl.cov, stringsAsFactors = FALSE)
    cntrl_revised_bed <- rbind(cntrl_split,cntrl_bed_single)
    
    results_list_control <- c(results_list_control, list(cntrl_revised_bed))
    results_list_dyz1<- c(results_list_dyz1, list(dyz1_revised_bed))
  }
  
  
  # Combine the coverage from the libraries if multiple samples!
  num_samples <- length(results_list_control)
  coverage_ratio <- NULL
  
  results_table_unique[mySubject,"Num_samples"] <- num_samples
  
  
  if(num_samples == 1){
    
    # Calculate coverage ratio
    coverage_ratio <-  mean(results_list_dyz1[[1]]$coverage) / mean(results_list_control[[1]]$coverage)
    
    # Add results to file
    results_table_unique[mySubject,"Single_copy_depth"] <- mean(results_list_control[[1]]$coverage)
    results_table_unique[mySubject,"DYZ1_depth"] <- mean(results_list_dyz1[[1]]$coverage)
    results_table_unique[mySubject,"DYZ1_norm_coverage"] <- mean(results_list_dyz1[[1]]$coverage) / mean(results_list_control[[1]]$coverage)
    
  } else if(num_samples == 2){
    merged_control <- merge(x=results_list_control[[1]], by = "start", y=results_list_control[[2]], all=TRUE)
    merged_control[is.na(merged_control)] <- 0
    merged_control$sum_coverage <- merged_control$coverage.x + merged_control$coverage.y
    
    merged_dyz1 <- merge(x=results_list_dyz1[[1]], by = "start", y=results_list_dyz1[[2]], all=TRUE)
    merged_dyz1[is.na(merged_dyz1)] <- 0
    merged_dyz1$sum_coverage <- merged_dyz1$coverage.x + merged_dyz1$coverage.y
    
    # Add results to file
    results_table_unique[mySubject,"Single_copy_depth"] <- mean(merged_control$sum_coverage)
    results_table_unique[mySubject,"DYZ1_depth"] <- mean(merged_dyz1$sum_coverage)
    results_table_unique[mySubject,"DYZ1_norm_coverage"] <- mean(merged_dyz1$sum_coverage) / mean(merged_control$sum_coverage)
    
  } else if(num_samples == 3){
    merged_control <- merge(x=results_list_control[[1]], by = "start", y=results_list_control[[2]], all=TRUE)
    merged_control_2 <- merge(x=merged_control, by = "start", y=results_list_control[[3]], all=TRUE)
    merged_control_2[is.na(merged_control_2)] <- 0
    merged_control_2$sum_coverage <- merged_control_2$coverage.x + merged_control_2$coverage.y
    
    merged_dyz1 <- merge(x=results_list_dyz1[[1]], by = "start", y=results_list_dyz1[[2]], all=TRUE)
    merged_dyz1_2 <- merge(x=merged_dyz1, by = "start", y=results_list_dyz1[[3]], all=TRUE)
    merged_dyz1_2[is.na(merged_dyz1_2)] <- 0
    merged_dyz1_2$sum_coverage <- merged_dyz1_2$coverage.x + merged_dyz1_2$coverage.y
    
    # Add results to file
    results_table_unique[mySubject,"Single_copy_depth"] <- mean(merged_control_2$sum_coverage)
    results_table_unique[mySubject,"DYZ1_depth"] <- mean(merged_dyz1_2$sum_coverage)
    results_table_unique[mySubject,"DYZ1_norm_coverage"] <- mean(merged_dyz1_2$sum_coverage) / mean(merged_control_2$sum_coverage)
    
  } else if(num_samples == 4){
    merged_control <- merge(x=results_list_control[[1]], by = "start", y=results_list_control[[2]], all=TRUE)
    merged_control_2 <- merge(x=merged_control, by = "start", y=results_list_control[[3]], all=TRUE)
    merged_control_3 <- merge(x=merged_control_2, by = "start", y=results_list_control[[4]], all=TRUE)
    merged_control_3[is.na(merged_control_3)] <- 0
    merged_control_3$sum_coverage <- merged_control_3$coverage.x + merged_control_3$coverage.y
    
    merged_dyz1 <- merge(x=results_list_dyz1[[1]], by = "start", y=results_list_dyz1[[2]], all=TRUE)
    merged_dyz1_2 <- merge(x=merged_dyz1, by = "start", y=results_list_dyz1[[3]], all=TRUE)
    merged_dyz1_3 <- merge(x=merged_dyz1_2, by = "start", y=results_list_dyz1[[4]], all=TRUE)
    merged_dyz1_3[is.na(merged_dyz1_3)] <- 0
    merged_dyz1_3$sum_coverage <- merged_dyz1_3$coverage.x + merged_dyz1_3$coverage.y
    
    # Add results to file
    results_table_unique[mySubject,"Single_copy_depth"] <- mean(merged_control_3$sum_coverage)
    results_table_unique[mySubject,"DYZ1_depth"] <- mean(merged_dyz1_3$sum_coverage)
    results_table_unique[mySubject,"DYZ1_norm_coverage"] <- mean(merged_dyz1_3$sum_coverage) / mean(merged_control_3$sum_coverage)
  }
}

write.table(x = results_table_unique, file="1225_subjects_DYZ1_depth_summed.txt", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
