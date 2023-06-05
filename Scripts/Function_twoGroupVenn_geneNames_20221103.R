library("RBGL")
library("Vennerable")
twoGroupVenn_geneNames <- function(setA, setAname, setB, setBname, expGenes, filename){
  #setA = character vector of gene names
  #setAname = label for the venn diagram
  #setB = character vector of gene names
  #setBname = label for the venn diagram
  #expGenes = character vector of expressed genes
  #filename = output filename.pdf for the venn diagram picture
  
  mySet <- list(setA, setB)
  names(mySet) <- c(setAname,setBname)
  myVenn <- Venn(mySet)
  myOverlap <- length(intersect(setA, setB))
  print(paste0("overlap = ", as.character(myOverlap)))
  m <- length(setA)
  print(paste0("marked = ", as.character(m)))
  expGenes <- length(expGenes)
  print(paste0("expressed = ", as.character(expGenes)))
  nonMarked <- (expGenes - m)
  print(paste0("non marked = ", as.character(nonMarked)))
  k <- length(setB)
  print(paste0("number chosen in set b = ", as.character(k)))
  print(phyper(myOverlap -1,m, nonMarked, k, lower.tail=FALSE))
  
  pdf(myVenn, file=filename)
  #plot(myVenn, show = list(SetLabels = TRUE,Faces = FALSE, FaceText = NULL))
  plot(myVenn, show = list(SetLabels = TRUE,Faces = FALSE))
  #plot(myVenn, show = list(SetLabels = FALSE, Faces=FALSE, FaceText=NULL))
  dev.off()
}