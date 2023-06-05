library(RBGL)
library(Vennerable)

threeGroupVenn_geneNames <- function(setA, setAname, setB, setBname,setC, setCname, expGenes, filename){
  #setA, B, C should be names of genes
  #expGenes should be vector of gene names
  mySet <- list(setA, setB, setC)
  names(mySet) <- c(setAname,setBname,setCname)
  myVenn <- Venn(mySet)
  pdf(myVenn, file=filename)
  plot(myVenn, show = list(SetLabels = TRUE,Faces = FALSE, FaceText = NULL))
  plot(myVenn, show = list(SetLabels = TRUE,Faces = FALSE))
  plot(myVenn, show = list(SetLabels = FALSE, Faces=FALSE, FaceText=NULL))
  dev.off() 
  
  #calculating hypergeometric test for each overlap
  expGenes <- length(expGenes)
  
  print(paste0("Overlap of ",setAname," and ",setBname))
  overlapAB <- length(intersect(setA, setB))
  print(paste0("overlap = ", as.character(overlapAB)))
  m_A <- length(setA)
  print(paste0("marked = ", as.character(m_A)))
  print(paste0("expressed = ", as.character(expGenes)))
  nonMarked_A <- (expGenes - m_A)
  print(paste0("non marked = ", as.character(nonMarked_A)))
  k_B <- length(setB)
  print(paste0("number chosen in set B = ", as.character(k_B)))
  print(phyper(overlapAB -1,m_A, nonMarked_A, k_B, lower.tail=FALSE))
  
  print(paste0("Overlap of ",setBname," and ",setCname))
  overlapBC <- length(intersect(setB, setC))
  print(paste0("overlap = ", as.character(overlapBC)))
  m_B <- length(setB)
  print(paste0("marked = ", as.character(m_B)))
  print(paste0("expressed = ", as.character(expGenes)))
  nonMarked_B <- (expGenes - m_B)
  print(paste0("non marked = ", as.character(nonMarked_B)))
  k_C <- length(setC)
  print(paste0("number chosen in set C = ", as.character(k_C)))
  print(phyper(overlapBC -1,m_B, nonMarked_B, k_C, lower.tail=FALSE))
  
  print(paste0("Overlap of ",setAname," and ",setCname))
  overlapAC <- length(intersect(setA, setC))
  print(paste0("overlap = ", as.character(overlapAC)))
  print(paste0("marked = ", as.character(m_A)))
  print(paste0("expressed = ", as.character(expGenes)))
  print(paste0("non marked = ", as.character(nonMarked_A)))
  print(paste0("number chosen in set c = ", as.character(k_C)))
  print(phyper(overlapAC -1,m_A, nonMarked_A, k_C, lower.tail=FALSE))
}