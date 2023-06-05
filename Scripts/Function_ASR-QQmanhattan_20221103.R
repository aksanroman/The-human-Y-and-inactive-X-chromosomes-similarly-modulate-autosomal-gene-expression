#Creates a manhattan plot
ASRQQmanhattan <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", 
                      col=c("gray10", "gray60"), colhighlight=c("red","black"), chrlabs=NULL,
                      guideline=FALSE, highlight=NULL, logp=TRUE, annotatePval = NULL, annotateTop = TRUE, myLim=NULL, myOneChr = "_", ...) {
  
  # Not sure why, but package check will warn without this.
  CHR=BP=P=index=NULL
  
  # Check for sensible dataset
  ## Make sure you have chr, bp and p columns.
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
  ## warn if you don't have a snp column
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  ## make sure chr, bp, and p columns are numeric.
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  
  # Create a new data.frame with columns called CHR, BP, and P.
  d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]])
  
  # If the input data frame has a SNP column, add it to the new data frame you're creating.
  if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
  
  # Set positions, ticks, and labels for plotting
  ## Sort and keep only values where is numeric.
  #d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  #d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
  if (logp) {
    d$logp <- -log10(d$P)
  } else {
    d$logp <- d$P
  }
  d$pos=NA
  
  
  # Fixes the bug where one chromosome is missing by adding a sequential index column.
  d$index=NA
  ind = 0
  for (i in unique(d$CHR)){
    ind = ind + 1
    d[d$CHR==i,]$index = ind
  }
  
  # This section sets up positions and ticks. Ticks should be placed between chromosomes. The a new pos column is added that keeps a running
  # sum of the positions of each successive chromsome. For example:
  # chr bp pos
  # 1   1  1
  # 1   2  2
  # 2   1  3
  # 2   2  4
  # 3   1  5
  nchr = length(unique(d$CHR))
  if (nchr==1) { ## For a single chromosome
    ## Uncomment the next two lines to plot single chr results in Mb
    options(scipen=999)
    d$pos=d$BP/1e6
    # d$pos=d$BP
    ticks=floor(length(d$pos))+1
    xlabel = paste0("Chromosome ", myOneChr, " position (Mb)")
    labs = ticks
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$BP
      } else {
        lastbase=lastbase+25000000+tail(subset(d,index==i-1)$BP, 1)
        d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
      }
      # Old way: assumes SNPs evenly distributed
      # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
      # New way: doesn't make that assumption
      ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
    }
    xlabel = 'Chromosome'
    #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
    labs <- unique(d$CHR)
  }
  
  # Initialize plot
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  
  # The old way to initialize the plot
  # plot(NULL, xaxt='n', bty='n', xaxs='i', yaxs='i', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
  #      xlab=xlabel, ylab=expression(-log[10](italic(p))), las=1, pch=19, ...)
  
  
  # The new way to initialize the plot.
  ## See http://stackoverflow.com/q/23922130/654296
  ## First, define your default arguments
  if(is.null(myLim)){
  def_args <- list(xaxt='n', bty='n', xaxs='i', yaxs='i', las=2, pch=19,
                   xlim=c(xmin,xmax), ylim=c(floor(min(d$logp)),ceiling(max(d$logp))),
                   xlab=xlabel, ylab=expression(-log[10](italic(p))))
  }
  else{
    def_args <- list(xaxt='n', bty='n', xaxs='i', yaxs='i', las=2, pch=19,
                     xlim=c(xmin,xmax), ylim=myLim,
                     xlab=xlabel, ylab=expression(-log[10](italic(p))))
  }
  ## Next, get a list of ... arguments
  #dotargs <- as.list(match.call())[-1L]
  dotargs <- list(...)
  ## And call the plot function passing NA, your ... arguments, and the default
  ## arguments that were not defined in the ... arguments.
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  
  # If manually specifying chromosome labels, ensure a character vector and number of labels matches number chrs.
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs)==length(labs)) {
        labs <- chrlabs
      } else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    } else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  
  # Add an axis. 
  if (nchr==1) { #If single chromosome, ticks and labels automatic.
    axis(1, ...)
  } else { # if multiple chrs, use the ticks and labels you created above.
    axis(1, at=ticks, labels=labs,las=2,...)
  }
  
  # Create a vector of alternatiting colors
  col=rep(col, max(d$CHR))
  colhighlight=rep(colhighlight, max(d$CHR))
  
  # Add points to plot depending on highlight status
  if (!is.null(highlight)) { #if we do want to highlight...
    if (any(!(highlight %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight=d[which(d$SNP %in% highlight), ] #set the highlight points
    if (nchr==1) {
      #plot normal points
      with(d, points(pos, logp, pch=19, col=col[1],...))
      with(d.highlight, points(pos, logp, pch=19, col=colhighlight[1], ...))
    } else {
      # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
      icol=1
      for (i in unique(d$index)) {
        #plot normal points
        with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=19, ...)) 
        #plot highlight points
        with(d.highlight[d.highlight$index==unique(d$index)[i], ], points(pos, logp, col=colhighlight[icol], pch=19, ...))
        icol=icol+1
      }
    }
  } else {
    if (nchr==1) {
      with(d, points(pos, logp, pch=19, col=col[1], ...)) 
    } else {
      # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
      icol=1
      for (i in unique(d$index)) {
        with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=19, ...))
        icol=icol+1
      }
    }
  }
  
  # Add lines
  if (guideline==TRUE){
    abline(h=1, col="gray19", lty=2)
       abline(h=-1, col="gray19", lty=2)
  } 
 
  
  
  
  # Highlight top SNPs
  if (!is.null(annotatePval)) {
    # extract top SNPs at given p-val
    topHits = subset(d, P <= annotatePval)
    par(xpd = TRUE)
    # annotate these SNPs
    if (annotateTop == FALSE) {
      with(subset(d, P <= annotatePval), 
           textxy(pos, -log10(P), offset = 0.625, labs = topHits$SNP, cex = 0.45), ...)
    }
    else {
      # could try alternative, annotate top SNP of each sig chr
      topHits <- topHits[order(topHits$P),]
      topSNPs <- NULL
      
      for (i in unique(topHits$CHR)) {
        
        chrSNPs <- topHits[topHits$CHR == i,]
        topSNPs <- rbind(topSNPs, chrSNPs[1,])
        
      }
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }  
  par(xpd = FALSE)
}

ASRQQmanhattan_blackBg <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", 
                           col=c("gray10", "gray60"), colhighlight=c("red","black"), chrlabs=NULL,
                           guideline=FALSE, highlight=NULL, logp=TRUE, annotatePval = NULL, annotateTop = TRUE, myLim=NULL, myOneChr = "_", ...) {
  
  # Not sure why, but package check will warn without this.
  CHR=BP=P=index=NULL
  
  # Check for sensible dataset
  ## Make sure you have chr, bp and p columns.
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
  ## warn if you don't have a snp column
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  ## make sure chr, bp, and p columns are numeric.
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  
  # Create a new data.frame with columns called CHR, BP, and P.
  d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]])
  
  # If the input data frame has a SNP column, add it to the new data frame you're creating.
  if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
  
  # Set positions, ticks, and labels for plotting
  ## Sort and keep only values where is numeric.
  #d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  #d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
  if (logp) {
    d$logp <- -log10(d$P)
  } else {
    d$logp <- d$P
  }
  d$pos=NA
  
  
  # Fixes the bug where one chromosome is missing by adding a sequential index column.
  d$index=NA
  ind = 0
  for (i in unique(d$CHR)){
    ind = ind + 1
    d[d$CHR==i,]$index = ind
  }
  
  # This section sets up positions and ticks. Ticks should be placed between chromosomes. The a new pos column is added that keeps a running
  # sum of the positions of each successive chromsome. For example:
  # chr bp pos
  # 1   1  1
  # 1   2  2
  # 2   1  3
  # 2   2  4
  # 3   1  5
  nchr = length(unique(d$CHR))
  if (nchr==1) { ## For a single chromosome
    ## Uncomment the next two lines to plot single chr results in Mb
    options(scipen=999)
    d$pos=d$BP/1e6
    # d$pos=d$BP
    ticks=floor(length(d$pos))+1
    xlabel = paste0("Chromosome ", myOneChr, " position (Mb)")
    labs = ticks
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$BP
      } else {
        lastbase=lastbase+25000000+tail(subset(d,index==i-1)$BP, 1)
        d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
      }
      # Old way: assumes SNPs evenly distributed
      # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
      # New way: doesn't make that assumption
      ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
    }
    xlabel = 'Chromosome'
    #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
    labs <- unique(d$CHR)
  }
  
  # Initialize plot
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  
  # The old way to initialize the plot
  # plot(NULL, xaxt='n', bty='n', xaxs='i', yaxs='i', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
  #      xlab=xlabel, ylab=expression(-log[10](italic(p))), las=1, pch=19, ...)
  
  
  # The new way to initialize the plot.
  ## See http://stackoverflow.com/q/23922130/654296
  ## First, define your default arguments
  if(is.null(myLim)){
    def_args <- list(xaxt='n', yaxt='n', bty='n', xaxs='i', yaxs='i', las=2, pch=19,
                     xlim=c(xmin,xmax), ylim=c(floor(min(d$logp)),ceiling(max(d$logp))),
                     xlab=xlabel, ylab=expression(-log[10](italic(p))), col.lab = "white", col.axis = "white")
  }
  else{
    def_args <- list(xaxt='n',yaxt='n', bty='n', xaxs='i', yaxs='i', las=2, pch=19,
                     xlim=c(xmin,xmax), ylim=myLim,
                     xlab=xlabel, ylab=expression(-log[10](italic(p))), col.lab = "white", col.axis = "white")
  }
  ## Next, get a list of ... arguments
  #dotargs <- as.list(match.call())[-1L]
  dotargs <- list(...)
  ## And call the plot function passing NA, your ... arguments, and the default
  ## arguments that were not defined in the ... arguments.
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  
  # If manually specifying chromosome labels, ensure a character vector and number of labels matches number chrs.
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs)==length(labs)) {
        labs <- chrlabs
      } else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    } else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  
  # Add an axis. 
  if (nchr==1) { #If single chromosome, ticks and labels automatic.
    axis(1, col="white",...)
    axis(2, col="white",...)
  } else { # if multiple chrs, use the ticks and labels you created above.
    axis(1, at=ticks, labels=labs,las=2,col="white",...)
    axis(2, col="white",...)
  } 
  
  # Create a vector of alternatiting colors
  col=rep(col, max(d$CHR))
  colhighlight=rep(colhighlight, max(d$CHR))
  
  # Add points to plot depending on highlight status
  if (!is.null(highlight)) { #if we do want to highlight...
    if (any(!(highlight %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight=d[which(d$SNP %in% highlight), ] #set the highlight points
    if (nchr==1) {
      #plot normal points
      with(d, points(pos, logp, pch=19, col=col[1],...))
      with(d.highlight, points(pos, logp, pch=19, col=colhighlight[1], ...))
    } else {
      # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
      icol=1
      for (i in unique(d$index)) {
        #plot normal points
        with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=19, ...)) 
        #plot highlight points
        with(d.highlight[d.highlight$index==unique(d$index)[i], ], points(pos, logp, col=colhighlight[icol], pch=19, ...))
        icol=icol+1
      }
    }
  } else {
    if (nchr==1) {
      with(d, points(pos, logp, pch=19, col=col[1], ...)) 
    } else {
      # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
      icol=1
      for (i in unique(d$index)) {
        with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=19, ...))
        icol=icol+1
      }
    }
  }
  
  # Add lines
  if (guideline==TRUE){
    abline(h=1, col="gray19", lty=2)
    abline(h=-1, col="gray19", lty=2)
  } 
  
  
  
  
  # Highlight top SNPs
  if (!is.null(annotatePval)) {
    # extract top SNPs at given p-val
    topHits = subset(d, P <= annotatePval)
    par(xpd = TRUE)
    # annotate these SNPs
    if (annotateTop == FALSE) {
      with(subset(d, P <= annotatePval), 
           textxy(pos, -log10(P), offset = 0.625, labs = topHits$SNP, cex = 0.45), ...)
    }
    else {
      # could try alternative, annotate top SNP of each sig chr
      topHits <- topHits[order(topHits$P),]
      topSNPs <- NULL
      
      for (i in unique(topHits$CHR)) {
        
        chrSNPs <- topHits[topHits$CHR == i,]
        topSNPs <- rbind(topSNPs, chrSNPs[1,])
        
      }
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }  
  par(xpd = FALSE)
} 
 