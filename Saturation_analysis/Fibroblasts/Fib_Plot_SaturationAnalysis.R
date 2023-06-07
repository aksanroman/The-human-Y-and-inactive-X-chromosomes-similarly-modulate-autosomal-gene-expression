#Making plots
suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("ggplot2"))


set.seed(seed = 1)
#read in all of the results
results10 <- data.frame(Sample_num = rep(10),read.delim(file="fib_saturationResults_100iterations_10.txt"))
results20 <- data.frame(Sample_num = rep(20),read.delim(file="fib_saturationResults_100iterations_20.txt"))
results30 <- data.frame(Sample_num = rep(30),read.delim(file="fib_saturationResults_100iterations_30.txt"))
results40 <- data.frame(Sample_num = rep(40),read.delim(file="fib_saturationResults_100iterations_40.txt"))
results50 <- data.frame(Sample_num = rep(50),read.delim(file="fib_saturationResults_100iterations_50.txt"))
results60 <- data.frame(Sample_num = rep(60),read.delim(file="fib_saturationResults_100iterations_60.txt"))
results70 <- data.frame(Sample_num = rep(70),read.delim(file="fib_saturationResults_100iterations_70.txt"))
results80 <- data.frame(Sample_num = rep(80),read.delim(file="fib_saturationResults_100iterations_80.txt"))
results90 <- data.frame(Sample_num = rep(90),read.delim(file="fib_saturationResults_100iterations_90.txt"))
results99 <- data.frame(Sample_num = rep(100),read.delim(file="fib_saturationResults_100iterations_99.txt"))

#make the table
saturationResults <- NULL
saturationResults <- rbind(saturationResults,results10, results20, results30, results40, results50, results60, results70, results80, results90,results99)

#Plot the results
data_summary <- function(x){
        # x <- gnomad_pLI_curated_NPX
        m <- median(x, na.rm=TRUE)
        ymin <- quantile(x, na.rm=TRUE)[2]
        ymax <- quantile(x, na.rm=TRUE)[4]
        test <- c("y"=m,"ymin"=ymin,"ymax"=ymax)
        names(test) <- c("y","ymin","ymax")
        return(test)
}

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

pdf(file = "saturation_plots_autosomes_100iterations.pdf", width=2.5, height=2.5)

ggplot(data=saturationResults, mapping = aes(x=as.factor(Sample_num), y=x_auto_sig)) + 
  geom_violin(size=0.25, scale = "width", col="#00000000", fill="00000040") +
  stat_summary(fun.data=data_summary, size=0.25) +
  theme_classic(base_size = 8) + 
  theme(
    axis.text = element_text(color = "black"), 
    axis.ticks = element_line(color = "black", size = 0.25), 
    axis.line = element_line(size = 0.25),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "none",
  ) +
  labs(
    title="X-responsive autosomal genes, LCLs",
    y="Number of genes (FDR<0.05)",
    x="Sample size"
  ) +
  scale_y_continuous(limits = c(0,1000))

ggplot(data=saturationResults, mapping = aes(x=as.factor(Sample_num), y=y_auto_sig)) + 
  geom_violin(size=0.25, scale = "width", col="#00000000", fill="00000040") +
  stat_summary(fun.data=data_summary, size=0.25) +
  theme_classic(base_size = 8) + 
  theme(
    axis.text = element_text(color = "black"), 
    axis.ticks = element_line(color = "black", size = 0.25), 
    axis.line = element_line(size = 0.25),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "none",
  ) +
  labs(
    title="Y-responsive autosomal genes, LCLs",
    y="Number of genes (FDR<0.05)",
    x="Sample size"
  ) +
  scale_y_continuous(limits = c(0,500))


dev.off()
