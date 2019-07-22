#----- PCAdapt -----#
# vignette: https://cran.r-project.org/web/packages/pcadapt/vignettes/pcadapt.html

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#args <- c("../01-SNPfilters/01-Diplodus/dip_all_filtered.vcf", "dip")
#args <- c("../01-SNPfilters/02-Mullus/mul_all_filtered.vcf", "mul")

# load libraries
library(pcadapt)
library(vcfR)

# define arguments
path_to_file <- args[1]
speciesCode <- args[2]

# import vcf file to pcadapt format
dat <- read.pcadapt(path_to_file, type = "vcf")

# choose number K of principal components
x <- pcadapt(input = dat, K = 20)
# set plot output file name
pdf(file=paste0("RPlots_", args[2],".pdf"))
# scree plot of principal components
plot(x, option = "screeplot", K = 10)
# score plot of principal components
plot(x, option = "scores")
plot(x, option = "scores", i = 3, j = 4) # check PC axis 3 and 4

# compute outliers detection with K = 4
x <- pcadapt(dat, K = 4)
# Manhattan plot
#jpeg(filename = paste0(args[2],"_manhattan_pcadpt",".jpeg"))
plot(x, option = "manhattan")
#dev.off()
# QQ plot
plot(x , option = "qqplot", threshold = 0.05)
# histogram of p-values
#jpeg(filename = paste0(args[2],"_histP_pcadpt",".jpeg"))
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
#dev.off()
# histogram of test statistic
plot(x, option = "stat.distribution")

# choose cutoff for outlier detection
# Bonferroni correction
padj <- p.adjust(x$pvalues, method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
print(paste(length(outliers), "outliers detected"),quote=0)

# extract corresponding SNP names and positions
vcf <- read.vcfR(path_to_file, verbose=F)
loci <- as.data.frame(vcf@fix[,1:2])
outlier_loci <- loci[outliers,]

# output positions table
write.table(outlier_loci, file = paste0("outl_pos_pcadpt_",speciesCode,".txt"), sep = "\t", quote = F, row.names = F, col.names = F)
print(paste0("outlier positions table exported to outl_pos_pcadpt_",speciesCode), quote = 0)
