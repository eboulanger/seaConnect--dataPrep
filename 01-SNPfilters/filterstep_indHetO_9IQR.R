#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(gridExtra)

start <- Sys.time()
print(paste("run started at",start))

# Read homozygosity file
dat <- read.csv(args[1], sep = "\t")
# dat <- read.csv("~/Documents/project_SEACONNECT/seaConnect--dataPrep/01-SNPfilters/02-Mullus/DP3g95maf05.FIL.HFis.indHo.csv", sep = "\t")
# dat <- read.csv("~/Documents/project_SEACONNECT/seaConnect--dataPrep/01-SNPfilters/01-Diplodus/DP3g95maf05.FIL.HFis.indHo.csv", sep = "\t")
# calculate heterozygosity
dat$O.HET. <- dat$N_SITES - dat$O.HOM.

# in interactive mode: visualise heterozygosity to determine 
# - if individuals need to be removed
# - adapted boxplot range
pdf(file="plot_indv_HETo.pdf")
   plot(dat$O.HET., main = "observed heterozygosity per individual")
   boxplot(dat$O.HET., main = "observed heterozygosity per individual")
   boxplot(dat$O.HET., range = 6, main = "o het with interquartal range = 6")

# extract extreme individuals
#OutVals <- boxplot(dat$O.HET., range = 9)$out
OutVals <- boxplot(dat$O.HET., range = 6)$out
plot.new()
grid.table(dat[dat$O.HET. %in% OutVals,])
dev.off()

OutVals.ind <- as.character(dat$INDV[which(dat$O.HET. %in% OutVals)])
print(paste("# Individuals to remove =", length(OutVals.ind)))
print(cbind(OutVals.ind, OutVals))

indToKeep <- setdiff(as.character(dat$INDV),OutVals.ind)
print(paste("# Individuals to keep =", length(indToKeep)))

# output individuals list
write.table(indToKeep, args[3], sep = "\t", quote = F, row.names = F, col.names = F)

# duration of script run
end <- Sys.time()
print(paste("run ended at", end))
end - start
