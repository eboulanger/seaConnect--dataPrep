#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

start <- Sys.time()
print(paste("run started at",start))

# Read homozygosity file
dat <- read.csv(args[1], sep = "\t")
# calculate heterozygosity
dat$O.HET. <- dat$N_SITES - dat$O.HOM.

# in interactive mode: visualise heterozygosity to determine 
# - if individuals need to be removed
# - adapted boxplot range
#   plot(dat$O.HET.)
#   boxplot(dat$O.HET.)
#   boxplot(dat$O.HET., range = 9)

# extract non-extreme individuals
OutVals <- boxplot(dat$O.HET., range = 9)$out

OutVals.ind <- as.character(dat$INDV[which(dat$O.HET. %in% OutVals)])
print(paste("# Individuals to remove =", length(OutVals.ind)))
print(cbind(OutVals.ind, OutVals))

indToKeep <- setdiff(as.character(dat$INDV),OutVals.ind)
print(paste("# Individuals to keep =", length(indToKeep)))

# output individuals list
write.table(indToKeep, args[2], sep = "\t", quote = F, row.names = F, col.names = F)

# duration of script run
end <- Sys.time()
print(paste("run ended at", end))
end - start
