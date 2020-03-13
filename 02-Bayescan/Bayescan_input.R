#----- BayeScan -----#
## produce Bayescan input file
# setwd("~/Documents/project_SEACONNECT/seaConnect--dataPrep/02-Bayescan/")

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#args <- c("../01-SNPfilters/01-Diplodus/dip_all_filtered.vcf","dip")
#args <- c("../01-SNPfilters/02-Mullus/mul_all_filtered.vcf","mul")

start <- Sys.time()

# load libraries
library(vcfR)
library(adegenet)
library(radiator)

# define arguments
inputVcf <- args[1]
speciesCode <- args[2]
outputBayeScan <- paste0(speciesCode,"_bayesc_input_2pop")
  
# import vcf file to genind
vcfFile <- read.vcfR(inputVcf, verbose = F)
!is.null(vcfFile)
genind <- vcfR2genind(vcfFile)
!is.null(genind)

# use DAPC groups as pop
grp <- find.clusters(genind, max.n.clust = 40, n.pca = round(nrow(genind$tab)/3))
2
grp$size 
pop(genind) <- grp$grp
pop(genind) <- paste0("pop", pop(genind))

# convert to bayescan format
tidy_genind <- tidy_genomic_data(genind)
write_bayescan(tidy_genind, pop.select = levels(pop(genind)), filename = outputBayeScan)

end <- Sys.time()
print(paste("run ended at", end))
end - start
