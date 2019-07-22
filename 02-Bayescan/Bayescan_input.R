#----- BayeScan -----#
## produce Bayescan input file

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
outputBayeScan <- paste0(speciesCode,"_bayesc_input_3pop")
  
# import vcf file to genind
vcfFile <- read.vcfR(inputVcf, verbose = F)
!is.null(vcfFile)
genind <- vcfR2genind(vcfFile)
!is.null(genind)

# use DAPC groups as pop
grp <- find.clusters(genind, max.n.clust = 40)
400
2
grp$size 
pop(genind) <- grp$grp

# convert to bayescan format
genomic_converter(genind, output="bayescan", filename = outputBayeScan)

end <- Sys.time()
print(paste("run ended at", end))
end - start
