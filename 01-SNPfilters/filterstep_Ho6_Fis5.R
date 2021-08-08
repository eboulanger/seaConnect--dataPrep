#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# load libraries
library(adegenet)
library(polysat)
library(pegas)
library(vcfR)
library(hierfstat)

start <- Sys.time()
print(paste("run started at",start))

# Read vcf
dat <- read.vcfR(args[1], verbose=F)

# Convert to genind format
genind <- vcfR2genind(dat)

# Assign a sampling site to each indiv
ind <- rownames(genind@tab)
pop <-sub("(.*?)_.*", "\\1", ind)
pop(genind) <- pop

#### Filter SNPs

hf <- genind2hierfstat(genind)
stats <- basic.stats(hf)

# Filter SNPs for Ho and Fis
loc_H.6 <- stats$perloc[which(stats$perloc$Ho<0.6),]
loc_H.6_Fis.5 <- loc_H.6[which(loc_H.6$Fis<0.5 & loc_H.6$Fis>-0.5),]
locnames <-rownames(loc_H.6_Fis.5)
print(paste("# SNPs to keep =", length(locnames)))

# extract loci to remove from the dataset
toremove <-which(locNames(genind) %in% locnames ==F)
print(paste("# SNPs to remove =", length(toremove)))
genind_filter<-genind[loc=-toremove]
posi_id <- levels(genind_filter@loc.fac)

# get positions of loci to keep
# needed to subset vcf file later
positions <- data.frame(row.names = 1:length(posi_id))
for (i in 1:length(posi_id)) {
  positions[i,1] <- paste(strsplit(posi_id, "_")[[i]][1:2], collapse = "_")
}
for (i in 1:length(posi_id)) {
  positions[i,2] <- strsplit(posi_id, "_")[[i]][3]
}

# output positions table
write.table(positions, file = args[2], sep = "\t", quote = F, row.names = F, col.names = F)

# duration of script run
end <- Sys.time()
print(paste("run ended at", end))
end - start