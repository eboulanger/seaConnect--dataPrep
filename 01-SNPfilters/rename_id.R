#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# read id file
id <- read.table(paste0("id_",args,".txt"), stringsAsFactors=F)

# correct names

if(args[1] == "mul") {
  
# mullus
cell <- gsub("_.*","", id$V1)
cell <- paste0("C", cell)
ind <- gsub(".*_","", id$V1)
ind <- gsub("E.*", "", ind)
newid <- paste0(cell,"i", ind) 

} else {
  
# dip
cell <- gsub("_.*","", id$V1)
ind <- gsub(".*i","", id$V1)
ind <- gsub("E.*", "", ind)
newid <- paste0(cell, "i", ind) 
}

# check
bothid <- cbind(id, newid)
print(bothid)

# export
#write.table(newid, file = paste0("new_id_", args,".txt"), quote = F, row.names = F, col.names = F)
write.table(bothid, file = paste0("new_id_", args,".txt"), sep = " ", quote = F, row.names = F, col.names = F)

# create popmap
popmap <- cbind(newid, cell)
colnames(popmap) <- c("INDIVIDUALS", "STRATA")
print(popmap)
nind <- nrow(popmap)
write.table(popmap, file = paste0(args,"_population_map_",nind,"ind.txt"), sep = "\t", quote = F, row.names = F)
