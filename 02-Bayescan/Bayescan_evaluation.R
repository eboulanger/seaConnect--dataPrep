setwd("~/Documents/project_SEACONNECT/seaConnect--dataPrep/02-Bayescan/")
#----- BayeScan -----#
## evaluate model convergence and extract results
## script adapted from http://evomicsorg.wpengine.netdna-cdn.com/wp-content/uploads/2016/01/BayeScan_BayeScEnv_exercises.pdf

# load libraries
library(coda)
library(vcfR)

# set arguments
spCode <- "mul" # species
runNum <- 1 # run number
inputName <- paste0(spCode,"_bayesc_2pop_run",runNum) # 2pop or 3pop analysis (only for mullus)
#vcfPath <- "../01-SNPfilters/02-Mullus/mul_all_filtered.vcf"
#vcfPath <- "../01-SNPfilters/01-Diplodus/dip_all_filtered.vcf"

## model convergence ####
chain1 <- read.table(paste0(inputName,"_ch1.sel"),header=TRUE)
#chain<-chain[-c(1)] # remove first column (step number)
# create a MCMC object with the correct thinning interval
chain1 <- mcmc(chain1,thin=10)
# Plot
plot(chain1)
# Summary statistics
summary(chain1)
# Test correlation between sampled parameters to verify that the sample is large enough
autocorr.diag(chain1)
# Effective size of the sample
effectiveSize(chain1)

# Statistic test for convergence
#### Geweke's convergence diagnistic
#The diagnostic reports the z- scores for each parameter. 
#With ?? = 0.05, We reject H0 (equality of means => convergence) if z < -1.96 or z > +1.96.
geweke <- geweke.diag(chain1, frac1=0.1, frac2=0.5)
geweke
which(abs(geweke$z)>1.96) # 0 parameter does not converge

#### Heidelberg and Welch's test
heidel.diag(chain1, eps=0.1, pvalue=0.05) # all parameters passed

#### Test comparing two chains
chain2 <- read.table(paste0(inputName,"_ch2.sel"),header=TRUE)
#chain2<-chain2[-c(1)]
chain2 <- mcmc(chain2,thin=10)
plot(chain2)

combined = mcmc.list(chain1,chain2) 
plot(combined)
gelman.diag(combined) # should be close to 1
gelman.plot(combined,ask)

## Bayescan results ####

## Functions defined in the script 'plot_R.r'in the BayeScan folder
results<-plot_bayescan(paste0(inputName,"_ch2_fst.txt"),FDR=0.05)

# detect outliers
outliers_bayescan <-results$outliers # position of outliers
results$nb_outliers #283 outliers for mullus, 126 for diplodus

# extract corresponding SNP names and positions
vcf <- read.vcfR(vcfPath, verbose=F)
loci <- as.data.frame(vcf@fix[,1:2])
outlier_loci <- loci[results$outliers,]

# output positions table
write.table(outlier_loci, file = paste0("outl_pos_bayesc_",spCode,"_run",runNum,".txt"), sep = "\t", quote = F, row.names = F, col.names = F)
print(paste0("outlier positions table exported to outl_pos_bayesc_", spCode,"_run", runNum,".txt"), quote = 0)


