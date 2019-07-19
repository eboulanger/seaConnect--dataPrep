# Outlier detection
## with BayeScan

### step 1: convert vcf files to Bayescan .txt files

This script will load your vcf file, determine the population identifier for each individual, and return a bayescan .txt inputfile
set arguments : argument 1 = input file (vcf), argument 2 = species code
script is currently set to detect and assign two populations. Run the script interactively if you want to determine and assign other values of K

#### for mullus :
```
Rscript --vanilla Bayescan_input.R ../01-SNPfilters/02-Mullus/mullus_all_filtered.vcf mul
```

#### for diplodus :
```
Rscript --vanilla Bayescan_input.R ../01-SNPfilters/01-Diplodus/diplodus_all_filtered.vcf dip
```

### step 2: determine outliers with Bayescan

download and compile Bayescan from [here](http://cmpg.unibe.ch/software/BayeScan/download.html) and copy the executable to your local bin folder:
```
cp ~/programmes/BayeScan2.1/binaries/BayeScan2.1_macos64bits /usr/local/bin/
```
#### run the BayeScan model for mullus and diplodus data
multiple runs with different parameters and different input data, with two chains per run to compare convergence later
```
bash run_bayescan.sh
```
### step 3: verify convergence
Run interactive R script called `Bayescan_evaluation.R`
run 1 seems to get "best" results for diplodus. Neither runs detects outliers for mullus.
extract outlier list and export loci positions for later subsetting
