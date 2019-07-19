# SNP filtering

script adapted from [ddocent tutorial](tutorial: https://www.ddocent.com/filtering/)

## run the script
move to the correct directory and make the output directories
```
cd ~/Documents/project_SEACONNECT/seaConnect--dataPrep/01-SNPfilters/
mkdir 01-Diplodus
mkdir 02-Mullus
```
set the species-specific arguments for the script to run on
$1 = input file (specify path if not in current directory)
$2 = output folder
$3 = species code
```
bash filtering.sh ../00-rawData/01-Diplodus/sar_ddocent.vcf 01-Diplodus dip
bash filtering.sh ../00-rawData/02-Mullus/mullus.vcf 02-Mullus mul
```

## SNP filtering results for Mullus surmuletus					
### ddocent pipeline + additions					
				
| filtering step | filter for                                        | individuals retained | SNPs retained | run time (sec) | output |
| -------------- | ------------------------------------------------- | :------------------: | :-----------: | :------------: | ------ |
| step 0         | ddocent output data                               | 431                  | 49027         |                | mullus.vcf
| step 1         | call below 50%, mac < 3, min quality score = 30   | 431                  | 49027         | 71.00          | mullus.g5mac3.recode.vcf
| step 2         | min mean depth genotype call = 3                  | 431                  | 49027         | 74.00          | mullus.g5mac3dp3.recode.vcf
| step 3         | individuals > 50% missing data                    | 424                  | 49027         | 70.00          | mullus.g5mac3dplm.recode.vcf
| step 4         | sites > 5% missing data, maf 0.05, min meanDP = 5 | 424                  | 18727         | 14.00          | DP3g95maf05.recode.vcf
| step 5         | filter for allele balance                         |                      | 17965         |                | DP3g95maf05.fil1.vcf
| step 6         | filter out sites with reads from both strands     |                      | **SKIP**      |	             | 
| step 7         | ration mapping qualities ref vs alternate alleles |                      | 17546         |                | DP3g95maf05.fil3.vcf
| step 8         | paired status                                     |                      | 17546         |                | DP3g95maf05.fil4.vcf
| step 9         | remove sites quality score < 1/4 depth            |                      | 17546         |                | DP3g95maf05.fil5.vcf
| step 10        | depth x quality score cutoff	                     | 424                  | 15466         |	             | 
| step 11        | He > 0.6 & Fis > 0.5 & Fix < -0.5                 | 424                  | 15232         | 25 min         | DP3g95maf05.FIL.HFis.recode.vcf
| step 12        | rename                                            |                      |               |                | mul_all_filtered.vcf

## SNP filtering results for Diplodus sargus				
### ddocent pipeline + additions					
				
| filtering step | filter for                                        | individuals retained | SNPs retained | run time (sec) | output |
| -------------- | ------------------------------------------------- | :------------------: | :-----------: | :------------: | ------ |
| step 0         | ddocent output data                               | 297                  | 13362         |                | diplodus.vcf
| step 1         | call below 50%, mac < 3, min quality score = 30   | 297                  | 13362         | 14.00          | g5mac3.recode.vcf
| step 2         | min mean depth genotype call = 3                  | 297                  | 13362         | 14.00          | g5mac3dp3.recode.vcf
| step 3         | individuals > 50% missing data                    | 297                  | 13362         | 16.00          | g5mac3dplm.recode.vcf
| step 4         | sites > 5% missing data, maf 0.05, min meanDP = 5 | 297                  | 10389         | 11.00          | DP3g95maf05.recode.vcf
| step 5         | filter for allele balance                         | 297                  | 10202         |                | DP3g95maf05.fil1.vcf
| step 6         | filter out sites with reads from both strands     | SKIP                 | **SKIP**      |                | SKIP
| step 7         | ration mapping qualities ref vs alternate alleles | 297                  | 9689          |                | DP3g95maf05.fil3.vcf
| step 8         | paired status                                     | 297                  | 9689          |                | DP3g95maf05.fil4.vcf
| step 9         | remove sites quality score < 1/4 depth            | 297                  | 9688          |                | DP3g95maf05.fil5.vcf
| step 10        | depth x quality score cutoff	                     | 297                  | 8325          | 11.00          | 
| step 11        | He > 0.6 & Fis > 0.5 & Fix < -0.5                 | 297                  | 8206          | 27 min         | DP3g95maf05.FIL.HFis.recode.vcf
| step 12        | rename                                            |                      |               |                | dip_all_filtered.vcf 
