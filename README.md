# seaConnect--dataPrep

Scripts to prepare GBS data for analysis. 
- filtering steps for dDocent SNP output
- outlier detection
- file conversions, subsetting and renaming 

## Dependencies
You will need to install the following software:  
- [VCFtools](https://vcftools.github.io)
- [BCFtools](https://samtools.github.io/bcftools/)
- [PLINK v1.9](https://www.cog-genomics.org/plink/1.9/)
- [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/)

You will need to have the following R packages:  
adegenet, polysat, pegas, vcfR, hierfstat, coda, radiator, pcadapt


## 01-SNPfiltering 

script adapted from [ddocent tutorial](https://www.ddocent.com/filtering/) and additions

### run the script
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
The final step consists of removing individuals with (extreme) outlier heterozygosity values. Here we define outliers as falling outside 6 * interquartile range. 
This theshold value is set in accordance to our data and it is advised to look at the outputted figures to validate this choice for your data.

### SNP filtering results for Mullus surmuletus					
				
| filtering step | filter for                                        | individuals retained | SNPs retained | run time (sec) | output |
| -------------- | ------------------------------------------------- | :------------------: | :-----------: | :------------: | ------ |
| step 0         | ddocent output data                               | 526                  |               |                | mullus.vcf
| step 1         | call below 50%, mac < 3, min quality score = 30   | 526                  | 549881        | 71.00          | mullus.g5mac3.recode.vcf
| step 2         | min mean depth genotype call = 3                  | 482                  | 549881        | 74.00          | mullus.g5mac3dp3.recode.vcf
| step 3         | individuals > 50% missing data                    | 482                  | 549881        | 70.00          | mullus.g5mac3dplm.recode.vcf
| step 4         | sites > 5% missing data, maf 0.05, min meanDP = 5 | 482                  | 3657          | 14.00          | DP3g95maf05.recode.vcf
| step 5         | filter for allele balance                         | 482                  | 3421          |                | DP3g95maf05.fil1.vcf
| step 6         | filter out sites with reads from both strands     | 482                  | **SKIP**      |	             | 
| step 7         | ration mapping qualities ref vs alternate alleles | 482                  | 3284          |                | DP3g95maf05.fil3.vcf
| step 8         | paired status                                     | 482                  | 3284          |                | DP3g95maf05.fil4.vcf
| step 9         | remove sites quality score < 1/4 depth            | 482                  | 3263          |                | DP3g95maf05.fil5.vcf
| step 10        | depth x quality score cutoff	                     | 482                  | 2857          |	             | 
| step 11        | He > 0.6 & Fis > 0.5 & Fix < -0.5                 | 482                  | 2794          | 25 min         | DP3g95maf05.FIL.HFis.recode.vcf
| step 12        | remove extreme outliers individual O HET          | 482                  | 2794          | 23.00          | DP3g95maf05.FIL.HFis.indHet.recode.vcf
| step 13        | rename                                            |                      |               |                | mul_all_filtered_origid.vcf

### SNP filtering results for Diplodus sargus								
				
| filtering step | filter for                                        | individuals retained | SNPs retained | run time (sec) | output |
| -------------- | ------------------------------------------------- | :------------------: | :-----------: | :------------: | ------ |
| step 0         | ddocent output data                               | 297                  | 13362         |                | sar_ddocent_.vcf
| step 1         | call below 50%, mac < 3, min quality score = 30   | 297                  | 13362         | 14.00          | g5mac3.recode.vcf
| step 2         | min mean depth genotype call = 3                  | 297                  | 13362         | 14.00          | g5mac3dp3.recode.vcf
| step 3         | individuals > 50% missing data                    | 297                  | 13362         | 16.00          | g5mac3dplm.recode.vcf
| step 4         | sites > 5% missing data, maf 0.05, min meanDP = 5 | 297                  | 10389         | 11.00          | DP3g95maf05.recode.vcf
| step 5         | filter for allele balance                         | 297                  | 10202         |                | DP3g95maf05.fil1.vcf
| step 6         | filter out sites with reads from both strands     |                      | **SKIP**      |                | SKIP
| step 7         | ration mapping qualities ref vs alternate alleles | 297                  | 9689          |                | DP3g95maf05.fil3.vcf
| step 8         | paired status                                     | 297                  | 9689          |                | DP3g95maf05.fil4.vcf
| step 9         | remove sites quality score < 1/4 depth            | 297                  | 9688          |                | DP3g95maf05.fil5.vcf
| step 10        | depth x quality score cutoff	                     | 297                  | 8325          | 11.00          | 
| step 11        | He > 0.6 & Fis > 0.5 & Fix < -0.5                 | 297                  | 8206          | 27 min         | DP3g95maf05.FIL.HFis.recode.vcf
| step 12        | remove extreme outliers individual O HET          | 297                  | 8206          |                | DP3g95maf05.FIL.HFis.indHet.recode.vcf
| step 13        | rename                                            |                      |               |                | dip_all_filtered_origid.vcf 

### rename individuals for conventional naming system and create population map
this is adapted for our case, review R script or skip according to naming needs
```
bash renaming.sh 01-Diplodus/ dip
bash renaming.sh 02-Mullus/ mul
```
output:
- dip_all_filtered.vcf
- dip_popualtion_map_297ind.txt
- mul_all_filtered.vcf
- mul_population_map_467ind.txt

## 03-PCAdapt

detect outliers with [PCAdapt](https://bcm-uga.github.io/pcadapt/articles/pcadapt.html) in 
R and retrieve loci positions from vcf file

### Run the PCAdapt.R script

This script will load your vcf file, detect outlier loci using the package pcadapt and
output a list of outlier loci positions that can be used to subset your vcf file

EDIT: must run interactively to choose K for each species

set arguments :  
$1 = input file (vcf)  
$2 = species code  

#### for mullus :
```
Rscript --vanilla PCAdapt.R ../01-SNPfilters/02-Mullus/mul_all_filtered.vcf mul
```
see how many outliers were detected :
```
wc -l outl_pos_pcadpt_mul.txt
```
291 outliers detected for mullus

#### for diplodus :
```
Rscript --vanilla PCAdapt.R ../01-SNPfilters/01-Diplodus/dip_all_filtered.vcf dip

wc -l outl_pos_pcadpt_dip.txt
```
413 outliers detected for diplodus

## 04-finalPrep

Final steps to get neutral and adaptive SNP sets and correct file formats for both species.

This script subsets the filtered vcf file from 01-SNPfilters by outlier positions detected 
in 03-PCAdapt.  
It also subsets the same vcf file for the remaining neutral positions and applies a final 
filter for HWE.

Finally, the script converts the final adaptive and neutral .vcf files in .bed .bim .fam .raw
and .strct_in format necessary for downstream analyses.

The conversion to genepop format for use of GENODIVE (to calculate kinship) is done with
the PGDSpider GUI.   
input: neutral.vcf and population map, output: neutral.gen.txt  
outputted .gen.txt file: add information on first line (otherwise genodive won't recognise format)  

To run the script, set arguments:  
  $1 = input file (vcf)  
  $2 = species code  
  
#### for diplodus
```
mkdir 01-Diplodus
cp ../01-SNPfilters/01-Diplodus/dip_population_map_*.txt 01-Diplodus/
cp ../01-SNPfilters/01-Diplodus/dip_all_filtered.vcf 01-Diplodus/

bash outlier_positions.sh 01-Diplodus/dip_all_filtered.vcf dip
```

After HWE filter, 7655 neutral loci were retained.

#### for Mullus
```
mkdir 02-Mullus
cp ../01-SNPfilters/02-Mullus/mul_population_map_*.txt 02-Mullus/
cp ../01-SNPfilters/02-Mullus/mul_all_filtered.vcf 02-Mullus/

bash outlier_positions.sh 02-Mullus/mul_all_filtered.vcf mul
```
After HWE filter, 2462 neutral loci were retained.


#### clean up files to separate directories
```
mv dip_* 01-Diplodus/
mv mul_* 02-Mullus/
```

