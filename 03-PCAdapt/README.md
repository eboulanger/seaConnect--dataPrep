# detect outliers with PCAdapt in R and retrieve loci positions from vcf file

### Run the PCAdapt.R script

This script will load your vcf file, detect outlier loci using the package pcadapt and
output a list of outlier loci positions that can be used to subset your vcf file

set arguments : argument 1 = input file (vcf), argument 2 = species code

#### for mullus :
```
Rscript --vanilla PCAdapt.R ../01-SNPfilters/02-Mullus/mul_all_filtered.vcf mul
```
see how many outliers were detected :
```
wc -l outl_pos_pcadpt_mul.txt
```
2078 outliers detected for mullus

#### for diplodus :
```
Rscript --vanilla PCAdapt.R ../01-SNPfilters/01-Diplodus/dip_all_filtered.vcf dip

wc -l outl_pos_pcadpt_dip.txt
```
312 outliers detected for diplodus
