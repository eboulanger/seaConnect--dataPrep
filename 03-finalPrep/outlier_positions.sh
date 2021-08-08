## Final steps to get neutral and adaptive SNP set and correct file formats

# define arguments

inputFile=$1
spCode=$2

#### adaptive SNPs ####

# step 1: rename list of outlier positions
cp ../03-PCAdapt/outl_pos_pcadpt_"$spCode".txt "$spCode"_outl_pos.txt
# extract total number of unique outlier positions (for file naming)
nPosOut=`sort "$spCode"_outl_pos.txt | uniq -u | wc -l` 
echo ""$nPosOut" unique outlier loci in total"
# remove leading whitespaces (from here: https://www.cyberciti.biz/faq/bash-remove-whitespace-from-string/)
shopt -s extglob # turn it on
nPosOut="${nPosOut##*( )}"
shopt -u extglob # turn it off
# add number of outliers to file name
mv "$spCode"_outl_pos.txt "$spCode"_outl_pos_"$nPosOut".txt

## step 2: subset filtered vcf file by outlier positions

vcftools --vcf "$inputFile" --positions "$spCode"_outl_pos_"$nPosOut".txt --recode --recode-INFO-all --out "$spCode"_adaptive_"$nPosOut"
mv "$spCode"_adaptive_"$nPosOut".recode.vcf "$spCode"_adaptive_"$nPosOut".vcf

## step 3 : convert vcf files to necessary formats for downstream analyses

# convert .vcf files to PLINK 1 binary files
# writes .bed .bim and .fam files (ADMIXTURE will need .bed file)
plink --vcf "$spCode"_adaptive_"$nPosOut".vcf --out "$spCode"_adaptive_"$nPosOut" --allow-extra-chr --vcf-half-call missing

# convert to minor allele frequencies for RDA and other analyses
# writes .raw file
plink --bfile "$spCode"_adaptive_"$nPosOut" --recodeA --out "$spCode"_adaptive_"$nPosOut" --allow-extra-chr

# convert to STRUCTURE format for assignPop 
# writes .strct_in file
plink --bfile "$spCode"_adaptive_"$nPosOut" --recode structure --out "$spCode"_adaptive_"$nPosOut" --allow-extra-chr

#### neutral SNPs ####

## step 1: create a final list of all neutral loci positions

# get list of all positions original vcf file
grep -v "^##" "$inputFile" | cut -f1-2 | sed '1d' > "$spCode"_all_pos.txt

# remove previously identified outlier positions to only retain neutral ones
cat "$spCode"_all_pos.txt "$spCode"_outl_pos_"$nPosOut".txt | sort | uniq -u > "$spCode"_ntrl_pos_preHWE.txt

## step 2 : subset filtered vcf file by neutral outlier positions
vcftools --vcf "$inputFile" --positions "$spCode"_ntrl_pos_preHWE.txt --recode --recode-INFO-all --out "$spCode"_neutral_preHWE

## step 3 : apply HWE filter
vcftools --vcf "$spCode"_neutral_preHWE.recode.vcf --hwe 0.000001 --recode --recode-INFO-all --out "$spCode"_neutral
nPosNtrl=`grep  -c -v "^#" "$spCode"_neutral.recode.vcf`
mv "$spCode"_neutral.recode.vcf "$spCode"_neutral_"$nPosNtrl".vcf
mv "$spCode"_neutral.log "$spCode"_neutral_"$nPosNtrl".log

## step 4 : convert vcf files to necessary formats for downstream analyses

# convert .vcf files to PLINK 1 binary files
# writes .bed .bim and .fam files (ADMIXTURE will need .bed file)
plink --vcf "$spCode"_neutral_"$nPosNtrl".vcf --out "$spCode"_neutral_"$nPosNtrl" --allow-extra-chr --vcf-half-call missing

# convert to minor allele frequencies for RDA and other analyses
# writes .raw file
plink --bfile "$spCode"_neutral_"$nPosNtrl" --recodeA --out "$spCode"_neutral_"$nPosNtrl" --allow-extra-chr

# convert to STRUCTURE format for assignPop 
# writes .strct_in file
plink --bfile "$spCode"_neutral_"$nPosNtrl" --recode structure --out "$spCode"_neutral_"$nPosNtrl"  --allow-extra-chr

# convert to genepop format for use of GENODIVE (to calculate kinship)
# use PGDSpider GUI
# input: neutral.vcf and population map, output: neutral.gen.txt
# outputted .gen.txt file: add information on first line (otherwise genodive won't recognise format)


#### all SNPs ####

## convert vcf files to necessary formats for downstream analyses

# convert .vcf files to PLINK 1 binary files
# writes .bed .bim and .fam files (ADMIXTURE will need .bed file)
plink --vcf "$inputFile" --out "$spCode"_all_filtered --allow-extra-chr --vcf-half-call missing

# convert to minor allele frequencies for RDA and other analyses
# writes .raw file
plink --bfile "$spCode"_all_filtered --recodeA --out "$spCode"_all_filtered --allow-extra-chr

# convert to STRUCTURE format for assignPop 
# writes .strct_in file
plink --bfile "$spCode"_all_filtered --recode structure --out "$spCode"_all_filtered --allow-extra-chr
