## Final steps to get neutral and adaptive SNP set and correct file formats

# define arguments

inputFile=$1
spCode=$2

#### adaptive SNPs ####

## step 1: create a final list of all outlier loci positions

cat ../02-Bayescan/outl_pos_bayesc_"$spCode".txt ../03-PCAdapt/outl_pos_pcadpt_"$spCode".txt > "$spCode"_outl_pos_combi.txt

# how many are detected by both methods?
nPosOutDup=`sort "$spCode"_outl_pos_combi.txt | uniq -d | wc -l` 
echo ""$nPosOutDup" outlier loci detected with both methods"
# extract total number of unique outlier positions (for file naming)
nPosOut=`sort "$spCode"_outl_pos_combi.txt | uniq -u | wc -l` 
echo ""$nPosOut" unique outlier loci in total"
# remove leading whitespaces (from here: https://www.cyberciti.biz/faq/bash-remove-whitespace-from-string/)
shopt -s extglob # turn it on
nPosOut="${nPosOut##*( )}"
shopt -u extglob # turn it off

# extract file with unique positions
sort "$spCode"_outl_pos_combi.txt | uniq -u > "$spCode"_outl_pos_"$nPosOut".txt
#rm "$spCode"_outl_pos_combi.txt

## step 2: subset filtered vcf file by outlier positions

vcftools --vcf "$inputFile" --positions "$spCode"_outl_pos_"$nPosOut".txt --recode --recode-INFO-all --out "$spCode"_adaptive_"$nPosOut"
mv "$spCode"_adaptive_"$nPosOut".recode.vcf "$spCode"_adaptive_"$nPosOut".vcf

## step 3 : convert vcf files to necessary formats for downstream analyses

# convert .vcf files to .bep files for ADMIXTURE
# .vcf to .tped
vcftools --vcf "$spCode"_adaptive_"$nPosOut".vcf --plink-tped --out "$spCode"_adaptive_"$nPosOut"
# .tped to .bed
plink --tped "$spCode"_adaptive_"$nPosOut".tped --tfam "$spCode"_adaptive_"$nPosOut".tfam --make-bed --out "$spCode"_adaptive_"$nPosOut" --noweb

# convert to minor allele frequencies for RDA and other analyses
# .tped to .raw 
plink --tped "$spCode"_adaptive_"$nPosOut".tped --tfam "$spCode"_adaptive_"$nPosOut".tfam --recodeA --out "$spCode"_adaptive_"$nPosOut" --noweb

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

# convert .vcf files to .bep files for ADMIXTURE
# .vcf to .tped
vcftools --vcf "$spCode"_neutral_"$nPosNtrl".vcf --plink-tped --out "$spCode"_neutral_"$nPosNtrl"
# .tped to .bed
plink --tped "$spCode"_neutral_"$nPosNtrl".tped --tfam "$spCode"_neutral_"$nPosNtrl".tfam --make-bed --out "$spCode"_neutral_"$nPosNtrl" --noweb

# convert to minor allele frequencies for RDA and other analyses
# .tped to .raw 
plink --tped "$spCode"_neutral_"$nPosNtrl".tped --tfam "$spCode"_neutral_"$nPosNtrl".tfam --recodeA --out "$spCode"_neutral_"$nPosNtrl" --noweb
