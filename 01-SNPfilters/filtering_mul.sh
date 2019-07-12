#### Filtering steps by dDocent
#           tutorial: https://www.ddocent.com/filtering/

cd ~/Documents/project_SEACONNECT/seaConnect--dataPrep/01-SNPfilters/

mkdir 02-Mullus
cd 02-Mullus

FILE=../../00-rawData/02-Mullus/mullus.vcf
# PATH=$PATH:~/programmes/vcflib/bin/

# step 1: filter genotypes called below 50% (across all individuals), 
#         filter SNPs that have a minor allele count less than 3
#         filter SNPs with a minimum quality score of 30

vcftools --vcf $FILE --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out mullus.g5mac3

# step 2: minimum depth for a genotype call
#         minimum mean depth

vcftools --vcf mullus.g5mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out mullus.g5mac3dp3 

# step 3: asses individual levels of missing data and remove individuals accordingly

vcftools --vcf mullus.g5mac3dp3.recode.vcf --missing-indv
cat out.imiss

# make histogram of percentage missing data

awk '!/IN/' out.imiss | cut -f5 > totalmissing
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

# create list of individuals with more than 50% missing data

awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv

# remove individuals from list

vcftools --vcf mullus.g5mac3dp3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out mullus.g5mac3dplm

# step 4 : restrict data to variants called in a high percentage of individuals and filter by mean depth of genotypes

vcftools --vcf mullus.g5mac3dplm.recode.vcf --max-missing 0.95 --maf 0.05 --recode --recode-INFO-all --out DP3g95maf05 --min-meanDP 0

###
# FreeBayes - specific filters

# look at header
awk '/#/' DP3g95maf05.recode.vcf

# step 5 : filter for allele balance
#          this requires vcflib

vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" DP3g95p5maf05.recode.vcf > DP3g95p5maf05.fil1.vcf
vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" DP3g95maf05.recode.vcf > DP3g95maf05.fil1.vcf
# check how many were removed
awk '!/#/' DP3g95maf05.recode.vcf | wc -l
awk '!/#/' DP3g95maf05.fil1.vcf | wc -l

# step 6 : filter out sites that have reads from both strands
vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s DP3g95maf05.fil1.vcf > DP3g95maf05.fil2.vcf
awk '!/#/' DP3g95maf05.fil2.vcf | wc -l

# step 7 : look at ration of mapping qualities between reference and alternate alleles 
vcffilter -f "MQM / MQMR > 0.9 & MQM / MQMR < 1.05" DP3g95maf05.fil2.vcf > DP3g95maf05.fil3.vcf
awk '!/#/' DP3g95maf05.fil3.vcf | wc -l

# step 8 : check whether there is discrepancy in the properly paired status (SNAP IK NIET MAAR MAAKT BLIJKBAAR NIET UIT WANT GEEN SNPs WEGGEFILTERD)
vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s DP3g95maf05.fil3.vcf > DP3g95maf05.fil4.vcf
awk '!/#/' DP3g95maf05.fil4.vcf | wc -l

# step 9 : remove any locus that has a quality score below 1/4 of the depth
#          (read suggested blogposts to understand why) -> removes a lot of SNPs so maybe skip?
vcffilter -f "QUAL / DP > 0.25" DP3g95maf05.fil4.vcf > DP3g95maf05.fil5.vcf
awk '!/#/' DP3g95maf05.fil5.vcf | wc -l

# filter 10 : consists of multiple steps 
## 1 : create list of depth of each locus
cut -f8 DP3g95maf05.fil5.vcf | grep -oe "DP=[0-9]*" | sed 's/DP=//g' > DP3g95maf05.fil5.DEPTH
## 2 : create list of quality scores
awk '!/#/' DP3g95maf05.fil5.vcf | cut -f1,2,6 > DP3g95maf05.fil5.vcf.loci.qual
## 3 : calculate mean depth
MEAN_DP=`awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' DP3g95maf05.fil5.DEPTH`
## 4 : take the mean depth plus 3X the square of the mean 
MEAN_DP_SQ=`python -c "print int($MEAN_DP+3*($MEAN_DP**0.5))"`
echo $MEAN_DP_SQ
## 5 : paste depth and quality files together and find the loci above the cutoff that do not have quality scores 2 times the depth
paste DP3g95maf05.fil5.vcf.loci.qual DP3g95maf05.fil5.DEPTH | awk -v x=$MEAN_DP_SQ '$4 > x' | awk '$3 < 2 * $4' > DP3g95maf05.fil5.lowQDloci
## 6 : remove those sites and recalculate depth across loci
vcftools --vcf DP3g95maf05.fil5.vcf --site-depth --exclude-positions DP3g95maf05.fil5.lowQDloci --out DP3g95maf05.fil5
## 7 : cut depth scores, calculate average depth and plot
NINDIV=297
cut -f3 DP3g95maf05.fil5.ldepth > DP3g95maf05.fil5.site.depth
awk '!/D/' DP3g95maf05.fil5.site.depth | awk -v x=$NINDIV '{print $1/x}' > meandepthpersite
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
set xrange [10:150] 
unset label
set title "Histogram of mean depth per site"
set ylabel "Number of Occurrences"
set xlabel "Mean Depth"
binwidth=1
bin(x,width)=width*floor(x/width) + binwidth/2.0
set xtics 5
plot 'meandepthpersite' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
                                                                                                
## 8 : combine all to get new vcf
vcftools --vcf  DP3g95maf05.fil5.vcf --recode-INFO-all --out DP3g95maf05.FIL --max-meanDP 102.5 --exclude-positions DP3g95maf05.fil5.lowQDloci --recode 
# > After filtering, kept 297 out of 297 Individuals
# > Outputting VCF file...
# > After filtering, kept 5722 out of a possible 8871 Sites
# > Run Time = 6.00 seconds

# step 10 : HWE
#           recommends by pop but since I have mostly one large pop: apply to all

vcftools --vcf DP3g95maf05.FIL.recode.vcf --hwe 0.000001 --recode --recode-INFO-all --out DP3g95maf05.FIL.hwe

# step 11 : He

TO DO

# step 12 : Fis

TO DO

# step 3 : outliers

TO DO




