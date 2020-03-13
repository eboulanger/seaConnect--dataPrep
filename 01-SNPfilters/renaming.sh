# define arguments

FOLDER=$1
SPCODE=$2

cd $FOLDER

# extract sample / individual IDs
bcftools query -l "$SPCODE"_all_filtered_origid.vcf > id_"$SPCODE".txt

# rename in R and create population map (necessary for PGDSpider conversions later on)
Rscript --vanilla ../rename_id.R "$SPCODE"

# replace sample names in vcf file
bcftools reheader "$SPCODE"_all_filtered_origid.vcf --samples new_id_"$SPCODE".txt > "$SPCODE"_all_filtered.vcf