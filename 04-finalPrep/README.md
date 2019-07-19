## final steps to get neutral and adaptive SNP set

# adaptive SNPs

### step 1: create a final list of all outlier loci positions
```
cat ../02-Bayescan/outl_pos_bayesc_dip.txt ../03-PCAdapt/outl_pos_pcadpt_dip.txt > outl_pos_combi_dip.txt
uniq -d outl_pos_combi_dip.txt | wc -l
```
no diplodus outliers were detected by both methods

keep only unique positions:
```
uniq -u outl_pos_combi_dip.txt > outl_pos_dip_433.txt
```

### step 2: subset vcf file by outlier loci positions
```
vcftools --vcf ../01-SNPfilters/dip_all_filtered.vcf --positions outl_pos_dip_433.txt --recode-INFO-all --out diplodus_adaptive
mv diplodus_adaptive.recode.vcf diplodus_adaptive.vcf
```

# neutral SNPs

### step 1 : remove outlier loci

### step 2 : filter for HWE

step 11 : HWE
   recommends by pop but since I have mostly one large pop: apply to all

```
vcftools --vcf dip_all_filtered.vcf --hwe 0.000001 --recode --recode-INFO-all --out mul.FIL.HFis.hwe
```
