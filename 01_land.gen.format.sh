# author: Leonardo Caproni
# date: 08/2022
#-------------------------------------------------------------------------------
# Reference: DOI
#-------------------------------------------------------------------------------
# Description: 	Get bioclimatic variables at sampling points for both historical 
#				and projected climate data. It prepares genotypes matadata for 
#				downstram anlyses
#===============================================================================

##SNP pruning
#mkdir pruned_100_5_05
plink --vcf barley.snps.noWild.FINAL.AF005.436 --allow-extra-chr --make-founders --indep-pairwise 150 5 0.5 --out barley
plink --vcf barley.snps.noWild.FINAL.AF005.436  --allow-extra-chr --exclude barley.prune.out --recode vcf --out barley.snps.noWild.FINAL.AF005.436.pruned_150_5_05

# Format for GF and RDA
vcftools --vcf  barley.snps.noWild.FINAL.AF005.436.pruned_150_5_05.vcf --012 --out snp

cut -f2- snp.012 | sed 's/-1/NA/g' >snp.temp
tr -d '\t' <snp.012.pos | tr '\n' '\t' | sed 's/[[:space:]]*$//' >header
paste <(echo "ID" | cat - snp.012.indv) <(echo "" | cat header - snp.temp) > snp.pruned.forGF
rm header snp.temp

