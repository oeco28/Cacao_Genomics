for((i=1;i<=10;i++)); do
vcftools --gzvcf cacao.filtered.vcf.gz --bed /path_to_intergenic_beds/chr$i.subset2.intergenic.bed \
  --chr scaffold_$i --min-alleles 2 --max-alleles 2 --remove-indels --max-missing 1 \
  --recode --out 4_demo/chr$i.intergenic
done
