# data was merged prior to all following analyses per chromosome
vcf-merge pop1/chr1.filtered.vcf pop2/chr1.filtered.vcf ... |bgzip -c > combined/chr1.vcf.gz

# data was first thinned and MAF < 0.05 SNPs were removed 

vcf-concat combined/chr1.vcf.gz combined/chr2.vcf.gz combined/chr3.vcf.gz combined/chr4.vcf.gz combined/chr5.vcf.gz combined/chr6.vcf.gz combined/chr7.vcf.gz combined/chr8.vcf.gz combined/chr9.vcf.gz combined/chr10.vcf.gz |bgzip -c > combined/cacao.vcf.gz

cd combined/

vcftools --gzvcf cacao.vcf.gz --thin 5000 --maf 0.05 --min-alleles 2 --max-alleles 2 --remove-indels \
  --max-missing 1 --recode --out 4_structure/cacao.thin

de 4_structure/

# for Admixture

plink --file cacao.thin --recode12 --make-bed cacao.thin2
sed -i s'/^0/'$i'/g cacao.thin2

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ; do
  admixture --cv cacao.thin2.bed $K -l 500 -e 0.2 |tee log${K}.out ;
done

# results were analyzed in R and barplots showing results were plotted in R

#  for PCA we generated a new matrix containing the data recoded in 012

vcftools --gzvcf cacao.vcf.gz --thin 5000 --maf 0.05 --min-alleles 2 --max-alleles 2 --remove-indels \
  --max-missing 1 --recode --recode012 cacao.4pca
  
# PCA analyses were carried out using the cmdscale function in R. See plotting_structure.R


Fst estimations were performed in 5Kb windows using vcftools for all pairs of populations


vcftools --gzvcf cacao.vcf.gz --thin 5000 --maf 0.05 --min-alleles 2 --max-alleles 2 --remove-indels \
  --max-missing 1  --keep 4_fsts/pop1_pop2.ids --weir-fst-pop 4_fsts/pop1.ids --weir-fst-pop 4_fsts/pop2.ids 
  --out 4_fsts/pop1_pop2/cacao_pop1pop2.fst
  
