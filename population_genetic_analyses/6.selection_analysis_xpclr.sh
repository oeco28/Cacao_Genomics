# . for each chromosome
# data was filtered and generated

for((i=1;i<=10;i++)); do 
        vcftools --vcf cacao.snps.vcf --chr scaffold_$i --keep ./4_demo/Criollo_Curaray.ids \
        --min-alleles 2 --max-alleles 2 --max-maf 0.9 --maf 0.09 --max-missing 1 --recode \
        4_demo/Criollo_Curaray/chr$i.FILTERED ; done

# all data was reformated to meet the expected input for XPCLR.

for((i=2;i<=10;i++)); do
        ./XPCLR -xpclr criollo.chr$i curaray.chr$i chr$i.map cacao.chr$i.xpclr.out -w1 0.005 200 1000 1 -p0 0.95
;
done
