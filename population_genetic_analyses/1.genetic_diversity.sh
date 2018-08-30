# Estimation of genetic diversity (Watterson's Theta and pi) was done with vcftools and custom scripts per chromosome

for((i=1;i<=10;i++)); do
  vcftools --gzvcf chr$i.population.filtered.vcf.gz --remove-filtered-all --chr scaffold_$i --window-pi 1000 --out chr$i.1Kb.windowed.pi
  vcftools --gzvcf chr$i.population.filtered.vcf.gz --remove-filtered-all --chr scaffold_$i --SNPdensity 1000 --out chr$i.theta.1Kb
done

# Figures and downstream analyses were generated with custom scripts
