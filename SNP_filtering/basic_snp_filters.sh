
# variant filters were performed per population and for the Admixed individuals separately

for((i=1;i<=10;i++)); do
  java -Djava.io.tmpdir=/path_to_tempdir -Xmx12g -jar GenomeAnalysisTK.jar \
    -R /path_to_reference/matina_v1.1.fasta \
    -T VariantFiltration \
    -L scaffold_10 \
    -V raw_snps.vcf \
    --filterExpression "QD < 2.0 || FS > 50.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \ 
    --filterName "my_snp_filter" \ 
    -o /path_to_pop_vcf/chr$i.filtered.population.vcf \

  bgzip /path_to_pop_vcf/chr$i.filtered.population.vcf
  tabix /path_to_pop_vcf/chr$i.filtered.population.vcf.gz
done

