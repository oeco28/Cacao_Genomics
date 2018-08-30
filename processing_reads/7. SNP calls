# variant identification was performed per population and for the Admixed individuals separately

for((i=1;i<=10;i++)); do
  java -Djava.io.tmpdir=/path_to_tempdir -Xmx12g -jar GenomeAnalysisTK.jar \
    -R /path_to_reference/matina_v1.1.fasta \
    -T UnifiedGenotyper \
    -I path_to_final_bam/sampleID.recal.bam \
    -I path_to_final_bam/sampleID.recal.bam \
    ...
    .
    .
    -I path_to_final_bam/sampleID.recal.bam \
    -L scaffold_10 \
    --dbsnp /path_to_VQSR_ref/6kSNPs.vcf \
    --output_mode EMIT_ALL_SITES \
    --heterozygosity 0.006 \
    --indel_heterozygosity 1.25E-3 \
    -dcov 200 \
    -o /path_to_pop_vcf/chr$i.population.vcf \

done
