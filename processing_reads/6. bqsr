mkdir recalibration
mkdir final_bam/

samtools index loc_realigned/JA_5_5.realigned.bam
java -Xmx6g -jar GenomeAnalysisTK.jar -R /path_to_reference/matina_v1.1.fasta \
  -T BaseRecalibrator \
  -I loc_realigned/JA_5_5.realigned.bam \
  -o recalibration/JA_5_5.recal_data.grp 
  --knownSites /path_to_VQSR_ref/6kSNPs.vcf 
  --covariate ContextCovariate \
  --covariate RepeatLengthCovariate \
  --covariate CycleCovariate \
  --validation_strictness SILENT


java -Xmx6g -jar GenomeAnalysisTK.jar -R /path_to_reference/matina_v1.1.fasta \
  -T PrintReads \
  -I loc_realigned/JA_5_5.realigned.bam \
  -BQSR recalibration/JA_5_5.recal_data.grp \
  -o final_bam/JA_5_5.recal.bam 
  --validation_strictness SILENT

