mkdir loc_realigned
java -Xmx6g -jar GenomeAnalysisTK.jar -R /path_to_reference/matina_v1.1.fasta \
  -T RealignerTargetCreator \
  -I JA_5_5.F.bam \
  -o loc_realigned/JA_5_5.forIndelRealigner.intervals

java -Xmx6g -jar GenomeAnalysisTK.jar -R /path_to_reference/matina_v1.1.fasta \
  -T IndelRealigner \
  -targetIntervals loc_realigned/JA_5_5.forIndelRealigner.intervals \
  -I JA_5_5.F.bam \
  -o loc_realigned/JA_5_5.realigned.bam

