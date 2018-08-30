# genome index was created only once for all the tools (samtools, picard, etc)

samtools faidx path_to_reference/matina_v1.1.fasta
java  -Xmx8G -jar picard.jar CreateSequenceDictionary R=path_to_reference/matina_v1.1.fasta \
  O=path_to_reference/matina_v1.1.dict

samtools view -ubhSt path_to_reference/matina_v1.1.fasta.fai sample_ID.sam -o sample_ID.bam

java -Xmx6G -jar picard.jar CleanSam INPUT=sample_ID.bam \ 
  OUTPUT=sample_ID.clean.bam \ 
  VALIDATION_STRINGENCY=SILENT

java -Xmx6G -jar picard.jar FixMateInformation INPUT=sample_ID.clean.bam \
  OUTPUT=sample_ID.clean.fix.bam VALIDATION_STRINGENCY=SILENT

java -Xmx6G -jar picard.jar ValidateSamFile INPUT=sample_ID.clean.fix.bam \
  OUTPUT=sample_ID.validation VALIDATION_STRINGENCY=LENIENT

java -Xmx6G -jar picard.jar CollectAlignmentSummaryMetrics INPUT=sample_ID.clean.fix.bam \
  OUTPUT=sample_ID.metrics VALIDATION_STRINGENCY=LENIENT

java -Xmx4G -jar picard.jar SortSam INPUT=sample_ID.clean.fix.bam OUTPUT=sample_ID.clean.fix.sorted.bam \
SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT

java -Xmx4G -jar picard.jar MarkDuplicates INPUT=sample_ID.clean.fix.sorted.bam \
  OUTPUT=sample_ID.F.bam \
  VALIDATION_STRINGENCY=LENIENT \
  REMOVE_DUPLICATES=TRUE METRICS_FILE=sample_ID.DUP_METRICS.OUT

# If necessary we would correct issues with @RG information using the following two commands:
#  java -Xmx6G -jar picard.jar AddOrReplaceReadGroups INPUT=sample_ID.F.bam \ 
#  OUTPUT=sample_ID.F_2.bam SORT_ORDER=coordinate \
#  RGPL=Illumina RGLB=sample_ID RGPU=L1 RGSM=sample_ID_new RGID=sample_ID_new.number_library
# mv sample_ID.F_2.bam sample_ID.F.bam

samtools index sample_ID.F.bam

samtools idxstats sample_ID.F.bam > sample_ID.IDXSTATS
bamtools stats -in sample_ID.F.bam > sample_ID.STATS
