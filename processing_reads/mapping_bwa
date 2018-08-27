# Initially the index for the matina genome was created (only once)

bwa index -a bwtsw path_to_reference/matina_v1.1.fasta

bwa aln -n 0.06 -M 6 path_to_reference/matina_v1.1.fasta path_to_reads/reads.trimmed.1.fastq > reads_1.sai

bwa aln -n 0.06 -M 6 path_to_reference/matina_v1.1.fasta path_to_reads/reads.trimmed.2.fastq > reads_2.sai

bwa sampe -r '@RG\tID:sample_ID.librarynumber\tSM:sample_ID\tPL:Illumina\tLB:library_sample_ID' path_to_reference/matina_v1.1.fasta reads_1.sai reads_2.sai path_to_reads/reads.trimmed.1.fastq path_to_reads/reads.trimmed.2.fastq > sample_ID.sam
