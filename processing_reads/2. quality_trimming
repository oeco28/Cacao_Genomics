# The preamble of the scripts is not provided. Users should 

module load python/2.7.10

my_trim_galore="path_to_trim_galore"


$my_trim_galore --output_dir path_to_Trimmed_data_out \
                  --quality 28 \
                  --illumina \
                  --phred33 \
                  --fastqc_args "--nogroup --noextract --outdir ../fastqc" \
                  --stringency 6 \
                  --length 60 \
                  --clip_R1 10 \
                  --clip_R2 10 \
                  --paired \
                  path_to_fastq_files/read_1.fastq.gz \
                  path_to_fastq_files/read_2.fastq.gz
