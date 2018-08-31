#######
# psmc
########

samtools mpileup -C50 -uf /path_to_reference/matina_v1.1.fasta sampleID.recal.bam sampleID.2.recal.bam .. \ 
  |bcftools call -c - \
  |vcfutils.pl vcf2fq -d 10 -D 100 | gzip > ../PSMC/sampleID.fq.gz

fq2psmcfa -q20 sampleID.fq.gz > sampleID.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid.psmc diploid.psmcfa
psmc2history.pl diploid.psmc | utils/history2ms.pl > ms-cmd.sh
psmc_plot.pl diploid diploid.psmc

then we used summarized results across individuals from the same population to generate the Figures of the paper with plotting_psmc.R)

##########
#  SMC++
#  for each population we followed the default conditions given for smc++ and generated results similar to those presented here.
###########


mkdir analysis3
for((i=1;i<=10;i++)); do
        smc++ vcf2smc -c 80000 -d Catongo MATINA_Tica2 chr$i.recode.vcf.gz out_cl/chr$i.1.smc.gz scaffold_$i Amelonado:Catongo,MATINA_Tica2,Matina,REDAMEL_1_27,REDAMEL_1_31,SIAL169,SIAL70,SIAL84,SIC806,TRD86,mvP30 ;
        smc++ vcf2smc -c 80000 -d Matina REDAMEL_1_27 chr$i.recode.vcf.gz out_cl/chr$i.2.smc.gz scaffold_$i Amelonado:Catongo,MATINA_Tica2,Matina,REDAMEL_1_27,REDAMEL_1_31,SIAL169,SIAL70,SIAL84,SIC806,TRD86,mvP30 ;
        smc++ vcf2smc -c 80000 -d REDAMEL_1_31 SIAL169 chr$i.recode.vcf.gz out_cl/chr$i.3.smc.gz scaffold_$i Amelonado:Catongo,MATINA_Tica2,Matina,REDAMEL_1_27,REDAMEL_1_31,SIAL169,SIAL70,SIAL84,SIC806,TRD86,mvP30 ;
        smc++ vcf2smc -c 80000 -d SIAL70 SIAL84 chr$i.recode.vcf.gz out_cl/chr$i.4.smc.gz scaffold_$i Amelonado:Catongo,MATINA_Tica2,Matina,REDAMEL_1_27,REDAMEL_1_31,SIAL169,SIAL70,SIAL84,SIC806,TRD86,mvP30 ;
        smc++ vcf2smc -c 80000 -d SIAL84 SIC806 chr$i.recode.vcf.gz out_cl/chr$i.5.smc.gz scaffold_$i Amelonado:Catongo,MATINA_Tica2,Matina,REDAMEL_1_27,REDAMEL_1_31,SIAL169,SIAL70,SIAL84,SIC806,TRD86,mvP30 ;
        smc++ vcf2smc -c 80000 -d TRD86 mvP30 chr$i.recode.vcf.gz out_cl/chr$i.6.smc.gz scaffold_$i Amelonado:Catongo,MATINA_Tica2,Matina,REDAMEL_1_27,REDAMEL_1_31,SIAL169,SIAL70,SIAL84,SIC806,TRD86,mvP30 ;
done

smc++ estimate -o analysis3/ 7.16e-9 out_cl/chr*smc.gz
smc++ plot amelonado_7.1e9.cl.pdf --csv analysis3/model.final.json
