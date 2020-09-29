sample=$1

#3-1.Remove duplicate
java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar MarkDuplicates \
-I ${sample}.sorted.bam \
-O ${sample}.rmdup.bam \
-M ${sample}.rmdup.metrics \
--REMOVE_DUPLICATES=true

#3-2 Make index 
/BiO/apps/samtools/samtools index ${sample}.rmdup.bam

#4.Base Quality Score Recalibration
java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar BaseRecalibrator \
-R /BiO/data/reference/hg19.fasta -I ${sample}.rmdup.bam \
-L /BiO/data/target/target.bed --known-sites \
/BiO/data/DB/dbSnp151_chr.vcf.gz --known-sites \
/BiO/data/DB/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz \
-O ${sample}_recal_data.table

#5-1.Apply BQSR
java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar ApplyBQSR \
-bqsr ${sample}_recal_data.table \
-I ${sample}.rmdup.bam \
-O ${sample}.recal.bam

#5-2 Make index
/BiO/apps/samtools/samtools index ${sample}.recal.bam

#6.On-target coverage
java -Xmx8g -jar /BiO/apps/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
-T DepthOfCoverage -R /BiO/data/reference/hg19.fasta \
-I ${sample}.recal.bam -O ${sample}_target_cov \
-ct 1 -ct 5 -ct 10 -ct 20 -ct 30 -omitBaseOutput \
-L /BiO/data/target/target.bed

#7-1.Remove off-target reads
bedtools intersect -wa -a ${sample}.recal.bam -b /BiO/data/target/target.bed > ${sample}.target.bam

#7-2.Remove off-target reads
/BiO/apps/samtools/samtools index ${sample}.target.bam

#8-1-1.Variant Calling
java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar HaplotypeCaller \
-R /BiO/data/reference/hg19.fasta \
-I ${sample}.target.bam \
-O ${sample}.gatk.vcf

#8-1-2.Variant Calling
/BiO/apps/bcftools/bcftools mpileup -f \
/BiO/data/reference/hg19.fasta ${sample}.target.bam | \
/BiO/apps/bcftools/bcftools call -mv -Ov \
-o ${sample}.samt.vcf

#8-2-1.Make index
bgzip ${sample}.gatk.vcf
tabix -p vcf ${sample}.gatk.vcf.gz

#8-2-2.Make index
bgzip ${sample}.samt.vcf
tabix -p vcf ${sample}.samt.vcf.gz

#9-1.Consensus VCF 
vcf-isec -o -n +2 ${sample}.gatk.vcf.gz ${sample}.samt.vcf.gz > ${sample}.consensus.vcf

#10-1.Filteration
java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar VariantFiltration \
-R /BiO/data/reference/hg19.fasta \
-O ${sample}.consensus.filt.vcf --variant ${sample}.consensus.vcf \
--filter-expression 'DP < 10 || FS > 60.0' --filter-name 'LOWQUAL' 

#10-2.Filteration
cat ${sample}.consensus.filt.vcf | awk -F '\t' '($7!="LOWQUAL") {print}' | bgzip > ${sample}.final.vcf.gz

#10-3.Filteration
tabix -p vcf ${sample}.final.vcf.gz
