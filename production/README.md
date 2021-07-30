# Mapping, variant calling and quality control

## Table of contents
1. [Raw data](#raw)
2. [Mapping](#mapping)
3. [Variant calling](#variantcalling)
4. [Quality control](#qc)

## 01 - Raw data <a name="raw"></a>
### Reference genome
```
# Remove unplaced contigs
samtools faidx tdSchCurr1.primary.fa
grep -v 'UNPLACED' tdSchCurr1.primary.fa.fai | cut -f1 > include.list
seqtk subseq tdSchCurr1.primary.fa include.list > tdSchCurr1.chrom.fa

# Create indexes and a sequence dictionary for the reference genome
bwa index tdSchCurr1.chrom.fa
gatk CreateSequenceDictionary --REFERENCE tdSchCurr1.chrom.fa
```
## 02 - Mapping <a name="mapping"></a>
### Map sequence reads to reference genome
```
# Repeat as needed for each read set, for example:
bwa mem -t 6 tdSchCurr1.chrom.fa SAMPLE1_1.fastq.gz SAMPLE1_2.fastq.gz | samtools sort -@6 -o SAMPLE1.bam -
```
### Mark PCR duplicates
```
# Mark duplicates
gatk MarkDuplicates --INPUT SAMPLE1.bam --OUTPUT SAMPLE1.markdup.bam --METRICS_FILE SAMPLE1.metrics.txt

# Merge BAM files (from samples where there are multiple sets of FASTQs)
samtools merge -@ 6 SAMPLE1.markdup.merged.bam SAMPLE1.markdup.bam SAMPLE1b.markdup.bam

# Index all BAMs
samtools index SAMPLE1.markdup.merged.bam
```
### Calculate coverage
```
# Create makewindows input
cut -f1,2 tdSchCurr1.chrom.fa.fai > tdSchCurr1.chrom.txt

# Create 25 kb windows
bedtools makewindows -g tdSchCurr1.chrom.txt -w 25000 > tdSchCurr1.chrom.25kb.bed

# Calculate per-sample coverage
bedtools coverage -sorted -g tdSchCurr1.chrom.fa.fai -d -a tdSchCurr1.chrom.25kb.bed -b SAMPLE1.markdup.merged.bam \| datamash -g1,2,3 median 5 mean 5 sstdev 5 > SAMPLE1.cov
awk '{print $1,$2,$3,$4,$5,$6,FILENAME}' SAMPLE1.cov | sed 's/.cov//g' > SAMPLE1.recov
cat *.recov > all.recov.txt
```
## 03 - Variant calling <a name="variantcalling"></a>

### Per-sampling variant calling
```
samtools index SAMPLE1.markdup.merged.bam
gatk HaplotypeCaller --emit-ref-confidence GVCF -I SAMPLE1.markdup.merged.bam -R tdSchCurr1.chrom.fa -O SAMPLE1.g.vcf
```
### Combine all samples into a single gVCF and genotype
```
# Combine gVCFs
ls | grep 'g.vcf' > argument.list
gatk CombineGVCFs --arguments_file argument.list --reference tdSchCurr1.chrom.fa --output merged_all_samples.g.vcf

# Genotype
gatk GenotypeGVCFs --reference tdSchCurr1.chrom.fa --variant merged_all_samples.g.vcf --output merged_all_samples.vcf
```
## 04 - Quality control <a name="qc"></a>
### Calculate quality scores for all variant sites
```
# Produce a table of quality scores for each variant site
gatk VariantsToTable --variant merged_all_samples.vcf -F CHROM -F POS -F TYPE -F QD -F FS -F MQ -F MQRankSum -F ReadPosRankSum -F SOR -F InbreedingCoeff -R tdSchCurr1.chrom.fa --output cohort.genotyped.tbl
```
### Separate and filter SNPs
```
# Select SNPs
gatk SelectVariants -R tdSchCurr1.chrom.fa --variant merged_all_samples.vcf --select-type-to-include SNP --output merged_all_samples.SNPs.vcf

# Tag low-quality SNPs
gatk VariantFiltration \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS8" \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "FS > 60.0" --filter-name "FS60" \
--filter-expression "MQ < 40.0" --filter-name "MQ40" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQ12.5" \
--filter-expression "SOR > 3.0" --filter-name "SOR3" \
--variant merged_all_samples.SNPs.vcf \
-R tdSchCurr1.chrom.fa  \
--output merged_all_samples.SNPs.tagged.vcf

# Remove low-quality sites
gatk SelectVariants -R tdSchCurr1.chrom.fa --variant merged_all_samples.SNPs.tagged.vcf --exclude-filtered --output merged_all_samples.SNPs.filtered.vcf
```
### Separate and filter indels and mixed sites
```
# Select indels and mixed sites
gatk SelectVariants -R tdSchCurr1.chrom.fa --variant merged_all_samples.vcf --select-type-to-include SNP --output merged_all_samples.indels_mixed.vcf

# Tag low-quality indels and mixed sites
/lustre/scratch118/infgen/team133/db22/software/gatk-4.1.0.0/gatk VariantFiltration \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "FS > 200.0" --filter-name "FS200" \
--filter-expression "ReadPosRankSum < -20.0" --filter-name "RPRS20" \
--filter-expression "SOR > 10.0" --filter-name "SOR10" \
--variant merged_all_samples.indels_mixed.vcf \
-R tdSchCurr1.chrom.fa  \
--output merged_all_samples.indels_mixed.tagged.vcf

# Remove low-quality sites
gatk SelectVariants -R tdSchCurr1.chrom.fa --variant merged_all_samples.indels_mixed.tagged.vcf --exclude-filtered --output merged_all_samples.indels_mixed.filtered.vcf
```
### Recombine filtered variants
```
gatk MergeVcfs --INPUT merged_all_samples.SNPs.filtered.vcf --INPUT merged_all_samples.indels_mixed.filtered.vcf --OUTPUT merged_all_samples.filtered.vcf
```
### Remove low-quality samples and variants 
```
# Calculate per-individual missingness rate 
vcftools --vcf merged_all_samples.filtered.vcf.FL1.vcf --missing-indv --out missing_indv

# Filter out individuals with high rates of missing variant calls
awk '$6<0.55' missing_site.imiss | grep -v "MISS" | cut -f1  > retain.samples.list
vcftools --vcf merged_all_samples.filtered.vcf.FL1.vcf --keep retain.samples.list --recode-INFO-all --recode --out merged_all_samples.filtered.vcf.FL1.vcf  

# Calculate per-site missingness rate 
vcftools --vcf merged_all_samples.filtered.vcf --missing-site --out missing_site

# Filter out sites with high rates of missing variant calls
awk '$6<=0.05' missing_site.lmiss | grep -v "MISS" | cut -f1,2 > retain.variants.list
vcftools --vcf merged_all_samples.filtered.vcf.FL1.vcf --postions retain.variants.list --recode-INFO-all --recode --out merged_all_samples.filtered.vcf.FL2.vcf 
```
### Move final versions of VCFs to the analysis folder
```
mv merged_all_samples.filtered.vcf.FL2.vcf FREEZE.FULLFILTER.vcf
```
