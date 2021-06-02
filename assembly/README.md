# Genome assembly

## Table of contents
1. [Raw data](#raw)
2. [FALCON assembly](#falcon)
3. [Decontamination](#decon)
4. [Polishing](#polish1)
5. [Haplotype resolution](#purge_dups)
6. [Manual curation](#mcur)
7. [Polishing](#polish2)

## 01 - Raw data <a name="raw"></a>
### Convert BAM to FASTA
```
# For each bam file index (For this examples I'm using 2/6 renamed bam files)
pbindex subreads1.bam
pbindex subreads2.bam

# Extract reads
bam2fasta -o subreads1.out subreads1.bam
bam2fasta -o subreads2.out subreads2.bam

# Combine FASTA files
cat subreads1.out.fasta subreads2.out.fasta > all.subreads.fasta

# Merge BAM files (for polishing later)
samtools merge all.merged.bam subreads1.bam subreads2.bam

# Index merged BAM (for polishing later)
samtools index all.merged.bam
pbindex all.merged.bam
```
## 02 - FALCON assembly <a name="falcon"></a>
### Assembly
```
# Create the primary assembly
fc_run falcon.cfg

# Unzip the assembly
fc_unzip.py falcon_unzip.cfg

# Rename contigs in primary and haplotig files and merge for decontamination
sed 's/>/>PRIM/g' cns_p_ctg.fasta
sed 's/>/>HAP/g' cns_h_ctg.fasta
cat cns_p_ctg.fasta cns_h_ctg.fasta > assembly.fasta
```
## 03 - Decontamination <a name="decon"></a>
### Assembly
```
# Run Diamond
diamond blastx \
 --query assembly.fasta \
 --db uniprot_ref_proteomes.diamond.dmnd \
 --outfmt '6' \
 --sensitive \
 --max-target-seqs 1 \
 --evalue 1e-25 \
 --out all.renamed.fa.diamond.tsv \
 --threads 10

# Taxify diamond hits
blobtools taxify -f all.renamed.fa.diamond.tsv -s 0 -t 1 -m list.taxid

# Run BLAST
blastn \
 -query assembly.fasta \
 -db nt \
 -outfmt '6 qseqid staxids bitscore std' \
 -max_target_seqs 10 \
 -max_hsps 1 \
 -num_threads 10 \
 -evalue 1e-25 \
 -out all.renamed.fa.blastn.tsv

# Map subreads to assembly
minimap2 -t 8 -ax map-pb assembly.fasta all.subreads.fasta  > aln.sam
samtools sort -@ 12 -O BAM -o aln.sort.bam aln.sam

# Create map file
blobtools map2cov -i all.renamed.fa -b aln.sort.bam -o aln.cov

# Create blotools DB
blobtools create -i assembly.fasta -o blobtools -c aln.cov.aln.sort.bam.cov -t all.renamed.fa.blastn.tsv -t all.renamed.fa.diamond.tsv.taxified.out -x bestsumorder --nodes nodesDB.txt

# Create blobplot
blobtools plot -i blobtools.blobDB.json -x bestsumorder --out blobtools_plot

# Identify the contaminated scaffolds (inspect the output for potential contaminants)
grep -v "Platy" blobtools_view.blobtools.blobDB.bestsumorder.table.txt

# Using a list of contaminant scaffolds
fastqualselect.pl -f assembly.fasta -e contaminant.list > assembly.clean.fasta
```
## 04 - Polishing <a name="polish1"></a>
### Polish assembly
```
# Align subreads
pbmm2 align -j 12 assembly.clean.fasta all.subreads.bam pbmm2.mapped.A0.bam

# Polish
gcpp -j 12 --referenceFilename assembly.clean.fasta -o assembly.clean.A1.fasta pbmm2.mapped.A0.bam

# Align subreads
pbmm2 align -j 12 assembly.clean.A1.fasta all.subreads.bam pbmm2.mapped.A1.bam

# Polish
gcpp -j 12 --referenceFilename assembly.clean.fasta -o assembly.clean.A2.fasta pbmm2.mapped.A1.bam
```
## 05 - Haplotype resolution <a name="purge_dups"></a>
### 
```
# Subset to get the primary assembly
grep PRIM assembly.clean.A2.fasta | sed 's/|//g' > primary.list
fastqualselect.pl -f assembly.clean.A2.fasta -i primary.list > assembly.clean.A2.primary.fasta
fastqualselect.pl -f assembly.clean.A2.fasta -i hap.list > assembly.clean.A2.haplotypes.fasta

# Purge duplicates, using an file with list of FASTA subreads
pd_config.py -l purge assembly.clean.A2.primary.fasta input.fofn
run_purge_dups_edit.py config.json purge_dups/src/ purge

# Get the primary and haplotype assemblies and redo
pd_config.py -l purge assembly.clean.A2.primary.purged.fa input.fofn 
cat assembly.clean.A2.haplotypes.fasta assembly.clean.A2.primary.hap.fa > all1.hap.fa
run_purge_dups_edit.py config.json purge_dups/src/ purge
```
## 06 - Manual curation <a name="mcur"></a>
### 
```
# Using a renamed haplotype purged assembly and preads from FALCON
minimap2 -t 12 -ax asm20 primary_assembly.fa falcon.preads.fasta | samtools sort -@ 12 -O BAM -o preads.sort.bam
tg_index -o preads.sort preads.sort.bam
```
## 07 - Polishing <a name="polish2"></a>
### Polish assembly
```
# Merge primary and haplotype assemblies
cat primary_assembly.curated.fasta haplotype_assembly.fasta > all.assembly.fasta

# Align subreads to curated assembly 
pbmm2 align -j 12 all.assembly.fasta all.subreads.bam pbmm2.mapped.A2.bam

# Polish
gcpp -j 12 --referenceFilename all.assembly.fasta -o all.assembly.A3.fasta pbmm2.mapped.A2.bam
```
