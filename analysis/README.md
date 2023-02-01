# Analysis

## Table of contents
1. [Assembly QC](#AQC)
2. [Population Structure](#pca)
3. [Admixture statistics](#admix)
4. [Windowed admixture](#admix2)
5. [Mitochondrial genome analysis](#mito)
6. [Fixed differences](#fixed)
7. [Contamination QC](#cont)
8. [Repeat masking](#rep)

## 01 - Assembly QC <a name="AQC"></a>
### Run BUSCO
```
run_BUSCO.py -i tdSchCurr1.chrom.fa -o busco_euk_09 -l lineage_datasets/eukaryota_odb9/ -m genome -c 8 --long -t euk_09
```
### Hi-C plots
```
# Made indices 
bwa index tdSchCurr1.primary.renamed.fa
python generate_site_positions.py Arima tdSchCurr1 tdSchCurr1.primary.renamed.fa
awk '{print $1,$NF}' tdSchCurr1_Arima.txt > tdSchCurr1.CHROM.sizes

# Run Juicer
juicer.sh -g tdSchCurr1 -s Arima -y tdSchCurr1_Arima.txt -p tdSchCurr1.CHROM.sizes -z tdSchCurr1.primary.renamed.fa -t 8 -D software/juicer/
```
## 02 - Population structure <a name="pca"></a>
### Make a PCA plot
```
# Convert vcf to plink .bed format, select autosomes. 
plink2 --vcf FREEZE.FULLFILTER.vcf --chr CHR_1, CHR_2, CHR_3, CHR_4, CHR_5, CHR_6, CHR_7 --make-bed --allow-extra-chr --set-all-var-ids @_# --out autosomes_unfiltered

# Remove variants in strong linkage disequilibrium
plink2 --bfile autosomes_unfiltered --allow-extra-chr --set-all-var-ids @_# --indep-pairwise 50 10 0.15
plink2 --bfile autosomes_unfiltered --allow-extra-chr --set-all-var-ids @_# --extract autosomes_unfiltered.prune.in --out prunedData --make-bed
```
### Principal component analysis
```
plink2 --bfile prunedData --allow-extra-chr --set-all-var-ids @_# --pca
```
### Neighbour-joining phylogeny
```
# Produce a distance matrix
plink2 --bfile prunedData --allow-extra-chr --set-all-var-ids @_# --distance square 1-ibs
paste <( cut -f2 prunedData_tree.mdist.id) prunedData_tree.mdist | cat <(cut -f2 prunedData_tree.mdist.id | tr '\n' '\t' | sed -e '1s/^/\t/g') - > autosomes.mdist
```
### Maximum likelihood phylogeny
```
# Get a list of variants not in LD
cat autosomes_unfiltered.prune.in | tr '_' '\t' > keep.vars.txt

# Subset VCF file to only contain variants not in LD 
bcftools view --threads 6 -T keep.vars.txt -o subset_LD.vcf FREEZE.FULLFILTER.vcf

# Convert to phylip format
vcf2phylip.py -i subset_LD.vcf

# Remove and count invariant sites
ascbias.py -p subset_LD.min4.phy

# Perform phylogenetic inference
iqtree -s out.phy -m MFP --safe -T 2 -B 1000 --seqtype DNA --alrt 1000
```

### Admixture
```
#Fix scaffold names in bim file (ADMIXTURE accepts numerical scaffold names only)
sed -i 's/CHR_//g' prunedData.bim

# Produce random list of seeds
shuf -i 0-10000 | head -10 > seed.list

# Run ADMIXTURE of values of K:1-20, using 10 randomly generated seeds.
# There is no way of renaming admixture output based on seed value, so to avoid overwriting output files for each seed replicate, run in their own directory or batch run each seed one at a time. 
parallel --dry-run "admixture -j2 --seed={1} -B1000 prunedData.bed {2} --cv=10" :::: seed.list ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 > run_admix.sh

# For example:
mkdir 6127
cd 6127
admixture -j2 --seed=6127 -B1000 ../prunedData.bed 1 --cv=10
admixture -j2 --seed=6127 -B1000 ../prunedData.bed 2 --cv=10
# ... up to 20

# Add K values to ADMIXTURE output files
parallel --dry-run "sed -e 's/^/{1} /g' autosomes.{1}.Q" ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
cat *.Q > admixture_all.txt

# Get a table of CV scores, found in the stdout files (in our case *.o files)
cat *.o | grep CV | cut -f2 -d "=" | sed 's/)://g' | tr ' ' '\t' > cv_scores.txt
```
## 03 - Admixture statistics <a name="admix"></a>
### Count heterozygous variants
```
# Make a VCF with biallelic SNPs only
bcftools view --min-ac 2 --types snps --max-ac 2 -o FREEZE.FULLFILTER.biallelic_snps.vcf FREEZE.FULLFILTER.vcf

# Get a list of SNPs per site per sample
bcftools query -f '[%CHROM\t%POS\t%SAMPLE\t%GT\n]' FREEZE.FULLFILTER.biallelic_snps.vcf > all_snps.list
grep -e "0\/1" -e "0|1" all_snps.list > het_snps.list

# Make a bed file (50 kb windows)
samtools faidx tdSchCurr1.primary.fa
bedtools makewindows -g tdSchCurr1.primary.fa.fai -w 50000 > 50kb.bed

# Merge SNP list and bed file
bedtools intersect -wb -a all_snps.list -b 50kb.bed | sed 's/|/\//g' | cut -f4,5,6,7,8,9 | awk '{print "A",$2,$3,$4,$5,$6}' | tr ' ' '\t' | sort -k2,2 -k3,3 -k4,4 | datamash -g2,3,4,5 count 1 | tr '\t' '-' > all.count.bed

bedtools intersect -wb -a het_snps.list -b 50kb.bed | sed 's/|/\//g' | cut -f4,5,6,7,8,9 | awk '{print "A",$2,$3,$4,$5,$6}' | tr ' ' '\t' | sort -k2,2 -k3,3 -k4,4 | datamash -g2,3,4,5 count 1 | tr '\t' '-' > het.count.bed

# Combine all SNP counts and heterozygous SNPs only counts
join <(sort all.count.bed) <(sort het.count.bed) -e0 -a1 -o auto | sed 's/-/ /g' | sed 's/ / /g' > merged.hets.all.txt
```
### *f*3
```
# Make a zarr file
python make_zarr.py

# Per window
python f3_stats.py

# Overall
python f3_sum.py
```
### Patterson's D
```
# Calculate Patterson's D statistics for all trios (where ds.list is a two column file in the format [Sample name]\t[Population])
Dsuite Dtrios FREEZE.FULLFILTER.vcf ds.list
```
### NEWHYBRIDS
```
# Subset VCF to keeep only sites with no missingness
bcftools view -i 'F_MISSING<0.1' FREEZE.FULLFILTER.biallelic_snps.vcf > FREEZE.FULLFILTER.biallelic_snps.nomiss.vcf
bcftools query -f '[%CHROM\t%POS\t%SAMPLE\t%GT\n]' FREEZE.FULLFILTER.biallelic_snps.nomiss.vcf > all_snps.nomiss.list

# Make random subsets of all variants
parallel "cut -f1,2 all_snps.nomiss.list | sort | uniq | shuf > head -200 > subset.{}.txt" ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15

# Make a new subset VCF
parallel "bcftools view -T subset.{}.txt -o subset.{}.vcf FREEZE.FULLFILTER.biallelic_snps.nomiss.vcf" ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15

# Convert to NEWHYBRIDS format
parallel "java -Xmx1024m -Xms512m -jar PGDSpider2-cli.jar -inputfile subset.{}.vcf -inputformat VCF -outputfile subset.{}.newhybrids -outputformat NEWHYBRIDS -spid all.spid" ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15

# Run NEWHYBRIDS for each subset, e.g.:
newhybs --no-gui -d subset.1.newhybrids --num-sweeps 500000 --burn-in 150000 --seeds 15 22
```
## 04 - Windowed admixture <a name="admix2"></a>
### 50 kb window
```
# Prune variants and exclude samples we don't want to analyse
plink2 --vcf FREEZE.FULLFILTER.vcf --chr CHR_1, CHR_2, CHR_3, CHR_4, CHR_5, CHR_6, CHR_7, CHR_Z --make-bed --allow-extra-chr --set-all-var-ids @_# --out autosomes_unfiltered --keep keep.list

# Remove variants in strong linkage disequilibrium
plink2 --bfile autosomes_unfiltered --allow-extra-chr --set-all-var-ids @_# --indep-pairwise 50 10 0.15
plink2 --bfile autosomes_unfiltered --allow-extra-chr --set-all-var-ids @_# --extract autosomes_unfiltered.prune.in --out prunedData --make-bed

# Using the 50 kb windows from 03, make bed files for each 50 kb window in parallel, keeping only those samples you want to analyse (or use as reference populations)
parallel --dry-run --colsep '\t' "plink2 --bfile prunedData --chr {1} --from-bp {2} --to-bp {3} --make-bed --out {1}_{2}_{3} --allow-extra-chr" :::: 50kb.bed

# So for example the first command would be:
plink2 --bfile prunedData --chr 1 --from-bp 0 --to-bp 50000 --make-bed --out 1_0_50000 --allow-extra-chr
ls | grep bed | cut -f1 -d "." > bed.list

# Run ADMIXTURE
parallel "admixture -j1 {}.bed --cv 2" :::: bed.list

# Merge runs
parallel --dry-run "awk '{print \$1,\$2,FILENAME}' {} | sed 's/.2.Q//g' | tr '_' '\t' > {}.X" ::: *.Q
cat *.X > all_admix_50kb.txt
```
### 1 Mb windows
```
# Make windows
bedtools makewindows -g tdSchCurr1.primary.fa.fai -w 1000000 -s 500000 > 1Mb.bed

# Make bed files
parallel --dry-run --colsep '\t' "plink2 --bfile prunedData --chr {1} --from-bp {2} --to-bp {3} --make-bed --out {1}_{2}_{3} --allow-extra-chr" :::: 1Mb.bed
ls | grep bed | cut -f1 -d "." > bed.list

# Run ADMIXTURE
parallel "admixture -j1 {}.bed --supervised --cv 2" :::: bed.list

# Merge runs
parallel --dry-run "awk '{print \$1,\$2,FILENAME}' {} | sed 's/.2.Q//g' | tr '_' '\t' > {}.X" ::: *.Q
cat *.X > all_admix_1Mb.txt
```
## 05 - Mitochondrial genome analysis <a name="mito"></a>
```
# Convert to phylip format
vcf2phylip.py -i MITO.vcf

# Remove and count invariant sites
ascbias.py -p MITO.min4.phy

# Perform phylogenetic inference
iqtree -s out_MITO.phy -m MFP --safe -T 2 -B 1000 --alrt 1000
```
## 06 - Fixed differences <a name="fixed"></a>
```
# Print a table of genotypes per sample
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' FREEZE.FULLFILTER.biallelic_snps.vcf | sed 's/|/\//g' > query.txt

# Identify sites fixed homozygous ref or alt in *S. bovis* and *S. curassoni* samples
cat query.txt | awk '$5==$16' | awk '$16==$25' | awk '$7==$17' | awk '$17==$21' | awk '$5!="./."' | awk '$7!="./."' | awk '$5!=$7' | grep -v Z > fix.difs.txt
```
## 07 - Contamination QC <a name="cont"></a>
```
# Get genotypes and allelic depths for all heterozygous sites
bcftools query -f '[%CHROM\t%POS\t%SAMPLE\t%GT\t%AD\n]' FREEZE.FULLFILTER.biallelic_snps.vcf | sed 's/|/\//g' | awk '$4=="0/1"' | tr ',' '\t' > AD.all.txt
```
## 08 - Repeat Masking <a name="rep"></a>
### Analyse repeat content
```
# Model repeats
RepeatModeler -database tdSchCurr1 -pa 10 -LTRStruct -genomeSampleSizeMax 500000000

# Mask repeats
RepeatMasker -pa 10 -a -s -gff -lib tdSchCurr1-families.fa -dir custom/ -xsmall tdSchCurr1.primary.fa

# Make a bed file
bedtools makewindows -g tdSchCurr1.primary.fa.fai -w 1000000 > 1Mb.bed

# Pick specific types of repeat (e.g. Penelope repeats):
cat tdSchCurr1.primary.fa.out | grep Penelope | tr -s ' ' | sed 's/^ //g' | sed 's/ / /g' | cut -f5,6,7 > Penelope.bed
bedtools intersect -wo -a 1mb.bed -b <( grep -v UNPLACED Penelope.bed | grep -v MITO | sed 's/CHR_//g') | cut -f1,2,3,7 | datamash -g1,2,3 sum 4 | sed 's/$/ PENNY/g' > penny-count.temp

# Merge all repeat types
*-count.temp > rep_features.bed
```
