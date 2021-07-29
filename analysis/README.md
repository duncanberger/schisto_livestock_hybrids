# Analyses

## Table of contents
1. [Assembly QC](#AQC)
2. [Population Structure](#pca)
3. [Admixture statistics](#admix)
4. [Mitochondrial genome analysis](#mito)

```
Checklist:
2. PCA, ADMIXTURE, tree
3. ADmixture stats
4. f3, PAT D, Junction counting, NEWHYBRIDS
5. MITO calling
```
## 01 - Assembly QC <a name="AQC"></a>
### Run BUSCO
```
run_BUSCO.py -i tdSchCurr1.chrom.fa -o busco_euk_09 -l lineage_datasets/eukaryota_odb9/ -m genome -c 8 --long -t euk_09
```
### Chromosomal synteny
```
# Run PROMER
promer --mum -p promer Smansoni_v7.fa.masked.masked tdSchCurr1.chrom.fa

# Circos plot
circos -conf circos.conf -param image/radius=850p
```
### Hi-C plots
```
# Made indices 
bwa index tdSchCurr1.primary.renamed.fa
python generate_site_positions.py Arima tdSchCurr1 tdSchCurr1.primary.renamed.fa
awk '{print $1,$NF}' tdSchCurr1_Arima.txt > tdSchCurr1.CHROM.sizes

# Run Juicer
juicer.sh -g tdSchCurr1 -s Arima -y tdSchCurr1_Arima.txt -p tdSchCurr1.CHROM.sizes -z tdSchCurr1.primary.renamed.fa -t 8 -D software/juicer/
------------------------------------
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
