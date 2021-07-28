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
## 02 - Population structure <a name="pca"></a>
### Make a PCA plot
```
```
