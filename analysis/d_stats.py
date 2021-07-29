import random
random.seed(42)
import time
import numpy as np
np.random.seed(42)
import bcolz
import pandas as pd
import allel; print('scikit-allel', allel.__version__)
import zarr
zarr.__version__

callset_fn = 'autosomes.zarr'
callset = zarr.open_group(callset_fn, mode='r')
metadata_fn = 'metadata.txt'

df = pd.read_csv(metadata_fn, sep='\t')
g = allel.GenotypeChunkedArray(callset['calldata/GT'])
ac=g.count_alleles()
biallelic = (ac.max_allele() == 1)
ac_biallelic = g.subset(biallelic)

SCUR=ac_biallelic.count_alleles(subpop=[2,12,16])
SBOV=ac_biallelic.count_alleles(subpop=[0,11,20])
SHAE=ac_biallelic.count_alleles(subpop=[4,5,6,7,8])

for i, row in df.iterrows():
	globals()[str(row['sample'])+"_ac"]=ac_biallelic.count_alleles(subpop=[row['index']])

for i, row in  df.iterrows():
	zz=allel.moving_patterson_d(SCUR, SBOV, (globals()[str(row['sample'])+"_ac"]), SHAE, 300, start=0, stop=None, step=None)
	df = pd.DataFrame(data=zz.flatten())
	df['sample']=str(row['sample'])
	df.to_csv('d_values.csv', mode='a', header=False)
