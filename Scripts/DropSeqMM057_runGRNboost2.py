import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import genie3, grnboost2
from distributed import Client, LocalCluster
import time
import sys

zd = '/ddn1/vol1/staging/leuven/stg_00002/lcb/zkalender/scRNA_seq_melanoma/analysis_40k/20.analysis_with_sub_matrices/'
subset = "53.MM057_DropSeq"

indir = zd + '{0}/10.TextData/'.format(subset)
saveout = sys.stdout
logs = open(indir + 'arboreto.log', 'w')
sys.stdout = logs

local_cluster = LocalCluster(n_workers=20, threads_per_worker=1, memory_limit=8e9)
custom_client = Client(local_cluster)
print(custom_client)

#### get expression matrix 
ex_path = indir + 'log_CPM_matrix.tsv'
print('reading expression data: ' + ex_path)
ex_matrix = pd.read_csv(ex_path, index_col=0, sep='\t', skiprows=1, header = None).T
#### get tfs
tf_path = '/ddn1/vol1/staging/leuven/stg_00002/lcb/kspan/analyses/ThreeLines10xSCENIC2/hg19_allTFs.lst'
tf_names = load_tf_names(tf_path)
tf_names = list(set(tf_names).intersection(ex_matrix.columns))
print(len(tf_names))

#run grnboost2
outfile = indir + 'grnboost2.tsv'
print('grnboost2 results will be printed to: ' + outfile)
start_time = time.time()
network = grnboost2(expression_data=ex_matrix,
	tf_names=tf_names,
	client_or_address=custom_client,
	verbose = True)
print(time.time() - start_time, "seconds")
print(network.head())
network.to_csv(outfile , sep='\t', index=False, header=False)
sys.stdout = saveout
logs.close()

exit()

# packages in environment at /data/leuven/306/miniconda3/envs/arboreto: {
# Name                    Version                   Build  Channel
#arboreto                  0.1.5                      py_0    bioconda
#blas                      1.0                         mkl  
#bokeh                     0.13.0                   py36_0  
#ca-certificates           2018.03.07                    0  
#certifi                   2018.8.24                py36_1  
#click                     6.7                      py36_0  
#cloudpickle               0.5.5                    py36_0  
#cytoolz                   0.9.0.1          py36h14c3975_1  
#dask                      0.19.1                   py36_0  
#dask-core                 0.19.1                   py36_0  
#distributed               1.23.1                   py36_0  
#heapdict                  1.0.0                    py36_2  
#intel-openmp              2019.0                      118  
#jinja2                    2.10                     py36_0  
#libedit                   3.1.20170329         h6b74fdf_2  
#libffi                    3.2.1                hd88cf55_4  
#libgcc-ng                 8.2.0                hdf63c60_1  
#libgfortran-ng            7.3.0                hdf63c60_0  
#libstdcxx-ng              8.2.0                hdf63c60_1  
#locket                    0.2.0                    py36_1  
#markupsafe                1.0              py36h14c3975_1  
#mkl                       2019.0                      118  
#mkl_fft                   1.0.4            py36h4414c95_1  
#mkl_random                1.0.1            py36h4414c95_1  
#msgpack-python            0.5.6            py36h6bb024c_1  
#ncurses                   6.1                  hf484d3e_0  
#numpy                     1.14.5           py36h1b885b7_4  
#numpy-base                1.14.5           py36hdbf6ddf_4  
#openssl                   1.0.2p               h14c3975_0  
#packaging                 17.1                     py36_0  
#pandas                    0.23.4           py36h04863e7_0  
#partd                     0.3.8                    py36_0  
#pip                       10.0.1                   py36_0  
#psutil                    5.4.7            py36h14c3975_0  
#pyparsing                 2.2.0                    py36_1  
#python                    3.6.6                hc3d631a_0  
#python-dateutil           2.7.3                    py36_0  
#pytz                      2018.5                   py36_0  
#pyyaml                    3.13             py36h14c3975_0  
#readline                  7.0                  h7b6447c_5  
#scikit-learn              0.19.1           py36hedc7406_0  
#scipy                     1.1.0            py36hd20e5f9_0  
#setuptools                40.2.0                   py36_0  
#six                       1.11.0                   py36_1  
#sortedcontainers          2.0.5                    py36_0  
#sqlite                    3.24.0               h84994c4_0  
#tblib                     1.3.2                    py36_0  
#tk                        8.6.8                hbc83047_0  
#toolz                     0.9.0                    py36_0  
#tornado                   4.5.2            py36h1283b2a_0  
#wheel                     0.31.1                   py36_0  
#xz                        5.2.4                h14c3975_4  
#yaml                      0.1.7                had09818_2  
#zict                      0.1.3                    py36_0  
#zlib                      1.2.11               ha838bed_2  
#}
