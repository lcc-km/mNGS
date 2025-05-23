#!/usr/bin/python3
# coding=gbk
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='output kraken report list')
parser.add_argument('-i', help='the sample ID')
parser.add_argument('-p', help='path of the sample to be analyzed')
parser.add_argument('-c', help='cutoff of pident')
parser.add_argument('-b', help='cutoff of bitscore')
args = parser.parse_args()

sampleID = args.i
pwd = args.p
cutoff_pident = args.c
cutoff_pident = int(cutoff_pident)
cutoff_bitscore = args.b
cutoff_bitscore = int(cutoff_bitscore)

# sampleID = 'T7_raw'
# pwd = '/Volumes/mac_ExFAT/微创文件/3.mNGS/seq_data/T7'
# cutoff_pident  = 98
# cutoff_bitscore = 500

inqut_file = pwd + '/' + sampleID + ".blastout.ARO.txt"

d = pd.read_table(inqut_file,sep="\t",header=None,converters={6:str})
d.columns = ['qseqid','sseqid','pident','length','info1',
             'info2','evalue','bitscore']

d2 = d.sort_values(by=['qseqid','evalue','bitscore'],ascending=[True,True,True])
## 保留出index
d2['index'] = d2.index
### 排序后 转为字典，默认保留最后一个，
dict = d2.set_index(['qseqid']).T.to_dict()
R = []
for a in dict.values():
    R.append(a['index'])
### bitscore 最高的
d_r = d2[d2.index.isin(R)]

d_r = d_r[d_r['pident'] > cutoff_pident ]
d_r = d_r[d_r['bitscore'] > cutoff_bitscore ]
d_r = d_r[['sseqid','pident','bitscore']]
d_r.loc[:,'pident'] = d_r.loc[:,'pident'].astype(float)
d_r.loc[:,'pident'] = d_r.loc[:,'pident'].map(lambda x:('%.1f')%x)
d_r.loc[:,'bitscore'] = d_r.loc[:,'bitscore'].map(lambda x:('%.0f')%x)


# d_r2 = d_r['sseqid'].str.split('|',expand=True).loc[:,(2,3)]
# d_r2 = d_r2.reset_index(drop=True)
# d_r2.columns = ['id','name']

f_out = pwd + "/" + sampleID + ".args_ARO_result.txt"
d_r.to_csv(f_out, sep='\t', index=False)

