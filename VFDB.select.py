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
# cutoff_pident = 80
# cutoff_bitscore = 100

inqut_file = pwd + '/' + sampleID + ".VFDB.tsv"

#d_description = pd.read_excel(ref_description,header=1)
#d_head = pd.read_table(ref_head,header=None)

d = pd.read_table(inqut_file,sep="\t",header=None)
d.columns = ['qseqid','sseqid','pident','length','mismatch',
             'gapopen','qstart','qend','sstart','send',
             'evalue','bitscore']

d2 = d.sort_values(by=['qseqid','evalue','bitscore'],ascending=[True,True,True])
## 保留出index
d2['index'] = d2.index
### 排序后 转为字典，默认保留最后一个，
dict = d2.set_index(['qseqid']).T.to_dict()

R = []
for a in dict.values():
    R.append(a['index'])
### bitscore 最高的
d3 = d2[d2.index.isin(R)]

###
d_r = d3

d_r = d_r[['sseqid','pident','bitscore']]
d_r = d_r[d_r['pident'] > cutoff_pident ]
d_r = d_r[d_r['bitscore'] > cutoff_bitscore ]


d_r.loc[:,'pident'] = d_r.loc[:,'pident'].astype(float)
d_r.loc[:,'pident'] = d_r.loc[:,'pident'].map(lambda x:('%.1f')%x)
d_r.loc[:,'bitscore'] = d_r.loc[:,'bitscore'].map(lambda x:('%.1f')%x)


f_out = pwd + "/" + sampleID + ".VFDB_result.txt"
d_r.to_csv(f_out, sep='\t', index=False)
