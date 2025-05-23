#!/usr/bin/python3
import pandas as pd
import numpy as np
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import palettable
from palettable.cartocolors.sequential import DarkMint_4
from itertools import repeat

parser = argparse.ArgumentParser(description='output kraken report list')
parser.add_argument('-i', help='the sample ID')
parser.add_argument('-p', help='path of the sample to be analyzed')
parser.add_argument('-c', help='number of cutoff')
args = parser.parse_args()

sampleID = args.i
pwd = args.p
cut_off = args.c
cut_off = int(cut_off)

# sampleID = 'T7_raw'
# pwd = '/Volumes/mac_ExFAT/微创文件/3.mNGS/seq_data/mNGS/T7'
# cut_off = 5

inqut_file = pwd + "/" + sampleID + ".bracken.mpa-style.report.plot.txt"
classFORplot_in = pwd + "/" + sampleID + ".classFORplot.txt"
d = pd.read_table(inqut_file,sep="\t",header=None)

d.columns = ['k','k2','g','s','count']
d["percent"] = d["count"] / d["count"].sum()
d['RPM'] = (d['count'] * 1000000 ) / d['count'].sum()
d = d[d['RPM'] > cut_off ]

s_Bacteria = "k__Bacteria"
s_Archaea = "k__Archaea"
s_Fungi = "k__Fungi"
s_Viruses = "k__Viruses"

target_list_Bacteria = ['k','g','s','count','percent','RPM']
target_list_Fungi = ['k2','g','s','count','percent','RPM']

def select_data(data, taxonomy, taxonomy_colname, target_list):
    data = data[data[taxonomy_colname] == taxonomy]
    data = data[target_list]
    data.columns = ['k', 'g', 's', 'count', 'percent','RPM']
    # data_pass = data[data['RPM'] > cut_off ]
    data["percent_k"] = data["count"] / data["count"].sum()
    data_pass = data[data['percent_k'] >= cut_off / 1000]
    data_low = data[data['percent_k'] < cut_off / 1000]
    info = taxonomy + " non-significant"
    d_dict = {"k": taxonomy,
              "g": [info],
              's': [info],
              "count": [data_low['count'].sum()],
              "percent": [data_low['percent'].sum()],
              "percent_k": [data_low['percent_k'].sum()]}
    data_tmp = pd.DataFrame(d_dict)
    data2 = pd.concat([data_pass, data_tmp])
    # data2 = data_pass
    return data2

d_Bacteria = select_data(d, s_Bacteria,"k",target_list_Bacteria)
d_Viruses = select_data(d,s_Viruses,"k",target_list_Bacteria)
d_Fungi = select_data(d,s_Fungi,"k2",target_list_Fungi)
d_Archaea = select_data(d,s_Archaea,"k",target_list_Bacteria)
d_r = pd.concat([d_Bacteria,d_Fungi,d_Viruses,d_Archaea])

d_r =d_r [d_r['count'] > 0  ]
#d_r = d_r.dropna()
d_r.loc[:,'RPM'] = d_r.loc[:,'RPM'].map(lambda x:('%.0f')%x)
d_r.loc[:,'percent_k'] = d_r.loc[:,'percent_k'] * 100
d_r.loc[:,'percent_k'] = d_r.loc[:,'percent_k'].map(lambda x:('%.2f')%x)
d_r.loc[:,'percent'] = d_r.loc[:,'percent'] * 100
d_r.loc[:,'percent'] = d_r.loc[:,'percent'].map(lambda x:('%.2f')%x)

f_out = pwd + "/" + sampleID + ".kraken2_bracken_result.txt"
d_r.to_csv(f_out, sep='\t', index=False)

######  plot #######
for i in  list(set(d_r.k)) :
    d_r_k = d_r[d_r.k == i]  #  选择 物种大类
    d_r_k.loc[:,'lable'] = d_r_k.loc[:,'s'].astype(str) +" (" +d_r_k.loc[:,'percent_k'].astype(str) +"%)"  # 给 s__添加 丰都
    d_r_k.loc[:,'percent_k'] = d_r_k.loc[:,'percent_k'].astype(float)
    d_r_k = d_r_k.sort_values(by=['percent_k'],ascending=[False])
    d_r_k = d_r_k.reset_index(drop=True)
    L = len(d_r_k.s) ## 设置使用颜色的数量
    figure_size = (15, 15)  ##  设置图片大小
    plt.figure(figsize=figure_size, dpi=100)
    patches, texts = plt.pie(d_r_k['percent_k'],
        colors=plt.get_cmap(DarkMint_4.mpl_colormap)(np.linspace(0,1,L)), ###  等份取色号
                             )
    plt.legend(patches, list(d_r_k['lable']),#添加图例
          title= (i + ' composition ratio'),
          loc='upper center', bbox_to_anchor=(0.5,0.1),
          ncol=4,#控制图例中按照两列显示，默认为一列显示，
                )
    # plt.show()
    plot_output = pwd + "/plot/" + sampleID + "." + i + '.jpg'
    plt.savefig(plot_output,bbox_inches='tight', pad_inches=0.0, dpi=300)

#### plot K__
d2 = pd.read_table(classFORplot_in,sep="\t",header=None)
d2.columns =  ['percent','count','info','LEVEL','ID','taxonomy']
list_K = ['Bacteria','Fungi','Archaea','Viruses']
d2_K = d2[d2.taxonomy.isin(list_K)]
d2_K.loc[:,'lable'] = d2_K.loc[:,'taxonomy'].astype(str) +" (" +d2_K.loc[:,'percent'].astype(str) +"%)"  # 添加 丰都

patches, texts = plt.pie(d2_K['percent'],
    colors=plt.get_cmap(DarkMint_4.mpl_colormap)(np.linspace(0, 1, 4)) ###  等份取色号
    )
plt.legend(patches, list(d2_K['lable']),  # 添加图例
           title=('composition ratio'),
           loc='upper center', bbox_to_anchor=(0.5, 0.1),
           ncol=2,  # 控制图例中按照两列显示，默认为一列显示，
           )
#plt.show()
plot_output = pwd + "/plot/" + sampleID + ".classFORplot.jpg"
plt.savefig(plot_output, bbox_inches='tight', pad_inches=0.0, dpi=300)