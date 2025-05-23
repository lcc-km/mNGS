# 不是很复杂的mNGS有参分析流程
## 主要步骤：

1.下机fq文件的质控，用fastp完成

2.比对到人类基因组，选择未比对上的序列重新输出为fq文件 （如果是人类💩样本，其实感觉人类基因序列几乎没有，跳过还可以省点时间。其他样本看情况，去掉人源数据能减少数据量，加快后续流程的分析速度）

3.kraken2 快速获得物种分类信息

4.humann获得乱七八糟的各种注释，同时也有分类信息，不过和 kraken2 略有出入

5.megahit 组装一下序列，编成长点的 Contig / Scaffold

6.编好的Contig / Scaffold 用 diamond 比对到 VFDB(毒力因子)，CARD(耐药性)，args_oap(耐药性)的数据库，获得上述的注释。

## 确定文件名称（如下图）这个适配自己的文件名
<img width="143" alt="image" src="https://github.com/user-attachments/assets/1064ae53-bb0d-4a2d-a84c-deb08f8859da" />
0509T1作为程序运行的ID名，改变工作路径后，运行命令：

## 转化路径至 fq文件所在路径
`cd /your_file_path`
需确保 fq 文件名称为 `"$sampleID"_1.fastq.gz` 以及 `"$sampleID"_2.fastq.gz` 的格式 #SRA拆分后的默认格式

## 运行
```bash mNGS_reference.sh --sampleID 0509T1 --input_folder `pwd` --threads 16```

## 分类的文本结果 在这个文件夹中：
`$input_folder/"$sampleID"/3.annotation/taxonomy `

`$sampleID.bracken.mpa-style.report.txt
$sampleID.kraken2_bracken_result.txt`

## 还有一些图片结果在这个文件夹中：
`$input_folder/"$sampleID"/3.annotation/taxonomy/plot`
![image](https://github.com/user-attachments/assets/d4417e5b-b1b3-41a9-b2a6-421d0ced2708)
![image](https://github.com/user-attachments/assets/892e903f-5d5d-42d5-8191-a05813b25a5f)
![image](https://github.com/user-attachments/assets/86925952-c365-44c9-9e55-b4848b80d608)



