#!/bin/bash
#  nohup sudo bash /mnt/vol1/lucc/work/script/04.mNGS_16s/mNGS.sh --sampleID $sampleID --input_folder `pwd` --threads 32

# bash /mnt/vol1/lucc/work/script/04.mNGS_16s/mNGS.sh --sampleID SRR8849290 --input_folder `pwd` --threads 64 

sampleID=$2
input_folder=$4
threads=$6

mkdir -p $input_folder/"$sampleID"/1.fastq_QC
mkdir -p $input_folder/"$sampleID"/2.mapping
mkdir -p $input_folder/"$sampleID"/2.mapping/tmp
# mkdir -p $input_folder/"$sampleID"/2.mapping/megahit  ### megahit 自建文件夹
mkdir -p $input_folder/"$sampleID"/3.annotation
mkdir -p $input_folder/"$sampleID"/3.annotation/taxonomy
mkdir -p $input_folder/"$sampleID"/3.annotation/taxonomy/plot
mkdir -p $input_folder/"$sampleID"/3.annotation/humann
mkdir -p $input_folder/"$sampleID"/4.result

fastq_folder=$input_folder/"$sampleID"/1.fastq_QC
mapping_folder=$input_folder/"$sampleID"/2.mapping
tmp_folder=$input_folder/"$sampleID"/2.mapping/tmp
taxonomy_folder=$input_folder/"$sampleID"/3.annotation/taxonomy
plot_folder=$input_folder/"$sampleID"/3.annotation/taxonomy/plot
humann_folder=$input_folder/"$sampleID"/3.annotation/humann
megahit_folder=$input_folder/"$sampleID"/2.mapping/megahit
# args_oap_folder=$input_folder/"$sampleID"/3.annotation/args_oap
# args_oap_out_folder=$input_folder/"$sampleID"/3.annotation/args_oap_out
VFDB_foldeer=$input_folder/"$sampleID"/3.annotation
result_folder=$input_folder/"$sampleID"/4.result

fastp=/toolkit/fastp
bowtie2=/toolkit/bowtie2/bowtie2
samtools=/toolkit/samtools/bin/samtools
sambamba=/toolkit/sambamba
bedtools=/toolkit/bedtools2/bin/bedtools
megahit=/toolkit/megahit/build/megahit
quast=/toolkit/quast-master/quast.py
script=/mnt/vol1/lucc/work/script/other_20241016/04.mNGS_16s/
bowtie2_index=/mnt/vol1/database/ucsc/ref/hg38/bowtie2/hg38.fa
kraken2_db=/mnt/vol1/database/mNGS/kraken2/k2_pluspf
# args_oap_structure1=/mnt/vol1/database/mNGS/args_oap/1.SARG_v3.2_20220917_Full_structure.txt
VFDB_setA=/mnt/vol1/database/mNGS/VFDB/VFDB_setA_pro_db.dmnd
CARD_db=/mnt/vol1/database/mNGS/card/card-data_v4.0.0/card.db.dmnd
args_oap_db=/mnt/vol1/database/mNGS/args_oap/1.SARG_v3.2_20220917_Full_database.db.dmnd

# #运行fastp
# $fastp -w 16 -p -I $input_folder/"$sampleID"_1.fq.gz -O $fastq_folder/$sampleID.trim_R1.fq.gz -i $input_folder/"$sampleID"_2.fq.gz -o $fastq_folder/$sampleID.trim_R2.fq.gz  --json $fastq_folder/$sampleID.fastp.json --html $fastq_folder/$sampleID.html
# # ##运行fastqc
# python /mnt/vol1/lucc/work/script/other_20241016/get_fastq_QC_liunx.py $sampleID $fastq_folder

# # # 比对到人类基因组,获取未比对序列
# $bowtie2 -p $threads -x $bowtie2_index -1 $fastq_folder/$sampleID.trim_R1.fq.gz -2 $fastq_folder/$sampleID.trim_R2.fq.gz -S $mapping_folder/$sampleID.sam --un-conc-gz $mapping_folder/$sampleID.unmap.fq.gz


################################################################ 有参 ################################################################ ./install_kraken2.sh /toolkit/kraken2-2.1.5
################    分类学的定量   ##########################
############ kraken2
# ### 1.kraken2
# echo $sampleID
# /toolkit/kraken2/kraken2 --db $kraken2_db --threads $threads --report $taxonomy_folder/$sampleID.kraken2.report.tsv --output $taxonomy_folder/$sampleID.kraken2.output.tsv --use-names --paired $mapping_folder/$sampleID.unmap.fq.*.gz
# #2 使用Braken校正
# /toolkit/Bracken/bracken -d $kraken2_db -t $threads -r 250 -l S -i $taxonomy_folder/$sampleID.kraken2.report.tsv -o $taxonomy_folder/$sampleID.bracken.out -w $taxonomy_folder/$sampleID.bracken.report.txt
# # 格式转换成--use-mpa-style格式
# /toolkit/kraken2/KrakenTools/kreport2mpa.py -r $taxonomy_folder/$sampleID.bracken.report.txt -o $taxonomy_folder/$sampleID.bracken.mpa-style.report.txt
# #### our report
# /toolkit/microbiome_helper/metaphlan_to_stamp.pl $taxonomy_folder/$sampleID.bracken.mpa-style.report.txt > $taxonomy_folder/$sampleID.bracken.mpa-style.report.spf

# grep "s__" $taxonomy_folder/$sampleID.bracken.mpa-style.report.spf | awk '{print  $1,$2,$(NF-2),$(NF-1),$(NF)}' OFS="\t" > $taxonomy_folder/$sampleID.bracken.mpa-style.report.plot.txt

# grep -E "root|Bacteria|Archaea|Fungi|Viruses" $taxonomy_folder/$sampleID.bracken.report.txt | grep -vE "D1|D2" | sed 's/ //g' > $taxonomy_folder/$sampleID.classFORplot.txt
# python3 $script/kraken2_count.py -i $sampleID -p $taxonomy_folder -c 5

# echo " "
# echo "mNGS process completed successfully!" 
# cat $taxonomy_folder/$sampleID.bracken.mpa-style.report.plot.txt
# echo "################" 
# cat $taxonomy_folder/$sampleID.kraken2_bracken_result.txt
# ## krona
# /toolkit/kraken2/KrakenTools/kreport2krona.py -r $taxonomy_folder/$sampleID.bracken.report.txt -o $taxonomy_folder/$sampleID.bracken.krona.txt


########################   功能组成的定量    ##################################
# cat $mapping_folder/$sampleID.unmap.fq.*.gz > $mapping_folder/$sampleID.unmap.fq.gz

# echo "humann"  
# humann --threads $threads --input $mapping_folder/$sampleID.unmap.fq.gz --output $humann_folder --memory-use maximum --nucleotide-database /mnt/vol1/database/mNGS/HUMAnN/chocophlan.v201901_v31 --protein-database /mnt/vol1/database/mNGS/HUMAnN/uniref90_annotated_v201901b/

# # ### 输出 taxonomy
# merge_metaphlan_tables.py $humann_folder/$sampleID.unmap_humann_temp/$sampleID.unmap_metaphlan_bugs_list.tsv | sed 's/_metaphlan_bugs_list//g' > $humann_folder/$sampleID.taxonomy.tsv
    
# # ### 标准化相对丰度
# # echo "标准化相对丰度"
# bash $script/humann_annotation.sh --sampleID $sampleID --humann_folder $humann_folder
# ### 合并注释结果
# cat $humann_folder/*cpm_rename.tsv | grep -Ev "#|UNMAPPED|UNGROUPED|\|" > $humann_folder/"$sampleID"_merge_ann.txt

######################  组装 ############################
$megahit -1 $mapping_folder/$sampleID.unmap.fq.1.gz -2 $mapping_folder/$sampleID.unmap.fq.2.gz -m 0.5 -t $threads -o $megahit_folder

#####  VFDB毒力注释   
#/toolkit/diamond makedb --in /mNGS/VFDB/VFDB_setA_pro.fas --db /mnt/vol1/database/mNGS/VFDB/VFDB_setA_pro_db
/toolkit/diamond blastx --threads $threads --db $VFDB_setA --query $megahit_folder/final.contigs.fa --out $VFDB_foldeer/$sampleID.VFDB.tsv
## sclect
# python3 $script/VFDB.select.py -i $sampleID -p $VFDB_foldeer -c 98 -b 500

#####  CARD 耐药性注释注释   
#/toolkit/diamond makedb --in /mnt/vol1/database/mNGS/card/card-data_v4.0.0/card.faa --db /mnt/vol1/database/mNGS/card/card-data_v4.0.0/card.db
/toolkit/diamond blastx --threads $threads --db $CARD_db --query $megahit_folder/final.contigs.fa --out $VFDB_foldeer/$sampleID.CARD.tsv

#####  args_oap 耐药性注释注释   
# #/toolkit/diamond makedb --in /mnt/vol1/database/mNGS/args_oap/1.SARG_v3.2_20220917_Full_database.fasta --db  /mnt/vol1/database/mNGS/args_oap/1.SARG_v3.2_20220917_Full_database.db
/toolkit/diamond blastx --threads $threads --db $args_oap_db --query $megahit_folder/final.contigs.fa --out $VFDB_foldeer/$sampleID.args_oap.tsv



# cp $plot_folder/*jpg $result_folder
# cp $VFDB_foldeer/$sampleID.VFDB_result.txt $result_folder
# cp $args_oap_out_folder/$sampleID.args_ARO_result.txt $result_folder





# cp $megahit_folder/final.contigs.fa $args_oap_folder/final.contigs.fa
# args_oap stage_one -i $args_oap_folder -o $args_oap_out_folder -f fa -t $threads --database $args_oap_db
# args_oap stage_two -i $args_oap_out_folder -t $threads --database $args_oap_db --structure1 $args_oap_structure1
# ## sclect
# grep "ARO" $args_oap_structure1/blastout.txt > $args_oap_structure1/$sampleID.blastout.ARO.txt
# python3 $script/args_oap_select.py -i $sampleID -p $args_oap_out_folder -c 98 -b 500