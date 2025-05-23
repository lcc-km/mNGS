#!/usr/bin/python
# coding=gbk
import os, sys, json
import re
import argparse

parser = argparse.ArgumentParser(description='runing get fastq QC result')
parser.add_argument('-i', help='the sample ID')
parser.add_argument('-f', help='fastp output .fastp.json')
parser.add_argument('-p', help='path of fastp output')
args = parser.parse_args()
#
sampleID = args.i
file_in = args.f
pwd = args.p

#file_in = pwd + "/" + sampleID + ".fastp.json"
file_out = pwd + "/" + sampleID + ".fastp_QC_result.txt"

fw = open(file_out, "w")

def get_info(json,str1,str2,str3) :
    info =  round(json[str1][str2][str3],3)
    return info

def get_num_info(json,str1,str2,str3) :
    info =  round((json[str1][str2][str3])*100,3)
    return info

with open(file_in,encoding='UTF-8') as f_in:
        d=json.load(f_in)
        raw_reads = get_info(d,"summary","before_filtering","total_reads")
        raw_data_size = get_info(d, "summary", "before_filtering", "total_bases")
        print(raw_data_size)
        raw_q20_rate = get_num_info(d,"summary","before_filtering","q20_rate")
        raw_q30_rate = get_num_info(d,"summary","before_filtering","q30_rate")
        raw_gc_content = get_num_info(d,"summary","before_filtering","gc_content")

        clean_reads = get_info(d,"summary","after_filtering","total_reads")
        clean_data_size = get_info(d, "summary", "after_filtering", "total_bases")
        clean_q20_rate = get_num_info(d,"summary","after_filtering","q20_rate")
        clean_q30_rate = get_num_info(d, "summary", "after_filtering", "q30_rate")
        clean_gc_content = get_num_info(d, "summary", "after_filtering", "gc_content")
        fa_duplication = d["duplication"]["rate"] * 100
        raw_data_size = raw_data_size/(1000*1000*1000)
        clean_data_size = clean_data_size / (1000*1000*1000)
        depth = clean_data_size/3
        passed_filter_reads = round((int(clean_reads) / int(raw_reads)*100),3)
        N_count = d["filtering_result"]["too_many_N_reads"]
        N_rate = round((int(N_count) / int(raw_reads)*100),3)
        out = "%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n" %("sampleID",sampleID,"raw_reads",
                raw_reads,"raw_data_size(G)",raw_data_size,"raw_q20_rate %",raw_q20_rate,"raw_q30_rate %",raw_q30_rate,"raw_gc_content %",raw_gc_content,
                "clean_reads",clean_reads,"clean_data_size(G)",clean_data_size,"depth(X)",depth,"q20_rate %",clean_q20_rate,"q30_rate %",clean_q30_rate,"gc_content %",clean_gc_content,"fastq_duplication %",fa_duplication,
                "passed_filter_reads %",passed_filter_reads,"N_count",N_count, "N_rate %",N_rate)
fw.write(out)
fw.close()


