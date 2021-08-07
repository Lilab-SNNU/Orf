import yaml
import pickle
import re
import argparse
import subprocess
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import time

def classify(tmp_file_location,
             name,
             coverage_counts,
             rcrj, 
             circrnas,
             result_file_location):
             
    tmp = tmp_file_location+'/'+'tmp_file'
    circ = tmp_file_location+'/'+name+'.fa'
    print(circ)
    RCRJ = tmp_file_location+'/'+rcrj+'_RCRJ.fa'
    
    subprocess.call('''awk '$4>{} {}' {} > {}'''.format(coverage_counts,'{print $0}',tmp_file_location+'/'+rcrj, tmp_file_location+'/'+rcrj+'.junction_filter_result'), shell=True)
    
    #getfasta
    subprocess.call('bedtools getfasta -s -fi {} -bed {} -split | fold -w 60 > {}'
    .format(circ,tmp_file_location+'/'+rcrj+'.junction_filter_result',RCRJ),shell=True)

    subprocess.call('mv {} {}'.format(RCRJ, tmp_file_location+'/'+rcrj+'_translated_circ.fa'),shell=True)
    subprocess.call('''sed -i 's/()//g' {}'''.format(tmp_file_location+'/'+rcrj+'_translated_circ.fa'),shell=True)
    final_trans_file = tmp_file_location+'/'+rcrj+'_translated_circ.fa'
    seqs = SeqIO.parse(final_trans_file, 'fasta')
    id_dic = pickle.load(open(tmp_file_location+'/junction_name_dic','rb'))
    translated_circ_id_list = []
    for seq in seqs:
        for i in id_dic:
            start = seq.id.split(':')[-1].split('-')[0]
            end = seq.id.split(':')[-1].split('-')[1]
            if int(start) < int(i) < int(end):
                translated_circ_id_list.append(id_dic[i])
                break

    fastqid = rcrj.split('A')[0]
    han_reads_file = open(tmp_file_location + '/' + fastqid + 'Aligned.sortedByCoord.out.bam.merge_result.RCRJ_result.csv.junction_filter_result')
    reads_list = []
    for line in han_reads_file:
        line = line.replace('\n','')
        line = line.replace('\t','$')
        reads_list.append(line.split('$')[3])
    
            
    total_circ = SeqIO.parse(circrnas, 'fasta')
    final_trans_circ_seq = []
    j = 0
    for i in total_circ:
        if i.id in translated_circ_id_list:
            rec = SeqRecord(i.seq,
                            id = i.id + '|' + reads_list[j],
                            description='')
            final_trans_circ_seq.append(rec)
            j += 1
            print(i)
    fastqid = rcrj.split('A')[0]
    
    final_trans_circ_seq_name = tmp_file_location+'/'+fastqid+'_filter.fa'
    subprocess.call('mkdir -p {}'.format(result_file_location),shell=True)
    SeqIO.write(final_trans_circ_seq, final_trans_circ_seq_name, 'fasta')


    
def main():
    parse = argparse.ArgumentParser(description='This script helps to clean reads and map to genome')
    parse.add_argument('-y', dest="yaml", required=True)
    parse.add_argument('-map-only', dest='map_only', required=False)
    args = parse.parse_args()

    yamlfile = args.yaml
    file = open(yamlfile)
    fileload = yaml.full_load(file)
    name = fileload['genome_name']
    tmp_file_location = fileload['tmp_file_location']
    raw_read = fileload['raw_reads']
    result_file_location = fileload['result_file_location']
    coverage_counts = fileload['coverage_counts']
    circrnas = fileload['circrnas']
    result_file_location = fileload['result_file_location']
    

    rcrj_results = list(filter(lambda x:x[-15:] == 'RCRJ_result.csv', os.listdir(tmp_file_location)))
    for rcrj in rcrj_results:
    	    classify(tmp_file_location, 
                     name, 
                     coverage_counts, 
                     rcrj, 
                     circrnas,
                     result_file_location)
    
    

if __name__ == '__main__':
    main()
