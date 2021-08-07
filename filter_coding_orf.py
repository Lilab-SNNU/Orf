from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import yaml
import argparse
import time
import os

# reverse_complement the '-' strand
def seq_pretreat(genome_name, circrnas_file, tmp_file_path):
     
     han_raw_circ_file = open(circrnas_file)
     
     
     global true_raw_file
     true_raw_file = tmp_file_path +'/'+ genome_name + '_true.fa'
     data_circ = []

     for seq_record in SeqIO.parse(han_raw_circ_file, 'fasta'):
          id_list = seq_record.description.split('|')
          circ_site = id_list[1]
          strand = circ_site[-1]                                      #strand
          if strand == '-':
               seq_true = seq_record.seq.reverse_complement()
               rec_circ_true = SeqRecord(seq_true,
                                            id = seq_record.id,
                                            description = '')
               data_circ.append(rec_circ_true)
          elif strand == '+':
               rec_circ_true = SeqRecord(seq_record.seq,
                                         id = seq_record.id,
                                         description = '')
               data_circ.append(rec_circ_true)
               

     han_true = open(true_raw_file, 'w+')
     SeqIO.write(data_circ, han_true, 'fasta')
     han_true.close()
     han_raw_circ_file.close()

# get all the potential ORF within the single circ and the upstream sequence of these ORF 
def get_1_orf_up(lenth_need, tmp_file_name, tmp_file_path):

     han_circ_file = open(true_raw_file)
     up_records = []
     orf_records = []
     stop_codon = ['TAA', 'TGA', 'TAG']

     for seq_record in SeqIO.parse(han_circ_file, 'fasta'):

          # show the sequence information being searched
          print('\r{}'.format(seq_record.id), end = '')

          seq_double = seq_record.seq * 2
          id_list = seq_record.id.split('|')
          strand = id_list[1][-1]
          
          # the case where the start codon does not cross the junction
          for i in range(len(seq_record.seq)-2):
               if seq_double[i:i+3] == 'ATG':
                    for j in range(i+3, len(seq_double), 3):
                         if j >= len(seq_record.seq)-2 and seq_double[j:j+3] in stop_codon and 9 <= j+3-i <= len(seq_record.seq):

                              new_up_seq = str(seq_double[:i+len(seq_record.seq)][-lenth_need:])
                              new_orf_seq = str(seq_double[i:j+3])

                              # calculate the serial number(not the true site)
                              if len(seq_record.seq)-1 < j:
                                   j = j -len(seq_record.seq)
                                   
                              # map the serial number to the positive strand
                              if strand == '-':
                                   i = len(seq_record.seq) - 1 - i
                                   j = len(seq_record.seq) - 1 - j

                              rec_up = SeqRecord(Seq(new_up_seq),
                                                 id = seq_record.id + '|start{}-stop{}'.format(i+1, j+1) + '|single',
                                                 description = '')
                              up_records.append(rec_up)

                              rec_orf = SeqRecord(Seq(new_orf_seq),
                                                  id = seq_record.id + '|start{}-stop{}'.format(i+1, j+1) + '|single',
                                                  description = '')
                              orf_records.append(rec_orf)
                              break

                         elif seq_double[j:j+3] in stop_codon and j < len(seq_record.seq)-2:
                              break

                         elif seq_double[j:j+3] in stop_codon and j+3-i <9:
                              break

          # the case where the start codon crosses the junction
          for i in [len(seq_record.seq)-2, len(seq_record.seq)-1]:
               if seq_double[i:i+3] == 'ATG':
                    for j in range(i+3, len(seq_double), 3):
                         if seq_double[j:j+3] in stop_codon and 9 <= j+3-i <= len(seq_record.seq):

                              new_up_seq = str(seq_double[:i][-lenth_need:])
                              new_orf_seq = str(seq_double[i:j+3])

                              # calculate the serial number
                              if  len(seq_record.seq)-1 < j:
                                   j = j - len(seq_record.seq)
                              
                              # map the serial number to the positive strand
                              if strand == '-':
                                   i = len(seq_record.seq) - 1 - i
                                   j = len(seq_record.seq) - 1 - j
                              
                              rec_up = SeqRecord(Seq(new_up_seq),
                                                 id = seq_record.id + '|start{}-stop{}'.format(i+1, j+1) + '|single',
                                                 description = '')
                              up_records.append(rec_up)

                              rec_orf = SeqRecord(Seq(new_orf_seq),
                                                  id =seq_record.id + '|start{}-stop{}'.format(i+1, j+1) + '|single',
                                                  description = '')
                              orf_records.append(rec_orf)
                              break

                         elif seq_double[j:j+3] in stop_codon and j+3-i < 9:
                              break

     # record the orf and upstream sequence in the tmp file
     global upstream_raw_file
     upstream_raw_file = tmp_file_path +'/'+ tmp_file_name + '_up.fa'
     han_up = open(upstream_raw_file,'w+')
     SeqIO.write(up_records, han_up, 'fasta')

     global orf_raw_file
     orf_raw_file = tmp_file_path +'/'+ tmp_file_name + '_raw_orf.fa'
     han_raw_orf = open(orf_raw_file, 'w+')
     SeqIO.write(orf_records, han_raw_orf, 'fasta')

     han_up.close()
     han_raw_orf.close()
     han_circ_file.close()

# get all the potential ORF within the double circ and the upstream sequence of these ORF     
def get_2_orf_up(lenth_need, tmp_file_name, tmp_file_path):

     han_circ_file = open(true_raw_file)
     up_records = []
     orf_records = []
     stop_codon = ['TAA', 'TGA', 'TAG']

     for seq_record in SeqIO.parse(han_circ_file, 'fasta'):
          
          # show the sequence information being searched
          print('\r{}'.format(seq_record.id), end = '')

          seq_triple = seq_record.seq * 3
          id_list = seq_record.id.split('|')
          strand = id_list[1][-1]
          
          for i in range(len(seq_record.seq)):
               if seq_triple[i:i+3] == 'ATG':
                    for j in range(i+3, len(seq_triple), 3):
                         if seq_triple[j:j+3] in stop_codon and len(seq_record.seq)*2 >= j+3-i > len(seq_record.seq):

                              new_up_seq = str(seq_triple[:i+len(seq_record.seq)][-lenth_need:])
                              new_orf_seq = str(seq_triple[i:j+3])

                              # calculate the serial number
                              if len(seq_record.seq)-1 < j <= len(seq_record.seq)*2-1:
                                   j = j - len(seq_record.seq)
                              elif len(seq_record.seq)*2-1 < j:
                                   j = j - len(seq_record.seq) * 2

                              # map the serial number to the '+' strand 
                              if strand == '-':
                                   i = len(seq_record.seq) - 1 - i
                                   j = len(seq_record.seq) - 1 - j
                                   
                              rec_up = SeqRecord(Seq(new_up_seq),
                                                 id = seq_record.id + '|start{0}-stop{1}'.format(i+1, j+1) + '|double',
                                                 description ='')
                              up_records.append(rec_up)

                              rec_orf = SeqRecord(Seq(new_orf_seq),
                                                 id = seq_record.id + '|start{0}-stop{1}'.format(i+1, j+1) + '|double',
                                                 description ='')
                              orf_records.append(rec_orf)    
                              break

                         elif seq_triple[j:j+3] in stop_codon and j+3-i < len(seq_record.seq):
                              break

     global upstream_raw_file
     upstream_raw_file = tmp_file_path +'/'+ tmp_file_name + '_up.fa'
     han_up = open(upstream_raw_file,'a+')
     SeqIO.write(up_records, han_up, 'fasta')
     
     global orf_raw_file
     orf_raw_file = tmp_file_path +'/'+ tmp_file_name + '_raw_orf.fa'
     han_raw_orf = open(orf_raw_file, 'a+')
     SeqIO.write(orf_records, han_raw_orf, 'fasta')

     han_up.close()
     han_raw_orf.close()
     han_circ_file.close()

# give a ires score for each ORF
def IRES_score(tmp_file_name, tmp_file_path):

     global IRES_score_file
     IRES_score_file = tmp_file_path +'/'+ tmp_file_name + '_up_IRESscore.result'
     
     os.system('python requiredSoft/IRESfinder_final/IRESfinder.py -f {} -o {}'.format(upstream_raw_file, IRES_score_file))

# filtrate from all the potential ORF fasta by comparing with the IRESscore file
def filter_ires(ires_score, tmp_file_name, tmp_file_path):
     
     global orf_IRES_filter_file
     orf_IRES_filter_file = tmp_file_path + '/' + tmp_file_name + '_up_IRESfilter.fa'
     
     han_raw_orf = open(orf_raw_file)
     han_IRES_score_file = open(IRES_score_file)
     IRESfilter_records = []
     data_IRES_score = []
     dic_orf = {}

     for line in han_IRES_score_file:
          line = line.replace('\n', '')
          line = line.replace('\t', ',')
          data_IRES_score.append(list(line.split(',')))

     for item in data_IRES_score:
          if item[2] != 'Score' and eval(item[2]) >= ires_score:
               dic_orf[item[0]] = item[2][0:4]                                       # Slice directly without rounding, otherwise '1' will appear

     seqs = SeqIO.parse(han_raw_orf, 'fasta')
     while True:
          try:
               seq_record = next(seqs)
               if seq_record.id in dic_orf.keys():
                    rec = SeqRecord(seq_record.seq,
                                    id=seq_record.id + '|' +dic_orf[seq_record.id],
                                    description='')
                    IRESfilter_records.append(rec)
          except StopIteration:
               break

     han_IRES_filter_file = open(orf_IRES_filter_file, 'w+')
     SeqIO.write(IRESfilter_records, han_IRES_filter_file, 'fasta')
     han_raw_orf.close()
     han_IRES_score_file.close()
     han_IRES_filter_file.close()

#give an orf score for each ORF
def ORF_score(tmp_file_name, tmp_file_path):

     global ORF_score_file
     ORF_score_file = tmp_file_path +'/'+ tmp_file_name + '_ORFscore'
     
     os.system('python requiredSoft/CPC_final/bin/CPC2.py -i {} -o {}'.format(orf_IRES_filter_file, ORF_score_file))

#filtrate from all the potential ORF fasta by comparing with the score file
def filter_orf(orf_score, final_name, final_file_path):
     
     global orf_true_final_file
     orf_true_final_file = final_file_path +'/'+ final_name + '_true.fa'
     
     han_ORF_score_file = open(ORF_score_file + '.txt')
     han_IRES_filter_file = open(orf_IRES_filter_file)
     data_ORF_score = []
     dic_score_orf = {}
     ORFfilter_records = []

     for line in han_ORF_score_file:
          line = line.replace('\t', ',')
          line = line.replace('\n', '')
          data_ORF_score.append(line.split(','))

     for item in data_ORF_score:
          if item[6] != 'coding_probability' and eval(item[6]) >= orf_score:
               dic_score_orf[item[0]] = item[6][0:4]

     seqs = SeqIO.parse(han_IRES_filter_file, 'fasta')
     while True:
          try:
               seq_record = next(seqs)
               if seq_record.id in dic_score_orf.keys():
                    rec = SeqRecord(seq_record.seq,
                                    id=seq_record.id + '|' + dic_score_orf[seq_record.id],
                                    description='')
                    ORFfilter_records.append(rec)
          except StopIteration:
               break

     han_orf_true_final_file = open(orf_true_final_file, 'w+')
     SeqIO.write(ORFfilter_records, han_orf_true_final_file, 'fasta')
     han_ORF_score_file.close()
     han_IRES_filter_file.close()
     han_orf_true_final_file.close()

# Make a virtual sequence and reverse_complement the '-' strand
def reverse_complment(final_name, final_file_path):
     
     han_true = open(orf_true_final_file)
     virtual_file = final_file_path +'/'+ final_name + '_virtual.fa'
     data_circ = []

     for seq_record in SeqIO.parse(han_true, 'fasta'):
          id_list = seq_record.description.split('|')
          circ_site = id_list[1]
          strand = circ_site[-1]                                      #strand
          if strand == '-':
               seq_virtual = seq_record.seq.reverse_complement()
               rec_circ_virtual = SeqRecord(seq_virtual,
                                            id = seq_record.id,
                                            description = '')
               data_circ.append(rec_circ_virtual)
          elif strand == '+':
               rec_circ_true = SeqRecord(seq_record.seq,
                                         id = seq_record.id,
                                         description = '')
               data_circ.append(rec_circ_true)

     han_virtual = open(virtual_file,'w+')
     SeqIO.write(data_circ, han_virtual, 'fasta')
     han_true.close()
     han_virtual.close()

# Make the amino acid sequence file according to the true ORF file
def translate_orf(final_name, final_file_path):

     true_amino_acid_file = final_file_path +'/'+ final_name + '_true_translated.fa'
     han_true = open(orf_true_final_file)
     true_circ = []
     for seq_record in SeqIO.parse(han_true, 'fasta'):
          rna_circ = seq_record.seq.transcribe()
          amino_acid_seq = rna_circ.translate()
          rec = SeqRecord(amino_acid_seq,
                          id=seq_record.id,
                          description='')
          true_circ.append(rec)
     han_amino_acid = open(true_amino_acid_file, 'w+')
     SeqIO.write(true_circ, han_amino_acid, 'fasta')
     han_true.close()
     han_amino_acid.close()

def main():
     parse = argparse.ArgumentParser(description='This script helps to filter coding orf')
     parse.add_argument('-y', '--yaml', required=True, help='please input the yamlfile')
     args = parse.parse_args()
     
     yamlfile = args.yaml
     con_file = open(yamlfile)
     fileload = yaml.full_load(con_file)
     
     lenth_need = fileload['lenth_need']
     raw_reads = fileload['raw_reads']
     ires_score = fileload['ires_score']
     orf_score = fileload['orf_score']
     tmp_file_path = fileload['tmp_file_location']
     final_file_path = fileload['result_file_location']

     for item in raw_reads:
          
          genome_name = item.split('/')[-1][:-4]
          circrnas_file = fileload['tmp_file_location']+'/'+genome_name+'_filter.fa'

          # tmp file name and path
          tmp_file_name = genome_name + '_orf'
     
          # final file name and path
          final_name = genome_name + '_orf_filter_result'

     
          begin = time.perf_counter()
          #lenth_need, genome_name, circrnas_file, ires_score, orf_score, tmp_file_name, tmp_file_path, final_name, final_file_path = file_load()

          time0 = time.perf_counter()
          seq_pretreat(genome_name, circrnas_file, tmp_file_path)
          time_add0 = time.perf_counter()-time0
          print('seq_pretreat is finished! Spending time is {:.2f}s'.format(time_add0))

          time1 = time.perf_counter()
          get_1_orf_up(lenth_need, tmp_file_name, tmp_file_path)
          time_add1 = time.perf_counter()-time1
          print('\nget_1_orf_up is finished! Spending time is {:.2f}s'.format(time_add1))

          #time2 = time.perf_counter()
          #get_2_orf_up(lenth_need, tmp_file_name, tmp_file_path)
          #time_add2 = time.perf_counter()-time2
          #print('\nget_2_orf_up time is:{:.2f}'.format(time_add2))

          time3 = time.perf_counter()
          IRES_score(tmp_file_name, tmp_file_path)
          time_add3 = time.perf_counter()-time3
          print('IRES_score is finished! Spending time is {:.2f}s'.format(time_add3))

          time4 = time.perf_counter()
          filter_ires(ires_score, tmp_file_name, tmp_file_path)
          time_add4 = time.perf_counter()-time4
          print('filter_ires is finished! Spending time is {:.2f}s'.format(time_add4))

          time5 = time.perf_counter()
          ORF_score(tmp_file_name, tmp_file_path)
          time_add5 = time.perf_counter()-time5
          print('ORF_score is finished! Spending time is {:.2f}s'.format(time_add5))

          time6 = time.perf_counter()
          filter_orf(orf_score, final_name, final_file_path)
          time_add6 = time.perf_counter()-time6
          print('filter_orf is finished! Spending time is {:.2f}s'.format(time_add6))

          time7 = time.perf_counter()
          reverse_complment(final_name, final_file_path)
          time_add7 = time.perf_counter()-time7
          print('reverse_complment is finished! Spending time is {:.2f}s'.format(time_add7))

          time8 = time.perf_counter()
          translate_orf(final_name, final_file_path)
          time_add8 = time.perf_counter()-time8
          print('translate_orf is finished! Spending time is {:.2f}s'.format(time_add8))

          print('the whole filter_coding_orf is finished! Spending time is {:.2f}s'.format(int(time.perf_counter() - begin)))
          print('End')

if __name__ == '__main__':
     main()
     

