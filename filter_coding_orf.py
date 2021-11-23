from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import yaml
import argparse
import time
import os

# reverse_complement the '-' strand
def seq_pretreat(circrna_gtf, genome_name, circrnas_file, tmp_file_path):
     
     han_raw_circ_file = open(circrnas_file)
     han_circ_gtf = open(circrna_gtf)

     dic_strand = {}
     for line in han_circ_gtf:
          line = line.replace('\n', '')
          line_list = line.split('\t')
          dic_strand[line_list[0]] = line_list[1]
     
     global true_raw_file
     true_raw_file = tmp_file_path +'/'+ genome_name + '_true.fa'
     data_circ = []

     for seq_record in SeqIO.parse(han_raw_circ_file, 'fasta'):
          circ_name = seq_record.id.split('$')[0]
          strand = dic_strand[circ_name]
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
     han_circ_gtf.close()

# get all the potential ORF within the single circ and the upstream sequence of these ORF 
def get_1_orf_up(circrna_gtf, lenth_need, tmp_file_name, tmp_file_path):

     han_circ_gtf = open(circrna_gtf)
     han_circ_file = open(true_raw_file)
     up_records = []
     orf_records = []
     stop_codon = ['TAA', 'TGA', 'TAG']
     start_codon = ['ATC','ATA','ATT','CTG','GTG','TTG','ATG','AAG','AGG']

     dic_strand = {}
     for line in han_circ_gtf:
          line = line.replace('\n', '')
          line_list = line.split('\t')
          dic_strand[line_list[0]] = line_list[1]

     for seq_record in SeqIO.parse(han_circ_file, 'fasta'):

          # show the sequence information being searched
          #print('\r{}'.format(seq_record.id), end = '')

          seq_double = seq_record.seq * 2
          circ_name = seq_record.id.split('$')[0]
          strand = dic_strand[circ_name]
          
          # the case where the start codon does not cross the junction
          for i in range(len(seq_record.seq)-2):
               #if seq_double[i:i+3] == 'ATG':
               if seq_double[i:i+3] in start_codon:
                    for j in range(i+3, len(seq_double), 3):
                         if j >= len(seq_record.seq)-2 and seq_double[j:j+3] in stop_codon and 60 <= j+3-i <= len(seq_record.seq):

                              new_up_seq = str(seq_double[:i+len(seq_record.seq)][-lenth_need:])
                              new_orf_seq = str(seq_double[i:j+3])

                              # calculate the serial number(not the true site)
                              if len(seq_record.seq)-1 < j:
                                   j = j -len(seq_record.seq)
                                   
                              # map the serial number to the positive strand
                              if strand == '-':
                                   i = len(seq_record.seq) - 1 - i
                                   j = len(seq_record.seq) - 1 - j

                              #rec_up = SeqRecord(Seq(new_up_seq),
                              #                   id = seq_record.id + '|start{}-stop{}'.format(i+1, j+1) + '|single',
                              #                   description = '')

                              rec_up = SeqRecord(Seq(new_up_seq),
                                                 id = seq_record.id + '|start{}-stop{}'.format(i+1, j+1),
                                                 description = '')
                              up_records.append(rec_up)

                              #rec_orf = SeqRecord(Seq(new_orf_seq),
                              #                    id = seq_record.id + '|start{}-stop{}'.format(i+1, j+1) + '|single',
                              #                    description = '')
                              rec_orf = SeqRecord(Seq(new_orf_seq),
                                                  id = seq_record.id + '|start{}-stop{}'.format(i+1, j+1),
                                                  description = '')
                              orf_records.append(rec_orf)
                              break

                         elif seq_double[j:j+3] in stop_codon and j < len(seq_record.seq)-2:
                              break

                         elif seq_double[j:j+3] in stop_codon and j+3-i <9:
                              break

          # the case where the start codon crosses the junction
          for i in [len(seq_record.seq)-2, len(seq_record.seq)-1]:
               #if seq_double[i:i+3] == 'ATG':
               if seq_double[i:i+3] in start_codon:
                    for j in range(i+3, len(seq_double), 3):
                         if seq_double[j:j+3] in stop_codon and 60 <= j+3-i <= len(seq_record.seq):

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
                                                 id = seq_record.id + '|start{}-stop{}'.format(i+1, j+1),
                                                 description = '')
                              up_records.append(rec_up)

                              rec_orf = SeqRecord(Seq(new_orf_seq),
                                                  id =seq_record.id + '|start{}-stop{}'.format(i+1, j+1),
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
     han_circ_gtf.close()
     han_circ_file.close()

# give a ires score for the upstream sequence of each ORF
def IRES_score(tmp_file_name, tmp_file_path):

     global IRES_score_file
     IRES_score_file = tmp_file_path +'/'+ tmp_file_name + '_up_IRESscore.result'
     
     os.system('python requiredSoft/IRESfinder_final/IRESfinder.py -f {} -o {}'.format(upstream_raw_file, IRES_score_file))

# give a m6a score for the upstream sequence of each ORF
def M6A_score(tmp_file_name, tmp_file_path, pos_fa, neg_fa):
     
     global M6A_score_file
     M6A_score_file = tmp_file_path +'/'+ tmp_file_name + '_up_M6Ascore.result'
     M6A_score_map = tmp_file_path +'/'+ tmp_file_name + '_up_M6Ascore_map'
     
     if pos_fa == 'default_path' or neg_fa == 'default_path':
          model_path = './requiredSoft/DeepM6ASeq/trained_models/hs/cnn/256-128_10-5_0_0.5_Y_Y_0.0_0.01_256/'
          os.system('python requiredSoft/DeepM6ASeq/main_test.py -m cnn -infa {} -md {} -outfn {}'.format(upstream_raw_file, model_path, M6A_score_file))

     else:
          out_dir = tmp_file_path +'/trained_models/'
          model_path = out_dir + 'cnn/256-128_10-5_0_0.5_Y_Y_0.0_0.01_256/'
          os.system('python requiredSoft/DeepM6ASeq/main_train.py -m cnn -pos_fa {} -neg_fa {} -od {}'.format(pos_fa, neg_fa, out_dir))
          os.system('python requiredSoft/DeepM6ASeq/main_test.py -m cnn -infa {} -md {} -outfn {}'.format(upstream_raw_file, model_path, M6A_score_file))
    
# filtrate from all the potential ORF fasta by comparing with the IRESscore file
def filter_ires_M6A(lenth_need, ires_score, m6a_score, tmp_file_name, tmp_file_path, final_name, final_file_path, pos_fa, neg_fa):
     
     han_raw_orf = open(orf_raw_file)
     han_IRES_score_file = open(IRES_score_file)
     han_M6A_score_file = open(M6A_score_file)
     dic_IRES = {}
     dic_M6A = {}

     # save the ires score
     for line in han_IRES_score_file:
          line = line.replace('\n', '')
          line = line.split('\t')
          if line[0] != 'ID':
               dic_IRES[line[0]] = line[2]

     # save the m6a score
     for line in han_M6A_score_file:
          line = line.replace('\n', '')
          line = line.split('\t')
          dic_M6A[line[0]] = line[1]

     han_IRES_M6A_score_file = open(tmp_file_path + '/' + tmp_file_name + '_up_IRES_M6A_score','w+')
     han_filter_score_file = open(tmp_file_path + '/' + tmp_file_name + '_up_IRES_M6A_filter_score','w+')

     # merge the ires and m6a score, obtain the orf with higher score
     for item in dic_IRES.keys():
          han_IRES_M6A_score_file.write(item + '\t' + dic_IRES[item] + '\t' + dic_M6A[item] + '\n')
          if eval(dic_IRES[item]) > ires_score or eval(dic_M6A[item]) > m6a_score:
               han_filter_score_file.write(item + '\t' + dic_IRES[item] + '\t' + dic_M6A[item] + '\n')

     han_raw_orf.close()
     han_IRES_score_file.close()
     han_M6A_score_file.close()
     han_IRES_M6A_score_file.close()
     han_filter_score_file.close()

     # filter the best orf in each circ
     han_filter_score_file = open(tmp_file_path + '/' + tmp_file_name + '_up_IRES_M6A_filter_score')
     
     circ_id_list = []
     file_list = han_filter_score_file.readlines()
     han_filter_score_file.close()

     # save the circ_id
     for line in file_list:
          circ_id = line.split('\t')[0].split('|')[0]
          if circ_id not in circ_id_list:
               circ_id_list.append(circ_id)

     # find the best orf in each circ
     whole_filter_orf_id_list = []
     for i in circ_id_list:
          orf_info_dic = {}
          orf_score_dic = {}
          for line in file_list:
               line = line.replace('\n', '')
               line_list = line.split('\t')
               high_score = eval(max(line_list[1], line_list[2]))
               orf_name = line_list[0]
               orf_score_dic[line_list[0]] = [line_list[1], line_list[2]]
               orf_name_list = orf_name.split('|')
               circ_id = orf_name_list[0]
               site = orf_name_list[-1]
               site = site.replace('start', '')
               site = site.replace('stop', '')
               start_site = site.split('-')[0]
               stop_site = site.split('-')[1]
               if circ_id == i:
                    if high_score not in orf_info_dic.keys():
                         orf_info_dic[high_score] = [start_site, orf_name]
                    elif high_score in orf_info_dic.keys():
                         high_score += 0.0000000000000001 
                         orf_info_dic[high_score] = [start_site, orf_name]
          list_key = list(orf_info_dic.keys())
          list_key.sort(reverse=True)

          need_list = []
          if len(list_key) == 1:
               need_list.append(list_key[0])
               whole_filter_orf_id_list.append(orf_info_dic[list_key[0]][1])
               continue

          need_list.append(list_key[0])
          for key in list_key[1:]:
               j = 0
               for item in need_list:
                    if abs(eval(orf_info_dic[key][0]) - eval(orf_info_dic[item][0])) <= lenth_need:
                         j += 1
               if j == 0:
                    need_list.append(key)

          
          
          for i in need_list:
               whole_filter_orf_id_list.append(orf_info_dic[i][1])
     
     # write the filter result into the new file
     han_raw_orf = open(orf_raw_file)

     IRES_M6A_filter_records = []

     seqs = SeqIO.parse(han_raw_orf, 'fasta')
     while True:
          try:
               seq_record = next(seqs)
               if seq_record.id in whole_filter_orf_id_list:
                    rec = SeqRecord(seq_record.seq,
                                    id=seq_record.id + '|' + orf_score_dic[seq_record.id][0][0:4] + '|' + orf_score_dic[seq_record.id][1][0:4],
                                    description='')
                    IRES_M6A_filter_records.append(rec)
          except StopIteration:
               break

     global orf_true_final_file
     orf_true_final_file = final_file_path +'/'+ final_name + '_true.fa'

     han_IRES_M6A_filter_file = open(orf_true_final_file, 'w+')
     SeqIO.write(IRES_M6A_filter_records, han_IRES_M6A_filter_file, 'fasta')
     han_raw_orf.close()
     han_IRES_M6A_filter_file.close()

     #make the final upstream ires file
     han_upstream_raw_file = open(upstream_raw_file)
     
     filter_upstream_ires_records = []

     seqs = SeqIO.parse(han_upstream_raw_file, 'fasta')
     while True:
          try:
               seq_record = next(seqs)
               if seq_record.id in whole_filter_orf_id_list:
                    rec = SeqRecord(seq_record.seq,
                                    id=seq_record.id + '|' + orf_score_dic[seq_record.id][0][0:4] + '|' + orf_score_dic[seq_record.id][1][0:4],
                                    description='')
                    filter_upstream_ires_records.append(rec)
          except StopIteration:
               break

     ires_final_file = tmp_file_path + '/' + tmp_file_name+'_final_IRES_up.fa'
     han_final_ires = open(ires_final_file, 'w+')
     SeqIO.write(filter_upstream_ires_records, han_final_ires, 'fasta')

     han_upstream_raw_file.close()
     han_final_ires.close()

     #give a score for each base to draw a picture
     M6A_score_map = tmp_file_path +'/'+ tmp_file_name + '_up_M6Ascore_map'
     
     if pos_fa == 'default_path' or neg_fa == 'default_path':
          model_path = './requiredSoft/DeepM6ASeq/trained_models/hs/cnn/256-128_10-5_0_0.5_Y_Y_0.0_0.01_256/'
          os.system('python requiredSoft/DeepM6ASeq/saliency_map.py -m cnn -infa {} -md {} -outfn {}'.format(ires_final_file, model_path, M6A_score_map))

     else:
          out_dir = tmp_file_path +'/trained_models/'
          model_path = out_dir + 'cnn/256-128_10-5_0_0.5_Y_Y_0.0_0.01_256/'
          os.system('python requiredSoft/DeepM6ASeq/saliency_map.py -m cnn -infa {} -md {} -outfn {}'.format(ires_final_file, model_path, M6A_score_map))

# Make a virtual sequence and reverse_complement the '-' strand
def reverse_complment(circrna_gtf, final_name, final_file_path):
     
     han_true = open(orf_true_final_file)
     virtual_file = final_file_path +'/'+ final_name + '_virtual.fa'
     han_circ_gtf = open(circrna_gtf)

     data_circ = []
     dic_strand = {}

     for line in han_circ_gtf:
          line = line.replace('\n', '')
          line_list = line.split('\t')
          dic_strand[line_list[0]] = line_list[1]


     for seq_record in SeqIO.parse(han_true, 'fasta'):
          circ_name = seq_record.id.split('$')[0]
          strand = dic_strand[circ_name]
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
     han_circ_gtf.close()

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

def get_time():
     time_now = time.strftime('[%Y-%m-%d %H:%M:%S]', time.localtime())
     return time_now

def main():
     parse = argparse.ArgumentParser(description='This script helps to filter coding orf')
     parse.add_argument('-y', '--yaml', required=True, help='please input the yamlfile')
     args = parse.parse_args()
     
     yamlfile = args.yaml
     con_file = open(yamlfile)
     fileload = yaml.full_load(con_file)
     
     circrna_gtf = fileload['circrna_gtf']
     raw_reads = fileload['raw_reads']
     tmp_file_path = fileload['tmp_file_location']
     final_file_path = fileload['result_file_location']
     pos_fa = fileload['pos_fa']
     neg_fa = fileload['neg_fa']
     ires_score = 0.5
     m6a_score = 0.5
     lenth_need = 101

     if not os.path.exists(tmp_file_path):
          os.makedirs(tmp_file_path)     
     
     if not os.path.exists(final_file_path):
          os.makedirs(final_file_path)

     begin = time.perf_counter()

     for item in raw_reads:

          genome_name = item.split('/')[-1][:-4]
          circrnas_file = fileload['tmp_file_location']+'/'+genome_name+'_filter.fa'

          print_file = genome_name+'_filter.fa'
          print(get_time(), 'analyzing {} ...'.format(print_file))

          # tmp file name and path
          tmp_file_name = genome_name + '_orf'

          # final file name and path
          final_name = genome_name + '_orf_filter_result'

          
          seq_pretreat(circrna_gtf, genome_name, circrnas_file, tmp_file_path)
          print(get_time(), 'seq_pretreat is finished!')

          get_1_orf_up(circrna_gtf, lenth_need, tmp_file_name, tmp_file_path)
          print(get_time(), 'get_1_orf_up is finished!')

          IRES_score(tmp_file_name, tmp_file_path)
          print(get_time(), 'IRES_score is finished!')

          M6A_score(tmp_file_name, tmp_file_path, pos_fa, neg_fa)
          print(get_time(), 'M6A_score is finished!')

          filter_ires_M6A(lenth_need, ires_score, m6a_score, tmp_file_name, tmp_file_path, final_name, final_file_path, pos_fa, neg_fa)
          print(get_time(), 'filter_ires_M6A is finished!')

          reverse_complment(circrna_gtf, final_name, final_file_path)
          print(get_time(), 'reverse_complment is finished!')

          translate_orf(final_name, final_file_path)
          print(get_time(), 'translate_orf is finished!')

     print(get_time(), 'the whole filter_coding_orf is finished! Spending time is {:.2f}s'.format(int(time.perf_counter() - begin)))
     print('{:=^23}'.format('End'))

if __name__ == '__main__':
     main()
     

