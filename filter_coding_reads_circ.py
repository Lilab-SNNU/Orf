from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import yaml
import argparse
import time
import os
import random

def make_mrna_rrna_index(transcript_fasta, ribosome_fasta, thread, tmp_file_path, raw_reads, soft_log_folder):

    read_name = str(raw_reads).split('/')[-1].split('.')[0]
    transcript_name = str(transcript_fasta).split('/')[-1].split('.')[0]
    ribosome_name = str(ribosome_fasta).split('/')[-1].split('.')[0]

    print(get_time(), 'make rRNA index...')
    rRNA_index_log = soft_log_folder + '/' + read_name +'_rRNA_index_log'
    os.system('bowtie-build --threads {} {} {}/{} > {}'.format(thread, ribosome_fasta, tmp_file_path, ribosome_name, rRNA_index_log))

    print(get_time(), 'make transcript index...')
    trans_index_log = soft_log_folder + '/' + read_name +'_transcprit_index_log'
    os.system('bowtie-build --threads {} {} {}/{} > {}'.format(thread, transcript_fasta, tmp_file_path, transcript_name, trans_index_log))

def predict_adapter(raw_reads, ribotype, tmp_file_path):

    read_name = str(raw_reads).split('/')[-1].split('.')[0]
    fastq_reads = raw_reads

    if ribotype == 'sra':
        print(get_time(), '{} sra to fastq...'.format(read_name))
        os.system('fastq-dump {} -O {}'.format(raw_reads, tmp_file_path))
        fastq_reads = tmp_file_path + '/' + read_name + '.fastq'

    fasta_reads = tmp_file_path + '/' + fastq_reads.split('/')[-1].split('.')[0]+'.fasta'

    # .fa is easy to predict adapter_start
    print(get_time(), 'fastq to fasta...')
    os.system('sed -n \'1~4s/^@/>/p;2~4p\' {} > {}'.format(fastq_reads, fasta_reads))

    # refer to these reads
    refer_reads_num = 5

    # obtain the index of refer reads
    print(get_time(), 'obtain refer reads...')
    han_fasta = open(fasta_reads)
    fasta_lines_num = int(len(han_fasta.readlines()))
    han_fasta.close()

    refer_reads_index = []
    while len(refer_reads_index) < refer_reads_num:
        reads_index = random.randrange(1, fasta_lines_num, 2)
        if reads_index not in refer_reads_index:
            refer_reads_index.append(reads_index)

    # window scanning and predict adapters
    print(get_time(), 'analysis refer reads...')
    predict_adapter_window = 8
    adapter_start_rate_list = [0.8, 0.7, 0.6]
    tmp_read_file = tmp_file_path + '/' + 'tmp_read.txt'
    tmp_count_file = tmp_file_path + '/' + 'tmp_count.txt'
    global adapter_start_set
    adapter_start_set = set()

    for adapter_start_rate in adapter_start_rate_list:

        for i in refer_reads_index:

            os.system('sed -n \'{},{}p\' {} > {}'.format(i, i+1, fasta_reads, tmp_read_file))
            han_tmp_read = open(tmp_read_file)
            for seq in SeqIO.parse(han_tmp_read, 'fasta'):

                for site in range(0,len(seq.seq)-8):

                    tmp_seq = str(seq.seq[site:site+8])

                    os.system('grep -c {} {} > {}'.format(tmp_seq, fasta_reads, tmp_count_file))
 
                    han_tmp_count = open(tmp_count_file)

                    tmp_seq_count = int(han_tmp_count.read())
                    han_tmp_count.close()

                    if tmp_seq_count/fasta_lines_num*2 >= adapter_start_rate:
                        adapter_start_set.add(tmp_seq)
                        break
            han_tmp_read.close()
        if len(adapter_start_set) != 0:
            break
    if len(adapter_start_set) == 0:
        print(get_time(), 'no adapter in these reads')
    elif len(adapter_start_set) != 0:
        print(get_time(), 'adapter_start in reads rate: {}'.format(adapter_start_rate))
        print(get_time(), 'adapter_start is: {}'.format(adapter_start_set))

def deal_raw_reads(raw_reads, tmp_file_path, thread, transcript_fasta, ribosome_fasta, trimmomatic, ribotype, soft_log_folder):

    read_name = str(raw_reads).split('/')[-1].split('.')[0]

    if ribotype == 'sra':
        fastq_reads = tmp_file_path + '/' + read_name + '.fastq'
    else:
        fastq_reads = raw_reads
    fastq_reads_cutadapt = tmp_file_path + '/' + read_name

    # cut predicted adapters
    print(get_time(), 'cut predicted adapters...')
    if len(adapter_start_set) != 0:
        adapter_command = ''
        fastq_reads_cutadapt = tmp_file_path + '/' + read_name + '_cutadapt'
        cutadapt_log = soft_log_folder + '/' + read_name + '_cutadapt_log'
        for item in adapter_start_set:
            adapter_command += ' -a ' + item
        os.system('cutadapt -j {}{} -o {} {} > {}'.format(thread, adapter_command, fastq_reads_cutadapt+ '.fastq', fastq_reads, cutadapt_log))

    # Filter out low quality reads by Trimmomatic
    print(get_time(), 'filter out low quality reads...')
    fastq_reads_cutadapt_clean = fastq_reads_cutadapt + '_clean'
    trim_log = soft_log_folder + '/' + read_name + '_trim_log'
    os.system('java -jar {} SE -threads {} {} {} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:16 > {}'.format(trimmomatic, thread, fastq_reads_cutadapt+'.fastq', fastq_reads_cutadapt_clean+'.fastq', trim_log))

    # Map clean reads to transcript sequence by bowtie
    print(get_time(), 'clean reads from transcript...')
    unmapped_transcript_reads = fastq_reads_cutadapt_clean + '_untranscipt'
    transcript_name = str(transcript_fasta).split('/')[-1].split('.')[0]
    transcript_index = tmp_file_path + '/' + transcript_name
    log_trans = soft_log_folder + '/' + transcript_name + '.map_to_transcript.sam'
    os.system('bowtie -p {} --un {} --norc {} {} > {}'.format(thread, unmapped_transcript_reads+'.fastq', transcript_index, fastq_reads_cutadapt_clean+'.fastq', log_trans))

    # Map clean reads to ribosome sequence by bowtie
    print(get_time(), 'clean reads from rRNA...')
    global unmapped_rRNA_reads
    unmapped_rRNA_reads = unmapped_transcript_reads + '_unrRNA'
    ribosome_name = str(ribosome_fasta).split('/')[-1].split('.')[0]
    rRNA_index = tmp_file_path + '/' + ribosome_name
    log_out = soft_log_folder + '/' + ribosome_name + '.map_to_rRNA.sam'
    os.system('bowtie -p {} --un {} --norc {} {} > {}'.format(thread, unmapped_rRNA_reads+'.fastq', rRNA_index, unmapped_transcript_reads+'.fastq', log_out))
    global final_reads
    final_reads = unmapped_rRNA_reads+'.fastq'


def make_junction_seq(raw_circrnas, circrna_gtf, raw_reads, tmp_file_path):

    print(get_time(), 'make junction seq...')
    # transform the final fastq to fasta, easy to filter
    global final_fasta_reads
    final_fasta_reads = unmapped_rRNA_reads + '.fasta'
    os.system('sed -n \'1~4s/^@/>/p;2~4p\' {} > {}'.format(final_reads, final_fasta_reads))

    max_reads_lenth = 0
    global reads_lenth_dic
    reads_lenth_dic = {}

    han_fasta_reads = open(final_fasta_reads)
    for seq_record in SeqIO.parse(han_fasta_reads, 'fasta'):
        reads_lenth_dic[seq_record.id] = len(seq_record.seq)
        if len(seq_record.seq) > max_reads_lenth:
            max_reads_lenth = len(seq_record.seq)

    # get the biggest length of reads, assemble the junction file 
    global junc_site
    junc_site = max_reads_lenth

    read_name = str(raw_reads).split('/')[-1].split('.')[0]
    circ_name = read_name + '_' + str(raw_circrnas).split('/')[-1].split('.')[0]

    han_raw_circ = open(raw_circrnas)
    han_circ_gtf = open(circrna_gtf)

    # reverse the '-' circRNA
    dic_strand = {}
    for line in han_circ_gtf:
        line = line.replace('\n', '')
        line_list = line.split('\t')
        dic_strand[line_list[0]] = line_list[1]

    global true_raw_file
    true_raw_file = tmp_file_path +'/'+ circ_name + '_true.fa'
    global true_junc_file
    true_junc_file = tmp_file_path + '/' + circ_name + '_true_junc.fa'
    data_circ = []
    data_junc = []
    junc_dic = {}

    # make the true seq file and junction file
    for seq_record in SeqIO.parse(han_raw_circ, 'fasta'):
        strand = dic_strand[seq_record.id]
        if strand == '-':
            seq_true = seq_record.seq.reverse_complement()
            rec_circ_true = SeqRecord(seq_true,
                                      id = seq_record.id,
                                      description = '')
            data_circ.append(rec_circ_true)

            seq_true_junc = seq_true[-max_reads_lenth:]+seq_true[:max_reads_lenth]
            rec_junc_true = SeqRecord(seq_true_junc,
                                      id = seq_record.id,
                                      description = '')
            data_junc.append(rec_junc_true)
        elif strand == '+':
            rec_circ_true = SeqRecord(seq_record.seq,
                                      id = seq_record.id,
                                      description = '')
            data_circ.append(rec_circ_true)

            seq_true_junc = seq_record.seq[-max_reads_lenth:]+seq_record.seq[:max_reads_lenth]
            rec_junc_true = SeqRecord(seq_true_junc,
                                      id = seq_record.id,
                                      description = '')
            data_junc.append(rec_junc_true)
    han_true = open(true_raw_file, 'w+')
    SeqIO.write(data_circ, han_true, 'fasta')
    han_junc = open(true_junc_file, 'w+')
    SeqIO.write(data_junc, han_junc, 'fasta')

    han_true.close()
    han_junc.close()
    han_raw_circ.close()
    han_circ_gtf.close()
    han_fasta_reads.close()

def make_junc_index_align(raw_reads, raw_circrnas, thread, tmp_file_path, soft_log_folder):

    read_name = str(raw_reads).split('/')[-1].split('.')[0]
    circ_name = read_name + '_' + str(raw_circrnas).split('/')[-1].split('.')[0]

    # make junction index
    print(get_time(), 'make junction index...')
    junc_index_log = soft_log_folder + '/' + read_name +'_junc_index_log'
    os.system('./requiredSoft/STAR --runMode genomeGenerate --runThreadN {} --genomeDir {} --genomeFastaFiles {} > {}'.format(thread, tmp_file_path+'/', true_junc_file, junc_index_log))

    # map the junction to filter reads
    print(get_time(), 'map the filter reads to the junction of circRNA...')
    global out_name
    out_name = tmp_file_path +'/'+ circ_name
    junc_map_reads_log = soft_log_folder + '/' + read_name +'_junc_map_to_reads_log'
    os.system('./requiredSoft/STAR --runThreadN {} --outSAMtype BAM SortedByCoordinate --genomeDir {} --readFilesIn {} --outFileNamePrefix {} > {}'.format(thread, tmp_file_path+'/', final_fasta_reads, out_name, junc_map_reads_log))

def bamtobed():

    print(get_time(), 'transform the bamflie to bedfile')
    bamfile = out_name + 'Aligned.sortedByCoord.out.bam'
    global bedfile
    bedfile = bamfile + '.bamtobedresult.bed'
    os.system('bedtools bamtobed -bed12 -i {} > {}'.format(bamfile, bedfile))

def filter_coding_circ_and_reads_index(raw_reads,tmp_file_path):

    print(get_time(), 'filter coding ablity by ribo-seq data...')
    coding_circ_set = set()
    coding_circ_reads_site = {}
    han_bed = open(bedfile)

    # save the circ which reads cross the junction, save the reads site
    for line in han_bed:
        list_line = line.split('\t')
        circ_id = list_line[0]
        read_id = list_line[3]
        reads_strand = list_line[5]
        start_site = eval(list_line[1]) + 1
        stop_site = eval(list_line[1]) + reads_lenth_dic[read_id]
        if start_site <= junc_site and stop_site >junc_site:
            coding_circ_set.add(circ_id)
            if circ_id in coding_circ_reads_site.keys():
                coding_circ_reads_site[circ_id] += [start_site-junc_site-1, stop_site-junc_site]
            if circ_id not in coding_circ_reads_site.keys():
                coding_circ_reads_site[circ_id] = [start_site-junc_site-1, stop_site-junc_site]

    # write the reads site to new file
    read_name = str(raw_reads).split('/')[-1].split('.')[0]
    reads_cover_jun_file = tmp_file_path + '/' + read_name + '_read_cover_jun.txt'
    han_reads_cover = open(reads_cover_jun_file, 'w+')
    for key in coding_circ_reads_site.keys():
        han_reads_cover.write(key + '\t' +str(int(len(coding_circ_reads_site[key])/2)))
        for item in coding_circ_reads_site[key]:
            han_reads_cover.write('\t'+str(item))
        han_reads_cover.write('\t'+ str(int(len(coding_circ_reads_site[key])/2)) +'\n')

    # write the circRNA with ribo_seq data to new file
    han_raw_circ = open(true_raw_file)
    read_name = str(raw_reads).split('/')[-1].split('.')[0]
    coding_circ_file = tmp_file_path + '/' + read_name + '_true.fa'
    coding_circ_seq_list = []
    for seq_record in SeqIO.parse(han_raw_circ, 'fasta'):
        if seq_record.id in coding_circ_set:
            junc_reads_count = str(int(len(coding_circ_reads_site[seq_record.id])/2))
            rec_coding = SeqRecord(seq_record.seq,
                            id=seq_record.id+'$'+junc_reads_count,
                            description='')
            coding_circ_seq_list.append(rec_coding)
    han_coding_circ = open(coding_circ_file,'w+')
    SeqIO.write(coding_circ_seq_list, han_coding_circ, 'fasta')
    han_bed.close()
    han_reads_cover.close()
    han_coding_circ.close()
        
def get_time():
    nowtime = time.strftime('[%Y-%m-%d %H:%M:%S]', time.localtime())
    return nowtime

def main():
    parse = argparse.ArgumentParser(description='This script helps to make junction index for each circ!')
    parse.add_argument('-y', '--yaml', required=True, help='please input the config.yaml file')
    args = parse.parse_args()

    yamlfile = args.yaml
    con_file = open(yamlfile)
    fileload = yaml.full_load(con_file)

    raw_reads = fileload['raw_reads']
    transcript_fasta = fileload['transcript_fasta']
    ribosome_fasta = fileload['ribosome_fasta']
    tmp_file_path = fileload['tmp_file_location']
    trimmomatic = fileload['trimmomatic_jar']
    thread = fileload['thread']
    ribotype = fileload['ribotype']
    raw_circrnas = fileload['circrnas']
    circrna_gtf = fileload['circrna_gtf']

    if not os.path.exists(tmp_file_path):
         os.makedirs(tmp_file_path)     

    soft_log_folder = tmp_file_path + '/' + 'log'

    if not os.path.exists(soft_log_folder):
         os.makedirs(soft_log_folder)

    for item in raw_reads:
        raw_reads = item
        predict_adapter(raw_reads, ribotype, tmp_file_path)
        make_mrna_rrna_index(transcript_fasta, ribosome_fasta, thread, tmp_file_path, raw_reads, soft_log_folder)
        deal_raw_reads(raw_reads, tmp_file_path, thread, transcript_fasta, ribosome_fasta, trimmomatic, ribotype, soft_log_folder)
        make_junction_seq(raw_circrnas, circrna_gtf, raw_reads, tmp_file_path)
        make_junc_index_align(raw_reads, raw_circrnas, thread, tmp_file_path, soft_log_folder)
        bamtobed()
        filter_coding_circ_and_reads_index(raw_reads,tmp_file_path)
        

if __name__ == '__main__':
    main()

