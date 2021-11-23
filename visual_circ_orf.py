from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import yaml
import argparse
import os
import time
from PIL import Image, ImageDraw, ImageFont
import math


# Draw the each IRES secondary structure of coding potetial ORF
def draw_IRES_secondary_structure(IRES_folder, tmp_file_name, tmp_file_path):
    IRES_graph_path = IRES_folder + '/'
    IRES_file = tmp_file_path + '/' + tmp_file_name+'_final_IRES_up.fa'

    sh_script = '''
                cd {}
                RNAfold {}
                '''.format(IRES_graph_path,IRES_file)
    sh_file = tmp_file_path+'/'+tmp_file_name+'_IRES_sh.r'
    han_sh_file = open(sh_file, 'w+')
    han_sh_file.write(sh_script)
    han_sh_file.close()

    os.system('sh {}'.format(sh_file))

# Make m6a site index for each orf
def make_m6a_site_index(m6a_score, tmp_file_name, tmp_file_path):

    M6A_score_map = tmp_file_path +'/'+ tmp_file_name + '_up_M6Ascore_map'
    han_M6Ascore_map = open(M6A_score_map)
    m6a_index_file = tmp_file_path + '/' + tmp_file_name +'_filter_M6A_site_index'
    han_m6a_index_file = open(m6a_index_file, 'w+')

    global dic_file
    dic_file ={}
    file_list = han_M6Ascore_map.readlines()
    
    for i in range(0, len(file_list), 4):
        site_list = []
        orf_name = file_list[i][11:-1]
        score = file_list[i+1][7:-1]
        nt_str = file_list[i+2].replace('\t', '').replace('\n', '')
        score_list = file_list[i+3].replace('\n', '').split('\t')
        if eval(score) >= m6a_score:
            for j in range(len(nt_str)):
                if nt_str[j:j+2] == 'AC' and eval(score_list[j]) >= eval(max(score_list))/2 and eval(score_list[j+1]) >= eval(max(score_list))/2:
                    site_list.append([str(j+1), str(j+2), nt_str[j-2:j+3]])
            dic_file[orf_name] = site_list
        lenth_site = len(site_list)
        han_m6a_index_file.write(orf_name)
        for k in site_list:
            han_m6a_index_file.write('\t'+k[0]+'\t'+k[1]+'\t'+k[2])     # true site
        han_m6a_index_file.write('\n')

    han_m6a_index_file.close()
    han_M6Ascore_map.close()

# Make the samples index and the Matrix
def make_index(circrna_gtf, tmp_file_name, final_name, tmp_file_path, final_file_path):

    # samples index building
    han_circ_gtf = open(circrna_gtf,'r+')
    han_orf = open(final_file_path + '/' + final_name + '_virtual.fa','r+')
    orf_list = []                                           # save orf long_id
    dic_circ_gtf = {}                                       # save {database_id : samples}
    new_file_list = []                                      # save long_id and samples

    for seq_record in SeqIO.parse(han_orf, 'fasta'):
        orf_list.append(seq_record.id)

    for line in han_circ_gtf:
        line = line.replace('\n', '')
        line_list = line.split('\t')
        if line_list[0] != 'circRNA_id':
            dic_circ_gtf[line_list[0]] = line_list[4]

    han_orf.close()
    han_circ_gtf.close()

    for i in orf_list:
        j = i.split('$')[0]
        if j in dic_circ_gtf.keys():
            new_file_list.append([i,dic_circ_gtf[j]])

    name = ['circRNA','samples']
    df = pd.DataFrame(columns=name, data=new_file_list)
    df.to_csv(tmp_file_path+'/'+tmp_file_name+'_visual_index.csv', mode='w', index=False, sep='\t')

    # Matrix_building()
    global name_set
    name_set = set()                                          # save raw all samples

    for i in dic_circ_gtf.values():
        for item in i.split(','):
            name_set.add(item)
    name_list = list(name_set)
    name_list.sort()

    global num_column
    num_column = len(name_list)
    matrix_list = []

    for i in new_file_list:
        tmp_list = []                                      # save each data matrix
        tmp_list.append(i[0])
        for j in name_list:
            if j in i[1].split(','):
                tmp_list.append(1)
            elif j not in i[1].split(','):
                tmp_list.append(0)
        matrix_list.append(tmp_list)
    name_final_list = ['circRNA'] + name_list
    df_matrix = pd.DataFrame(columns=name_final_list, data=matrix_list)
    df_matrix.to_csv(tmp_file_path+'/'+tmp_file_name+'_Matrix.csv', mode='w', index=False, sep='\t')

# Draw the circ graph of each coding potential ORF
def make_circ_pic(circrna_gtf, CIRC_folder, lenth_need, tmp_file_name, tmp_file_path, final_name, final_file_path, genome_name):

    def get_angle(bp, length):
        return bp * 360 / length

    def coord(angle, center, radius):
        rad = math.radians(angle)
        x = int(center[0] + math.cos(rad) * radius)
        y = int(center[1] + math.sin(rad) * radius)
        return x, y

    def draw_circle(radius_i, radius_o, start_angle, stop_angle, color1, color2):
        lenth = size[0]
        width = size[1]
    
        x1 = lenth/2 - radius_o
        y1 = CENTER[1] - radius_o
        x2 = lenth/2 + radius_o
        y2 = CENTER[1] + radius_o
        DRAW.pieslice((x1, y1, x2, y2), start_angle, stop_angle, fill=color1)
    
        x3 = lenth/2 - radius_i
        y3 = CENTER[1] - radius_i
        x4 = lenth/2 + radius_i
        y4 = CENTER[1] + radius_i
        DRAW.pieslice((x3, y3, x4, y4), start_angle, stop_angle, fill=color2)

    def draw_arrow_tip(start, direction, radius_mid, color):
        p1 = coord(start + direction, CENTER, radius_mid)
        p2 = coord(start, CENTER, radius_mid-100)
        p3 = coord(start, CENTER, radius_mid+100)
        DRAW.polygon((p1, p2, p3), fill=color)
        

    han_orf = open(final_file_path + '/' + final_name + '_virtual.fa','r+')
    han_circ_gtf = open(circrna_gtf,'r+')

    dic_gtf = {}
    for line in han_circ_gtf:
        line = line.replace('\n', '')
        line_list = line.split('\t')
        dic_gtf[line_list[0]] = line_list[1:]

    not_in_gtf_list = []
    for seq_record in SeqIO.parse(han_orf, 'fasta'):

        site = seq_record.id.split('|')[-3]
        site = site.replace('start', '')
        site = site.replace('stop', '')
        site_list = site.split('-')

        circ_name = seq_record.id.split('$')[0]

        # new file name
        orf_id = circ_name + '|' + site
            
        strand = dic_gtf[circ_name][0]
        

        # in case that circ name not in the gtf
        if circ_name not in dic_gtf.keys():
            not_in_gtf_list.append(circ_name+'\n')
            continue
        circ_lenth = eval(dic_gtf[circ_name][1])
        orf_lenth = len(seq_record.seq)
        orf_start_site = eval(site_list[0])
        orf_stop_site = eval(site_list[1])
        ires_score = seq_record.id.split('|')[-2]
        m6a_score = seq_record.id.split('|')[-1]
        gene_id = dic_gtf[circ_name][2]

        if strand == '+':
            orf_draw_start_site = orf_start_site
            orf_draw_stop_site = orf_stop_site
            ires_draw_start_site = orf_start_site - lenth_need
            ires_draw_stop_site = orf_start_site - 1
        elif strand == '-':
            orf_draw_start_site = orf_stop_site           # not true site, draw start site
            orf_draw_stop_site = orf_start_site
            ires_draw_start_site = orf_start_site + 1
            ires_draw_stop_site = orf_start_site + lenth_need

        # set the canvas and the center
        size = (5000, 5000)
        CENTER = (2500, 2500)

        myseq = Image.new('RGB', size, 'white')
        DRAW = ImageDraw.Draw(myseq)

        # draw start angle and draw stop angle
        ORF_START, ORF_END = get_angle(orf_draw_start_site, circ_lenth), get_angle(orf_draw_stop_site, circ_lenth)
        IRES_START, IRES_END = get_angle(ires_draw_start_site, circ_lenth), get_angle(ires_draw_stop_site, circ_lenth)

        # caluate
        if lenth_need*360/circ_lenth >= 15:
            if strand == '+':
                draw_circle(1300, 1400, IRES_START, IRES_END-10, 'red', 'white')
                draw_arrow_tip(IRES_END-10, 10, 1350, 'red')
            elif strand == '-':
                draw_circle(1300, 1400, IRES_START+10, IRES_END, 'red', 'white')
                draw_arrow_tip(IRES_START+10, -10, 1350, 'red')
        else:
            draw_circle(1300, 1400, IRES_START, IRES_END, 'red', 'white')
 
        if orf_lenth*360/circ_lenth >= 15:
            if strand == '+':
                draw_circle(1100, 1200, ORF_START, ORF_END-10, 'orange', 'white')
                draw_arrow_tip(ORF_END-10, 10, 1150, 'orange')
            elif strand == '-':
                draw_circle(1100, 1200, ORF_START+10, ORF_END, 'orange', 'white')
                draw_arrow_tip(ORF_START+10, -10, 1150, 'orange')
        else:
            draw_circle(1100, 1200, ORF_START, ORF_END, 'orange', 'white')
            
             
        # draw the circRNA sequence
        draw_circle(800, 1000, 0, 360, 'lightgrey', 'white')

        # draw the junction site
        p_i = coord(0, CENTER, 800)
        p_o = coord(0, CENTER, 1000)
        DRAW.line((p_i, p_o), fill= 'white', width=30)

        # draw m6a site
        if seq_record.id in dic_file.keys():
            if dic_file[seq_record.id] != []:
                for item in dic_file[seq_record.id]:
                    m6a_start = eval(item[0])
                    if strand == '-':
                        m6a_site = ires_draw_stop_site - m6a_start
                    elif strand == '+':
                        m6a_site = ires_draw_start_site + m6a_start

                    m6a_angle_site = get_angle(m6a_site, circ_lenth)
                    m_i = coord(m6a_angle_site, CENTER, 800)
                    m_o = coord(m6a_angle_site, CENTER, 1000)
                    DRAW.line((m_i, m_o), fill= 'green', width=30)

        # rotate the picture so that the junction is facing up
        myseq = myseq.rotate(90)
        DRAW = ImageDraw.Draw(myseq)

        # Linux type font
        arial100 = ImageFont.truetype('LiberationSans-Regular.ttf', 100)

        # write the title of circRNA name
        arial200 = ImageFont.truetype('LiberationSans-Regular.ttf', 200)
        DRAW.text((1400, 4300), '{}({})'.format(orf_id, strand), fill='black', font = arial200)

        # write the annoation information
        DRAW.rectangle((3600, 300, 3800, 400), fill='lightgrey')
        DRAW.text((3900, 300), 'ref spliced circ(+)', fill='black', font = arial100)
        DRAW.rectangle((3600, 500, 3800, 600), fill='red')
        DRAW.text((3900, 500), 'IRES score : {}'.format(ires_score), fill='black', font = arial100)
        DRAW.rectangle((3600, 700, 3800, 800), fill='green')
        DRAW.text((3900, 700), 'm6A score : {}'.format(m6a_score), fill='black', font = arial100)
        

        DRAW.text((1900,2000), 'circRNA_lenth = {}'.format(circ_lenth), fill='black', font = arial100)
        DRAW.text((1900,2200), 'ORF_lenth = {}'.format(orf_lenth), fill='black', font = arial100)
        DRAW.text((1900,2400), 'gene_id = {}'.format(gene_id), fill='black', font = arial100)
        DRAW.text((1900,2600), 'ORF start site : {}'.format(orf_start_site), fill='black', font = arial100)
        DRAW.text((1900,2800), 'ORF stop site : {}'.format(orf_stop_site), fill='black', font = arial100)
        
        arial150 = ImageFont.truetype('LiberationSans-Regular.ttf', 150)
        DRAW.text((2280,1800), "3'", fill='black', font = arial150)
        DRAW.text((2600,1800), "5'", fill='black', font = arial150)

        new_name = orf_id
        myseq.save('{}/{}.png'.format(CIRC_folder, new_name))

    han_not = open(tmp_file_path + '/' + tmp_file_name + '_not_in_gtf_list','w+')
    han_not.writelines(not_in_gtf_list)
    han_not.close()

    han_orf.close()
    han_circ_gtf.close()

# Draw the word graph of each coding potential ORF
def make_circ_pic2(circrna_gtf, lenth_need, tmp_file_name, tmp_file_path, final_name, final_file_path, genome_name, WORD_folder):
    def get_angle(bp, lenth):
        return bp * 360 / lenth

    def coord(angle, center, radius):
        rad = math.radians(angle)
        x = int(center[0] + math.cos(rad) * radius)
        y = int(center[1] + math.sin(rad) * radius)
        return x, y

    def draw_circle(radius_i, radius_o, start_angle, stop_angle, color1, color2):
        lenth = size[0]
        width = size[1]
    
        x1 = lenth/2 - radius_o
        y1 = CENTER[1] - radius_o
        x2 = lenth/2 + radius_o
        y2 = CENTER[1] + radius_o
        DRAW.pieslice((x1, y1, x2, y2), start_angle, stop_angle, fill=color1)
    
        x3 = lenth/2 - radius_i
        y3 = CENTER[1] - radius_i
        x4 = lenth/2 + radius_i
        y4 = CENTER[1] + radius_i
        DRAW.pieslice((x3, y3, x4, y4), start_angle, stop_angle, fill=color2)
    
    def draw_arrow_tip(start, direction, radius_mid, color):
        p1 = coord(start + direction, CENTER, radius_mid)
        p2 = coord(start, CENTER, radius_mid-50)
        p3 = coord(start, CENTER, radius_mid+50)
        DRAW.polygon((p1, p2, p3), fill=color)

    # draw the word_graph
    han_true_ires = open(tmp_file_path +'/'+ tmp_file_name + '_up.fa', 'r+')
    han_true_orf = open(final_file_path + '/' + final_name + '_true.fa', 'r+')
    han_circ_raw = open(tmp_file_path + '/' + genome_name+'_filter.fa', 'r+')
    han_amino = open(final_file_path +'/'+ final_name + '_true_translated.fa', 'r+')
    han_reads_cover = open(tmp_file_path + '/' + genome_name + '_read_cover_jun.txt','r+')
    han_circ_gtf = open(circrna_gtf,'r+')

    dic_gtf = {}
    for line in han_circ_gtf:
        line = line.replace('\n', '')
        line_list = line.split('\t')
        dic_gtf[line_list[0]] = line_list[1:]


    orf_dic={}
    ires_dic={}
    circ_dic={}
    trans_dic={}

    for seq1 in SeqIO.parse(han_true_orf, 'fasta'):
        orf_dic[seq1.id] = str(seq1.seq)
    for seq2 in SeqIO.parse(han_circ_raw, 'fasta'):
        circ_dic[seq2.id] = str(seq2.seq)
    for seq3 in SeqIO.parse(han_amino, 'fasta'):
        trans_dic[seq3.id] = str(seq3.seq)
    for seq4 in SeqIO.parse(han_true_ires, 'fasta'):
        ires_dic[seq4.id] = str(seq4.seq)
    line_reads_list = []
    for line in han_reads_cover:
        line = line.replace('\n', '')
        line = line.split('\t')
        line_reads_list.append(line)    



    for id1 in orf_dic.keys():

        name_circ = id1.split('$')[0]
        str_seq = orf_dic[id1]
        for id2 in circ_dic.keys():
            circ_name = id2.split('$')[0]
            if name_circ == circ_name:
                str_circ = circ_dic[id2]
        for id3 in trans_dic.keys():
            if id3 == id1:
                str_trans = trans_dic[id3]
        for id4 in ires_dic.keys():
            if id4 == id1[0:-10]:
                str_ires = ires_dic[id4]

        site = id1.split('|')[-3]
        site = site.replace('start', '')
        site = site.replace('stop', '')
        site_list = site.split('-')
        orf_start_point = eval(site_list[0])
        orf_stop_point = eval(site_list[1])
        strand = dic_gtf[name_circ][0]

        # new file name
        orf_id = name_circ + '|' + site

        # sequence and start point
        if strand == '-':
            orf_start_site = orf_stop_point - 2
            ires_start_site = orf_start_point 
            strseq = str_seq[::-1]
            strtrans = str_trans[::-1]
            strires = str_ires[::-1]
        elif strand =='+':
            orf_start_site = orf_start_point
            ires_start_site = orf_start_point- lenth_need
            strseq = str_seq
            strtrans = str_trans
            strires = str_ires
        
        
        lenth0 = len(strtrans)
        lenth1 = len(strseq)
        lenth2 = len(strires)
        strcirc = str(str_circ)
        lenth3 = len(strcirc)        
        move_num = lenth3/4
        angle = 360/lenth3
        angle1 = 359/lenth3
        
        
        size = (5000, 5000)
        CENTER = (2500, 2500)

        myseq = Image.new('RGB', size, 'white')
        DRAW = ImageDraw.Draw(myseq)
        if strand == '+':
            start_arrow = get_angle(-move_num, lenth3)
            end_arrow = get_angle(-move_num+10, lenth3)
            draw_circle(1950, 2000, start_arrow, end_arrow, 'orange', 'white')
            draw_arrow_tip(end_arrow, 5, 1975, 'orange')
        elif strand == '-':
            start_arrow = get_angle(-move_num-10, lenth3)
            end_arrow = get_angle(-move_num, lenth3)
            draw_circle(1950, 2000, start_arrow, end_arrow, 'orange', 'white')
            draw_arrow_tip(start_arrow, -5, 1975, 'orange')

        p_i = coord((-move_num+0.5)*angle, CENTER, 2100)
        p_o = coord((-move_num+0.5)*angle, CENTER, 2250)
        DRAW.line((p_i, p_o), fill= 'green', width=20)

        if lenth3 <= 500:
            type_num = 40
        elif 500 < lenth3 <= 1500:
            type_num = 20
        elif 1500 < lenth3 <= 2500:
            type_num = 10
        elif 2500 < lenth3:
            type_num = 5
        arial100 = ImageFont.truetype('LiberationSans-Regular.ttf', type_num)

        for i in range(lenth0):
            if i+orf_start_site<=lenth3:
                p = coord((i*3+orf_start_site-move_num)*angle, CENTER, 2100)
                DRAW.text(p, strtrans[i], fill = 'green', font=arial100)
            elif i+orf_start_site>lenth3:
                p = coord((i*3+orf_start_site-(move_num))*angle, CENTER, 2100)
                DRAW.text(p, strtrans[i], fill = 'green', font=arial100)



        for i in range(lenth3):
            p = coord((i+1-move_num)*angle1, CENTER, 2200)
            if i < 10:
                colorcirc = 'black'
            else:
                colorcirc = 'black'
            DRAW.text(p, strcirc[i], fill = colorcirc, font=arial100)


        for i in range(lenth1):
            if i+orf_start_site<=lenth3:
                p = coord((i+orf_start_site-move_num)*angle, CENTER, 2300)
                DRAW.text(p, strseq[i], fill = 'blue', font=arial100)
            elif i+orf_start_site>lenth3:
                p = coord((i+orf_start_site-(move_num))*angle, CENTER, 2300)
                DRAW.text(p, strseq[i], fill = 'blue', font=arial100)


        for i in range(lenth2):
            if i+ires_start_site<=lenth3:
                p = coord((i+ires_start_site-move_num)*angle, CENTER, 2400)
                DRAW.text(p, strires[i], fill = 'red', font=arial100)
            elif i+ires_start_site>lenth3:
                p = coord((i+ires_start_site-(move_num))*angle, CENTER, 2400)
                DRAW.text(p, strires[i], fill = 'red', font=arial100)

        # draw m6a site
        if id1 in dic_file.keys():
            for item in dic_file[id1]:
                m6a_start = eval(item[0])
                strm6a = item[2]
                if strand == '+':
                    m6a_site = ires_start_site + m6a_start
                elif strand == '-':
                    m6a_site = ires_start_site + lenth2 - m6a_start
                    strm6a = strm6a[::-1]
                for i in range(len(strm6a)):
                    q = coord((i+m6a_site-move_num)*angle, CENTER, 2450)
                    DRAW.text(q, strm6a[i], fill = 'purple', font=arial100)

        # reads
        for line in line_reads_list:
            
            if line[0] == name_circ and int(line[-1]) <= 20:
                reads_list = line[2:-1]
                lenth_reads = len(reads_list)
                for it in range(0,lenth_reads,2):
                    n_small = int((it + 2)/2)
                    in_r = 1880 - 40*n_small
                    out_r = 1900 - 40*n_small
                    left_a = (int(reads_list[it])-move_num)*angle
                    right_a = (int(reads_list[it+1])-move_num+1)*angle
                    draw_circle(in_r, out_r, left_a, right_a, 'pink', 'white')

        new_name = orf_id
        myseq.save('{}/{}_word.png'.format(WORD_folder, new_name))



    han_true_ires.close()
    han_true_orf.close()
    han_circ_raw.close()
    han_amino.close()
    han_reads_cover.close()
    han_circ_gtf.close()

# Statistics the distribution of all result data
def UpsetR(final_name, tmp_file_name, tmp_file_path, final_file_path):
    
    Matrix_file = tmp_file_path+'/'+tmp_file_name+'_Matrix.csv'
    num_list = (1,1,1,1,1,1)
    final_file_Name = tmp_file_name+'_samples.pdf' 
    
    r_script = '''
    library(UpSetR)
    setwd("{}")
    pdf(file='{}',onefile=FALSE)
    example = read.csv("{}",header=TRUE,row.names=1,check.names = FALSE)
    upset(example, mb.ratio = c(0.55, 0.45), order.by = "freq",
    nsets = {}, number.angles = 0, point.size = 1.2, line.size = 0.5,
    mainbar.y.label = "Intersection",sets.x.label = "Frequency count", text.scale = c{})
    while(!is.null(dev.list())) dev.off()
    '''.format(final_file_path, final_file_Name, Matrix_file, num_column, num_list)

    r_file = tmp_file_path+'/'+tmp_file_name+'_rscript.r'
    han_r_file = open(r_file, 'w+')
    han_r_file.write(r_script)
    han_r_file.close()
    os.system('Rscript {}'.format(r_file))

# draw the express_analysis pdf
def express_analysis(raw_reads, tmp_file_path, final_file_path):
    
    orf1_name = raw_reads[0].split('/')[-1][:-4]
    orf2_name = raw_reads[1].split('/')[-1][:-4]

    han_orf1 = open(final_file_path + '/' + orf1_name + '_orf_filter_result' + '_virtual.fa','r+')
    han_orf2 = open(final_file_path + '/' + orf2_name + '_orf_filter_result' + '_virtual.fa','r+')
    han1_dic = {}
    han2_dic = {}

    for seq_record in SeqIO.parse(han_orf1, 'fasta'):

        seq_record.id = seq_record.id.replace('$', '|')

        circ_name = seq_record.id.split('|')[0]
        reads_count = seq_record.id.split('|')[-5]
        han1_dic[circ_name] = reads_count

    for seq_record in SeqIO.parse(han_orf2, 'fasta'):

        seq_record.id = seq_record.id.replace('$', '|')

        circ_name = seq_record.id.split('|')[0]
        reads_count = seq_record.id.split('|')[-5]
        han2_dic[circ_name] = reads_count

    union_set = han1_dic.keys()|han2_dic.keys()
    circ_list = list(union_set)
    circ_list.sort()
    file_list = []
    for item in circ_list:
        file_list.append(item+'\t')
        if item in han1_dic.keys():
            file_list.append(han1_dic[item]+'\t')
        else:
            file_list.append('0'+'\t')
        if item in han2_dic.keys():
            file_list.append(han2_dic[item]+'\n')
        else:
            file_list.append('0'+'\n')    

    index_file = tmp_file_path + '/two_samples_express.txt'
    han3 = open(index_file, 'w+')
    han3.writelines(['*\t', orf1_name+'\t', orf2_name+'\n'])
    han3.writelines(file_list)
    han_orf1.close()
    han_orf2.close()
    han3.close()


    output_pdf = final_file_path + '/express_analysis.pdf'
    output_csv = tmp_file_path + '/express_score.csv'
    r_script = '''
    library(edgeR)
    library(ggplot2)
    data2 <- read.csv('{}',sep='\t',row.names=1)
    counts <- data2[,c(1,2)]
    group <- c(1,2)
    y <- DGEList(counts=counts,group=group)
    keep <- rowSums(cpm(y)>1) >=1
    y <- y[keep,,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    y_bcv <- y
    bcv <- 0.4
    et <- exactTest(y_bcv,dispersion = bcv^2)
    genel <- decideTestsDGE(et,p.value=0.05,lfc=0)
    df <- y_bcv$table
    results <- cbind(y$counts,et$table,genel)
    summary(genel)
    write.csv(x=results,file='{}')
    cut_off_pvalue = 0.05
    cut_off_logFC = 1
    results$change = ifelse(results$PValue < cut_off_pvalue & abs(results$logFC)>= cut_off_logFC,ifelse(results$logFC>cut_off_logFC,'Up','Dpwn'),'Stable')
    pdf('{}')
    ggplot(results,aes(x = logFC, y = -log10(PValue), colour=change)) +
    geom_point(alpha=0.4, size=3.5) +
    scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",y="-log10 (p-value)")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),legend.position="right",legend.title = element_blank())
    dev.off()
    '''.format(index_file, output_csv, output_pdf)
    
    r_file = tmp_file_path+'/' + 'express_analysis_rscript.r'
    han_r_file = open(r_file, 'w+')
    han_r_file.write(r_script)
    han_r_file.close()
    os.system('Rscript {}'.format(r_file))


def get_time():
    time_now = time.strftime('[%Y-%m-%d %H:%M:%S]', time.localtime())
    return time_now

def main():

    parse = argparse.ArgumentParser(description='This script helps to visualize the circ')
    parse.add_argument('-y', '--yaml', required=True, help='please input the yamlfile')
    args = parse.parse_args()
     
    yamlfile = args.yaml
    con_file = open(yamlfile)
    fileload = yaml.full_load(con_file)

    circrna_gtf = fileload['circrna_gtf']
    raw_reads = fileload['raw_reads']
    tmp_file_path = fileload['tmp_file_location']
    final_file_path = fileload['result_file_location']

    lenth_need = 101
    m6a_score = 0.5

    begin = time.perf_counter()

    for item in raw_reads:
          
        genome_name = item.split('/')[-1][:-4]
        circrnas_file = fileload['tmp_file_location']+'/'+genome_name+'_filter.fa'

        print(get_time(), 'analyzing {} ...'.format(circrnas_file))

        IRES_folder = '{}/{}_IRES_graph'.format(final_file_path, genome_name)
        CIRC_folder = '{}/{}_CIRC_graph'.format(final_file_path, genome_name)
        WORD_folder = '{}/{}_WORD_graph'.format(final_file_path, genome_name)
        if not os.path.exists(IRES_folder):
            os.makedirs(IRES_folder)

        if not os.path.exists(CIRC_folder):
            os.makedirs(CIRC_folder)

        if not os.path.exists(WORD_folder):
            os.makedirs(WORD_folder)

        # tmp file name and path
        tmp_file_name = genome_name + '_orf'
     
        # final file name and path
        final_name = genome_name + '_orf_filter_result'
      
        draw_IRES_secondary_structure(IRES_folder, tmp_file_name, tmp_file_path)
        print(get_time(), 'draw_IRES_secondary_structure is finished!')

        make_m6a_site_index(m6a_score, tmp_file_name, tmp_file_path)
        print(get_time(), 'make_m6a_site_index is finished!')

        make_index(circrna_gtf, tmp_file_name, final_name, tmp_file_path, final_file_path)
        print(get_time(), 'make_index is finished!')

        if len(name_set) > 5:
            UpsetR(final_name, tmp_file_name, tmp_file_path, final_file_path)
            print(get_time(), 'UpsetR finished!')
        else:
            print(get_time(), 'the num of samples are too few to draw UpsetR!!!')

        make_circ_pic(circrna_gtf, CIRC_folder, lenth_need, tmp_file_name, tmp_file_path, final_name, final_file_path, genome_name)
        print(get_time(), 'make_circ_pic is finished!')
        
        make_circ_pic2(circrna_gtf, lenth_need, tmp_file_name, tmp_file_path, final_name, final_file_path, genome_name, WORD_folder)
        print(get_time(), 'make_circ_seq_pic is finished!')

    if len(raw_reads) == 2:
        time4 = time.perf_counter()
        express_analysis(raw_reads, tmp_file_path, final_file_path)
        time_add4 = time.perf_counter()-time4
        print(get_time(), 'express_analysis is finished!')
    print(get_time(), 'the whole visual_circ_orf is finished! Spending time is {:.2f}s'.format(int(time.perf_counter() - begin)))
    print('{:=^23}'.format('End'))

if __name__ == '__main__':
    main()







