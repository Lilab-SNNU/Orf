from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import yaml
import argparse
import os
import time
from PIL import Image, ImageDraw, ImageFont
import math

def file_load():
     
     parse = argparse.ArgumentParser(description='Welcome to this project')
     parse.add_argument('-y', '--yaml', required=True, help='please input the yamlfile')
     args = parse.parse_args()
     
     yamlfile = args.yaml
     con_file = open(yamlfile)
     fileload = yaml.full_load(con_file)
     
     lenth_need = fileload['lenth_need']
     circrnas_file = fileload['result_file_location']+'/'+fileload['genome_name']+'.fa'
     circrna_gtf = fileload['circrna_gtf']

     # tmp file name and path
     tmp_file_name = fileload['genome_name'] + '_orf'
     tmp_file_path = fileload['tmp_file_location']

     # final file name and path
     final_name = fileload['genome_name'] + '_orf_filter_result'
     final_file_path = fileload['result_file_location']

     return lenth_need, circrnas_file, circrna_gtf, tmp_file_name, tmp_file_path, final_name, final_file_path    

# Draw the each IRES secondary structure of coding potetial ORF
def draw_IRES_secondary_structure(tmp_file_name, tmp_file_path, final_name, final_file_path):

    # make the final IRES file
    han_orf = open(final_file_path + '/' + final_name + '_virtual.fa','r+')
    han_raw_IRES = open(tmp_file_path + '/' + tmp_file_name+'_up.fa')
    
    id_dic = {}
    final_ires_records = []
    for seq_record in SeqIO.parse(han_orf, 'fasta'):
        seq_id = seq_record.id.split('|')
        seq_name = '{}|{}|{}|{}|{}|{}'.format(seq_id[0],seq_id[1],seq_id[2],seq_id[3],seq_id[4],seq_id[5])
        id_dic[seq_name] = seq_record.id
    for sequence in SeqIO.parse(han_raw_IRES, 'fasta'):
        if sequence.id in id_dic.keys():
            rec = SeqRecord(sequence.seq,
                            id=str(id_dic[sequence.id]),
                            description='')
            final_ires_records.append(rec)
    han_final_ires = open(tmp_file_path + '/' + tmp_file_name+'final_IRES_up.fa', 'w+')
    SeqIO.write(final_ires_records, han_final_ires, 'fasta')
    
    han_orf.close()
    han_raw_IRES.close()
    han_final_ires.close()


    IRES_graph_path = final_file_path + '/IRES_graph'
    IRES_file = tmp_file_path + '/' + tmp_file_name+'final_IRES_up.fa'
    
    sh_script = '''
    cd {}
    mkdir IRES_graph
    cd {}
    RNAfold {}
    '''.format(final_file_path, IRES_graph_path, IRES_file)

    sh_file = tmp_file_path + '/' + tmp_file_name + '_shscript.sh'
    han_sh_file = open(sh_file,'w+')
    han_sh_file.write(sh_script)
    han_sh_file.close()
    os.system('sh {}'.format(sh_file))
    print(len(id_dic))

# Make the samples index and the Matrix
def make_index(circrna_gtf, tmp_file_name, final_name, tmp_file_path, final_file_path):

    # samples index building
    han_circ_gtf = open(circrna_gtf,'r+')
    han_orf = open(final_file_path + '/' + final_name + '_virtual.fa','r+')
    orf_list = []                                           # save orf long_id
    dic_circ_gtf = {}                                       # save {database_id : samples}
    new_file_list = []                                      # save long_id and samples

    global dic_circ_lenth
    dic_circ_lenth = {}
    for seq_record in SeqIO.parse(han_orf, 'fasta'):
        orf_list.append(seq_record.id)
    
    for line in han_circ_gtf:
        line = line.replace('\t','$')
        line = line.replace('\n','')
        line_list = line.split('$')
        line_str = line_list[7].replace(' ','')
        dic_circ_gtf[line_list[4]] = line_str
        dic_circ_lenth[line_list[4]] = line_list[6]
    
    for i in orf_list:
        j = i.split('|')[0]
        if j in dic_circ_gtf.keys():
            new_file_list.append([i,dic_circ_gtf[j]])
    
    name = ['circRNA','samples']
    df = pd.DataFrame(columns=name, data=new_file_list)
    df.to_csv(tmp_file_path+'/'+tmp_file_name+'_visual_index.csv', mode='w', index=False)

    # Matrix_building()
    name_ori = []                                          # save raw all samples

    for i in dic_circ_gtf.values():
        name_ori +=  i.split(',')
    name_s = set(name_ori[1:])
    name_list = list(name_s)

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
    df_matrix.to_csv(tmp_file_path+'/'+tmp_file_name+'_Matrix.csv', mode='w', index=False)

# Draw the circ graph of each coding potential ORF
def make_circ_pic(lenth_need, tmp_file_name, tmp_file_path, final_name, final_file_path):

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
        


    sh_script = '''
    cd {}
    mkdir CIRC_graph
    '''.format(final_file_path)
    sh_file = tmp_file_path + '/' + tmp_file_name + '_png_shscript.sh'
    han_sh_file = open(sh_file,'w+')
    han_sh_file.write(sh_script)
    han_sh_file.close()
    os.system('sh {}'.format(sh_file))



    han_orf = open(final_file_path + '/' + final_name + '_virtual.fa','r+')

    
    not_in_gtf_list = []
    for seq_record in SeqIO.parse(han_orf, 'fasta'):
        
             
         
        site = seq_record.id.split('|')[4]
        site = site.replace('start', '')
        site = site.replace('stop', '')
        site_list = site.split('-')

        strand = seq_record.id.split('|')[1][-1]
        circ_name = seq_record.id.split('|')[0]

        # in case that circ name not in the gtf
        
        if circ_name not in dic_circ_lenth.keys():
            not_in_gtf_list.append(circ_name+'\n')
            continue
        circ_lenth = eval(dic_circ_lenth[circ_name])
        orf_lenth = len(seq_record.seq)
        orf_start_site = eval(site_list[0])
        orf_stop_site = eval(site_list[1])
        ires_score = seq_record.id.split('|')[6]
        orf_score = seq_record.id.split('|')[7]
        transcript_id = seq_record.id.split('|')[2]

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
            #draw_arrow_tip(IRES_START-10, -10, 1350, 'red')
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

        # rotate the picture so that the junction is facing up
        myseq = myseq.rotate(90)
        DRAW = ImageDraw.Draw(myseq)

        # Linux type font
        arial100 = ImageFont.truetype('LiberationSans-Regular.ttf', 100)

        # write the title of circRNA name
        arial200 = ImageFont.truetype('LiberationSans-Regular.ttf', 200)
        DRAW.text((1400, 700), '{}({})'.format(circ_name, strand), fill='black', font = arial200)

        # write the annoation information
        DRAW.rectangle((3600, 300, 3800, 400), fill='lightgrey')
        DRAW.text((3900, 300), 'ref spliced circ(+)', fill='black', font = arial100)
        DRAW.rectangle((3600, 500, 3800, 600), fill='red')
        DRAW.text((3900, 500), 'IRES score : {}'.format(ires_score), fill='black', font = arial100)
        DRAW.rectangle((3600, 700, 3800, 800), fill='orange')
        DRAW.text((3900, 700), 'ORF score : {}'.format(orf_score), fill='black', font = arial100)
        DRAW.text((3900, 900), 'ORF start site : {}'.format(orf_start_site), fill='black', font = arial100)
        DRAW.text((3900, 1100), 'ORF stop site : {}'.format(orf_stop_site), fill='black', font = arial100)

        DRAW.text((1900,2200), 'circRNA_lenth = {}'.format(circ_lenth), fill='black', font = arial100)
        DRAW.text((1900,2400), 'ORF_lenth = {}'.format(orf_lenth), fill='black', font = arial100)
        DRAW.text((1900,2600), 'transcript_id = {}'.format(transcript_id), fill='black', font = arial100)

        arial150 = ImageFont.truetype('LiberationSans-Regular.ttf', 150)
        DRAW.text((2280,1800), "3'", fill='black', font = arial150)
        DRAW.text((2600,1800), "5'", fill='black', font = arial150)

        new_name = circ_name + '|' + seq_record.id.split('|')[4]
        myseq.save(final_file_path + '/CIRC_graph/{}.png'.format(new_name))

    han_not = open(tmp_file_path + '/' + tmp_file_name + '_not_in_gtf_list','w+')
    han_not.writelines(not_in_gtf_list)
    han_not.close()


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

def main():

    begin = time.perf_counter()
    lenth_need, circrnas_file, circrna_gtf, tmp_file_name, tmp_file_path, final_name, final_file_path = file_load()

    time0 = time.perf_counter()
    draw_IRES_secondary_structure(tmp_file_name, tmp_file_path, final_name, final_file_path)
    time_add0 = time.perf_counter()-time0
    print('draw_IRES_secondary_structure is finished! Spending time is {:.2f}s'.format(time_add0))

    time1 = time.perf_counter()
    make_index(circrna_gtf, tmp_file_name, final_name, tmp_file_path, final_file_path)
    time_add1 = time.perf_counter()-time1
    print('make_index is finished! Spending time is {:.2f}s'.format(time_add1))

    time2 = time.perf_counter()
    make_circ_pic(lenth_need, tmp_file_name, tmp_file_path, final_name, final_file_path)
    time_add2 = time.perf_counter()-time2
    print('make_circ_pic is finished! Spending time is {:.2f}s'.format(time_add2))

    time3 = time.perf_counter()
    UpsetR(final_name, tmp_file_name, tmp_file_path, final_file_path)
    time_add3 = time.perf_counter()-time3
    print('UpsetR finished! Spending time is {:.2f}'.format(time_add3))

    print('the whole visual_circ_orf is finished! Spending time is {:.2f}s'.format(int(time.perf_counter() - begin)))
    print('End')


if __name__ == '__main__':
    main()







