#!/usr/bin/python

'''
Created on 2009-11-25

@author: Administrator
'''

import os
import sys
import optparse

# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-g", "--gene-annotation", action = "store", type = "string", dest = "gene_annotation",
                  help = "gene annotation from UCSC")
parser.add_option("-e", "--exon-combination", action = "store", type = "string", dest = "exon_combination",
                  help = "exon combination by chromosome")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output",
                  help = "output: exon index by chromosome")

(options, args) = parser.parse_args()

if (options.gene_annotation is None or
    options.exon_combination is None or
    options.output is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_gene_isoform_pair = options.gene_annotation
fname_exon = options.exon_combination
output_filename = options.output

#####################################
###this program will build the table of non-repetive exons list for each gene
#fname_gene_isoform_pair='C:/python26/test file/GencodeAnnotationManual_chrY.gtf'
#fname_exon='C:/python26/test file/GencodeExonCombination_chrY.txt'
#output_filename='GencodeExonIndex_chrY.txt'

#fname_gene_isoform_pair=sys.argv[1]
#fname_exon=sys.argv[2]
#output_filename=sys.argv[3]

try:
    bk_pair=open(fname_gene_isoform_pair)
except:
    print prog_base + ": error: cannot open file " + fname_gene_isoform_pair
    sys.exit(1)
sh_pair=bk_pair.readlines()
row_number_pair=len(sh_pair)
bk_pair.close()

#JSH 2013-08-26
### create empty file if chr not present (required for PEGASUS workflow) ###
if(row_number_pair==0):
    print "ExonIndex.py: No chr data for " + fname_gene_isoform_pair
    try:
        open(output_filename,'w').close()
    except:
        print prog_base + ": error: cannot create file " + output_filename
        sys.exit(1)

    sys.exit(0)


try:
    bk_exon=open(fname_exon)
except:
    print prog_base + ": error: cannot open file " + fname_exon
    sys.exit(1)
sh_exon=bk_exon.readlines()
row_number_exon=len(sh_exon)
bk_exon.close()

#JSH 2013-08-26
### create empty file if chr not present (required for PEGASUS workflow) ###
if(row_number_exon==0):
    print "ExonIndex.py: No chr data for " + fname_exon
    try:
        open(output_filename,'w').close()
    except:
        print prog_base + ": error: cannot create file " + output_filename
        sys.exit(1)

    sys.exit(0)


exon_combination_sort=[]
original_form_sort=[]
exon_info_list=[]
current_gene_info=[]
gene_name_list=[]


##############################################################################
#######get gene list in chromosome
pre_transcript_id	='none'
chr = ''
for v in range(0,row_number_pair):
    temp=sh_pair[v][0:-1].split('\t')
    chr			=temp[0]
    type		=temp[2]
    start		=int(temp[3])-1
    end			=int(temp[4])
    strand	=temp[6]
    temp1		=temp[8]
    attributes		=temp1.split(';')
    temp2					=attributes[0].split('"')
    gene_id				=temp2[1]
    temp3					=attributes[1].split('"')
    transcript_id	=temp3[1]
    if transcript_id!=pre_transcript_id:
        if v!=0:
            original_form_sort.append([pre_gene_id,pre_transcript_id,pre_strand,exonCount,exonStarts,exonEnds])
        exonCount	=0
        exonStarts=''
        exonEnds  =''
    if type=='exon':
        exonCount	=exonCount+1
        exonStarts=exonStarts+str(start)+','
        exonEnds	=exonEnds+str(end)+','
    pre_gene_id=gene_id
    pre_transcript_id=transcript_id
    pre_strand=strand    
    #whole_info.append([temp1[1],temp1[3],temp1[8],temp1[9],temp1[10],temp1[12]])
#**************************save the last transcript information**************    
original_form_sort.append([pre_gene_id,transcript_id,pre_strand,exonCount,exonStarts,exonEnds])

original_form_sort.sort()

for v in range(1,row_number_exon):
    temp3=sh_exon[v].split('\t')
    temp3[4]=temp3[4][0:-1]
    exon_combination_sort.append(temp3)
exon_combination_sort.sort()

    
for v in range(0,len(original_form_sort)):    
    gene_name_temp1=original_form_sort[v][0]
    gene_name_list.append(gene_name_temp1)
    isoform_name=original_form_sort[v][1]
    exon_number=int(original_form_sort[v][3])
    strand=original_form_sort[v][2]
    exon_start_temp1=original_form_sort[v][4].split(',')
    exon_end_temp1=original_form_sort[v][5].split(',')

    del exon_end_temp1[len(exon_end_temp1)-1] #delete the last empty element
    del exon_start_temp1[len(exon_start_temp1)-1]
    

    exon_start_temp3=[]
    exon_end_temp3=[]
    location_isoform=0
    exon_location_isoform=[]
    serial_number_temp=[]
    
    ##########################################################################################
    ####read out the current gene's all exon information and then delete from the database
    if v==0 or gene_name_list[v]!=gene_name_list[v-1]:# a new gene appear
        current_gene_info=[]       
        for number1 in range(0,len(exon_combination_sort)):
            findit='False'
            if exon_combination_sort[number1][0]==gene_name_temp1:
                findit='True'
                temp4=exon_combination_sort[number1]
                temp5=temp4[4].split('),')
                temp7=[]
                for i in range(len(temp5)):
                    temp6=temp5[i].split(',')
                    temp6[0]=temp6[0][2:]
                    temp6[1]=temp6[1][1:]
                    if i==(len(temp5)-1):
                        temp6[1]=temp6[1][0:-2]
                    temp7=temp7+temp6
                current_gene_info=[gene_name_temp1]+temp7
                break #if find the gene's exon information ,then skip out
        if findit=='True':
            del exon_combination_sort[number1]
        
            
    if strand=='+':
        for number2 in range(0,exon_number):
            exon_start_temp3=int(exon_start_temp1[number2])
            exon_end_temp3=int(exon_end_temp1[number2])
            for number1 in range(0,(len(current_gene_info)-1)/2):
                if exon_start_temp3==int(current_gene_info[2*number1+1]) and exon_end_temp3==int(current_gene_info[2*number1+2]):
                    serial_number_temp=number1+1
                    break
            if number2==0:
                location_isoform=1
            else:
                location_isoform=location_isoform+(int(exon_end_temp1[number2-1]))-(int(exon_start_temp1[number2-1]))
            exon_location_isoform=exon_location_isoform+zip([location_isoform],[serial_number_temp],[int(exon_start_temp1[number2])],[int(exon_end_temp1[number2])])
    if strand=='-':
        x=range(0,exon_number)
        x.sort(reverse=1)
        for number2 in x:
            exon_start_temp3=int(exon_start_temp1[number2])
            exon_end_temp3=int(exon_end_temp1[number2])
            for number1 in range(0,(len(current_gene_info)-1)/2):
                if exon_start_temp3==int(current_gene_info[2*number1+1]) and exon_end_temp3==int(current_gene_info[2*number1+2]):
                    serial_number_temp=number1+1
                    break
        
            if number2==exon_number-1:
                location_isoform=1
            else:
                location_isoform=location_isoform+(int(exon_end_temp1[number2+1]))-(int(exon_start_temp1[number2+1]))
            exon_location_isoform=exon_location_isoform+zip([location_isoform],[serial_number_temp],[int(exon_start_temp1[number2])],[int(exon_end_temp1[number2])])
       
    exon_info_list.append([isoform_name]+[gene_name_temp1]+[strand]+[exon_number]+exon_location_isoform)

#############################################################################
######write the results to file
try:
    ws=open(output_filename,'w')
except:
    print prog_base + ": error: cannot open file " + output_filename
    sys.exit(1)

titles=['IsoformID'+'\t'+'GeneID'+'\t'+'StrandDirection'+'\t'+'ExonNumber'+'\t'+'ExonLocation(StartPointInIsoform,ExonSerialNumber,StartPointInChromo,EndPointInChromo)'+'\n']
ws.writelines(titles)

for i in range(0,len(exon_info_list)):
    ws.writelines(exon_info_list[i][0]+'\t'+exon_info_list[i][1]+'\t'+str(exon_info_list[i][2])+'\t'+str(exon_info_list[i][3])+'\t'+str(exon_info_list[i][4:])+'\n')

ws.close()

print 'ExonIndex Done for Chromosome %s' % (chr)
