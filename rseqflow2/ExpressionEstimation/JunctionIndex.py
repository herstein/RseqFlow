#!/usr/bin/python

'''
Created on 2009-10-30
Modified on 2012-07-19

@author: Administrator
@modify by Liu Lin
'''

import os
import sys
import optparse

# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-g", "--gene-annotation", action = "store", type = "string", dest = "gene_annotation",
		  help = "gene annotation from UCSC")
parser.add_option("-j", "--junction-combination", action = "store", type = "string", dest = "junction_combination",
		  help = "junction combination by chromosome")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output",
		  help = "output: junction index by chromosome")

(options, args) = parser.parse_args()

if (options.gene_annotation is None or
    options.junction_combination is None or
    options.output is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_gene_isoform_pair = options.gene_annotation
fname_junction = options.junction_combination
output_filename = options.output

#####################################
###this program will build the table of junction information list
###including isoform_id ,gene_id,absolute_location_isoform 

#fname_gene_isoform_pair='C:/python26/test file/GencodeAnnotationManual_chrY.gtf'
#fname_junction='C:/python26/test file/GencodeJunctionCombination_chrY.txt'
#output_filename='GencodeJunctionIndex_chrY.txt'

#fname_gene_isoform_pair=sys.argv[1]
#fname_junction=sys.argv[2]
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
    print "JunctionIndex.py: No data for " + fname_gene_isoform_pair
    try:
        open(output_filename,'w').close()
    except:
        print prog_base + ": error: cannot create file " + output_filename
        sys.exit(1)

    sys.exit(0)


try:
    bk_junc=open(fname_junction)
except:
    print prog_base + ": error: cannot open file " + fname_junction
    sys.exit(1)
sh_junc=bk_junc.readlines()
row_number_junc=len(sh_junc)
bk_junc.close()

#JSH 2013-08-26
### create empty file if chr not present (required for PEGASUS workflow) ###
if(row_number_junc==0):
    print "JunctionIndex.py: No data for " + fname_junction
    try:
        open(output_filename,'w').close()
    except:
        print prog_base + ": error: cannot create file " + output_filename
        sys.exit(1)

    sys.exit(0)


junction_combination_sort=[]
original_form_sort=[]
junction_info_list=[]
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
    #if gene_id==transcript_id:
    #    print "Gene ID:"+gene_id+'\n'
    #    print "Transcrit ID:"+transcript_id+'\n'
    #    print "Warning:gene ID and Transcript ID are the same, so please check the annotaion file." 
        #ys.exit(1)
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

for v in range(1,row_number_junc): 
    temp3=sh_junc[v].split('\t')
    if len(temp3)==4:   #there is junction for this isoform
        temp3[3]=temp3[3][0:-1]
    elif len(temp3)==3: #there is no junction for this isoform
        temp3[2]=temp3[2][0:-1]
    junction_combination_sort.append(temp3)
junction_combination_sort.sort()

for v in range(0,len(original_form_sort)):    
    gene_name_temp1=original_form_sort[v][0]
    gene_name_list.append(gene_name_temp1)
    isoform_name=original_form_sort[v][1]
    junction_number=int(original_form_sort[v][3])-1
    strand=original_form_sort[v][2]
    exon_start_temp1=original_form_sort[v][4].split(',')
    exon_end_temp1=original_form_sort[v][5].split(',')
    
    del exon_end_temp1[len(exon_end_temp1)-1] #delete the last empty element
    del exon_start_temp1[len(exon_start_temp1)-1]

    junction_start_temp3=[]
    junction_end_temp3=[]
    location_isoform=0
    junction_location_isoform=[]
    serial_number_temp=[]
    
    ##########################################################################################
    ####read out the current gene's all junction information and then delete from the database
    if v==0 or gene_name_list[v]!=gene_name_list[v-1]:
        current_gene_info=[]       
        for number1 in range(0,len(junction_combination_sort)):
            findit='False'
            if junction_combination_sort[number1][0]==gene_name_temp1:
                findit='True'
                temp4=junction_combination_sort[number1]
                temp5=temp4[3].split(')')
                temp7=[]
                for i in range(len(temp5)-1):
                    temp6=temp5[i].split(',')
                    temp6[0]=temp6[0][1:]
                    temp6[1]=temp6[1][1:]
                    temp7=temp7+temp6
                current_gene_info=[gene_name_temp1]+temp7
                break #if find the gene's junction information ,then skip out
        if findit=='True':
            del junction_combination_sort[number1]
    ############################################################################################  
    if junction_number==0:
        junction_info_list.append([isoform_name]+[gene_name_temp1]+[0,0])
        continue     
    ###################
    if strand=='+':
        for number2 in range(0,junction_number):
            junction_start_temp3=int(exon_end_temp1[number2])
            junction_end_temp3=int(exon_start_temp1[number2+1])
            for number1 in range(0,(len(current_gene_info)-1)/2):  
                if junction_start_temp3==int(current_gene_info[2*number1+1]) and junction_end_temp3==int(current_gene_info[2*number1+2]):
                    serial_number_temp=number1+1
                    break
            location_isoform=location_isoform+(int(exon_end_temp1[number2]))-(int(exon_start_temp1[number2]))
            junction_location_isoform=junction_location_isoform+zip([location_isoform],[serial_number_temp],[int(exon_end_temp1[number2])],[int(exon_start_temp1[number2+1])])
    if strand=='-':
        x=range(1,junction_number+1);
        x.sort(reverse=1)
        for number2 in x:
            junction_start_temp3=int(exon_end_temp1[number2-1])
            junction_end_temp3=int(exon_start_temp1[number2])
            for number1 in range(0,(len(current_gene_info)-1)/2):  
                if junction_start_temp3==int(current_gene_info[2*number1+1]) and junction_end_temp3==int(current_gene_info[2*number1+2]):
                    serial_number_temp=number1+1
                    break
            location_isoform=location_isoform+(int(exon_end_temp1[number2]))-(int(exon_start_temp1[number2]))
            junction_location_isoform=junction_location_isoform+zip([location_isoform],[serial_number_temp],[int(exon_end_temp1[number2-1])],[int(exon_start_temp1[number2])])
        
    junction_info_list.append([isoform_name]+[gene_name_temp1]+[strand]+[junction_number]+junction_location_isoform)
#############################################################################
######write the results to file

try:
    ws=open(output_filename,'w')
except:
    print prog_base + ": error: cannot open file " + output_filename
    sys.exit(1)

titles=['IsoformID'+'\t'+'GeneID'+'\t'+'StrandDirection'+'\t'+'JunctionNumber'+'\t'+'JunctionLocation(StartPointInIsoform,JunctionSerialNumber,StartPointInChromosome,EndPointInChromosome)'+'\n']
ws.writelines(titles)

for i in range(0,len(junction_info_list)):
    if junction_info_list[i][2]!=0:
        ws.writelines(junction_info_list[i][0]+'\t'+junction_info_list[i][1]+'\t'+str(junction_info_list[i][2])+'\t'+str(junction_info_list[i][3])+'\t'+str(junction_info_list[i][4:])+'\n')

ws.close()

print 'JunctionIndex Done for Chromosome %s' % (chr)

