#!/usr/bin/env python

'''
Created on 2009-10-30
Modified on 2012-07-19

@author: Administrator
@modify by Liu Lin
'''

#####################################
#this program will get the unrepeated junctions list for each gene 

import os
import sys
import optparse

# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-g", "--gene-annotation", action = "store", type = "string", dest = "gene_annotation",
		  help = "gene annotation from UCSC")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output",
		  help = "output: junction combination by chromosome")

(options, args) = parser.parse_args()

if (options.gene_annotation is None or
    options.output is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

file_exon_raw = options.gene_annotation
output_filename = options.output

#file_exon_raw='c:/Python26/test file/GencodeAnnotationManual_chrY.gtf'
#output_filename='GencodeJunctionCombination_chrY.txt'

#file_exon_raw=sys.argv[1]
#output_filename=sys.argv[2]

try:
    bk=open(file_exon_raw)
except:
    print prog_base + ": error: cannot open file " + file_exon_raw
    sys.exit(1)
sh=bk.readlines()
row_number=len(sh)
bk.close()

#JSH 2013-08-26
### create empty file if chr not present (required for PEGASUS workflow) ###
if(row_number==0):
    print "JunctionCombination.py: No data for " + file_exon_raw
    try:
        open(output_filename,'w').close()
    except:
        print prog_base + ": error: cannot create file " + output_filename
        sys.exit(1)

    sys.exit(0)


genename_repeat=[]
genename=[]
gene_info=[]


##############################################################################
#######get gene list in chromosome
pre_transcript_id	='none'
chr = ''
for v in range(0,row_number):
    temp=sh[v][0:-1].split('\t')
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
    #    print "Gene ID:"+gene_id
    #    print "Transcript ID:"+transcript_id+
    #    print "Gene ID and Transcript ID are the same, so please check the annotaion file." 
        #sys.exit(1)
    if transcript_id!=pre_transcript_id:
        if v!=0:
            genename_repeat.append(pre_gene_id)
            #transcripts_list.append(pre_transcript_id)
            gene_info.append([pre_transcript_id,pre_strand,exonCount,exonStarts,exonEnds,pre_gene_id])
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
genename_repeat.append(pre_gene_id)
#transcripts_list.append(pre_transcript_id)
gene_info.append([transcript_id,pre_strand,exonCount,exonStarts,exonEnds,pre_gene_id])

genename=list(set(genename_repeat))
 ########get gene list in chromosome        
###################################################################################        


junction_info=[]#all possible junctions in one gene 
    
for v in range(0,len(genename)):
    genename1=genename[v]
    junction_start=[]
    junction_end=[]
    

#################################################
#find gene isoform and put junction number,starts,ends points to list
    for number1 in range(0,len(gene_info)):
        if gene_info[number1][5]==genename1:
            gene_name=genename1
            junction_number=int(gene_info[number1][2])-1
            exon_start_temp1=gene_info[number1][3].split(',')
            exon_end_temp1=gene_info[number1][4].split(',')
            strand=gene_info[number1][1]
            junction_start_temp2=[]
            junction_end_temp2=[]
            for number2 in range(0,junction_number):
                junction_start_temp2.append(int(exon_end_temp1[number2]))
                junction_end_temp2.append(int(exon_start_temp1[number2+1]))
            junction_start=junction_start+junction_start_temp2
            junction_end=junction_end+junction_end_temp2
            

######################################################################
###########combine the isoforms information together to get all junctions
    
    junction_start_end_repeat_delete=list(set(zip(junction_start,junction_end)))
    junction_info.append([genename1]+[strand]+[len(junction_start_end_repeat_delete)]+junction_start_end_repeat_delete)
    
#############################################################################
######write the results to file

try:
    ws=open(output_filename,'w')
except:
    print prog_base + ": error: cannot open file " + output_filename
    sys.exit(1)

titles=['GeneID'+'\t'+'strand'+'\t'+'JunctionNumber'+'\t'+'JunctionStart,JunctionEnd'+'\n']
ws.writelines(titles)
	
for i in range(len(junction_info)):
    insertline=str(junction_info[i][0])+'\t'+str(junction_info[i][1])+'\t'+str(junction_info[i][2])+'\t'
    for j in range(3,len(junction_info[i])):
        insertline=insertline+str(junction_info[i][j])
    ws.writelines(insertline+'\n')

ws.close()

print 'JunctionCombination Done for Chromosome %s' % (chr)

