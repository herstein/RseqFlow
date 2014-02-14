#!/usr/bin/python

'''
Created on 2009-10-30
Modified on 2012-07-18

@author: Administrator
@modify by Liu Lin
'''

import os
import sys
import optparse
import copy
# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-g", "--gene-annotation", action = "store", type = "string", dest = "gene_annotation",
		  help = "gene annotation from UCSC")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output",
		  help = "output: exon combination by chromosome")

(options, args) = parser.parse_args()

if (options.gene_annotation is None or
    options.output is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname = options.gene_annotation
outputFilename = options.output

#fname=sys.argv[1]
#outputFilename=sys.argv[2]

#fname='c:/Python26/test file/GencodeAnnotationManual_chrY.gtf'
#outputFilename='GencodeExonCombination_chrY.txt'

try:
    bk=open(fname)
except:
    print prog_base + ": error: cannot open file " + fname
    sys.exit(1)
sh=bk.readlines()
row_number=len(sh)
bk.close()

#JSH 2013-08-26
### create empty file if chr not present (required for PEGASUS workflow) ###
if(row_number==0):
    print "ExonCombination.py: No data for " + fname
    try:
        open(outputFilename,'w').close()
    except:
        print prog_base + ": error: cannot create file " + outputFilename
        sys.exit(1)

    sys.exit(0)

 
genename_repeat=[]
genename=[]
whole_info=[]

##############################################################################
#-------------------get gene list in chromosome------------------------------#

#************************initial*******************#
pre_transcript_id	='none'
pre_gene_id='none'
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
    if transcript_id!=pre_transcript_id:
        if v!=0:
            genename_repeat.append(pre_gene_id)
            #transcripts_list.append(pre_transcript_id)
            whole_info.append([pre_transcript_id,pre_strand,exonCount,exonStarts,exonEnds,pre_gene_id])
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
whole_info.append([transcript_id,pre_strand,exonCount,exonStarts,exonEnds,pre_gene_id])

genename=list(set(genename_repeat))
########get gene list in chromosome        
###################################################################################        

gene_info=[]# gene and isoform information
exon_info=[]#all possible exon in one gene 
exon_number_include_within=[]     
    
for v in range(0,len(genename)):
    genename1=genename[v]
    exon_start=[]
    exon_end=[]
    whole_length_allexons=0
#################################################
#find gene isoform and put exon number,starts,ends to list
    for number1 in range(0,len(whole_info)):
        if whole_info[number1][5]==genename1:
            gene_name=genename1	
            strand=whole_info[number1][1]
            exon_number=whole_info[number1][2]
            exon_start_temp1=whole_info[number1][3].split(',')
            exon_end_temp1=whole_info[number1][4].split(',')
            exon_start_temp2=[]
            exon_end_temp2=[]
            for number2 in range(0,int(exon_number)):
                exon_start_temp2.append(int(exon_start_temp1[number2]))
                exon_end_temp2.append(int(exon_end_temp1[number2]))
            exon_start=exon_start+exon_start_temp2
            exon_end=exon_end+exon_end_temp2
#find gene isoform and put exon number,starts,ends to list
####################################################

######################################################################
###########combine the isoforms information together to get all exons
    
    exon_start_end_repeat_disorder=list(set(zip(exon_start,exon_end)))
    exon_start_end_repeat_disorder.sort()
    exon_start_end_repeat_delete=exon_start_end_repeat_disorder
    exon_number_include_within.append([genename1,strand,len(exon_start_end_repeat_delete),exon_start_end_repeat_delete])
    wait_delete_temp=[] #find the exons within other exons and record it
    for numi in range(0,len(exon_start_end_repeat_delete)):
        for numj in range(numi+1,len(exon_start_end_repeat_delete)):
            if (exon_start_end_repeat_delete[numi][0]<=exon_start_end_repeat_delete[numj][0])&(exon_start_end_repeat_delete[numi][1]>=exon_start_end_repeat_delete[numj][1]):
                wait_delete_temp.append(numj)
            if (exon_start_end_repeat_delete[numi][0]==exon_start_end_repeat_delete[numj][0])&(exon_start_end_repeat_delete[numi][1]<=exon_start_end_repeat_delete[numj][1]):
                wait_delete_temp.append(numi)


    wait_delete=list(set(wait_delete_temp))#get the un-repeated serial number of exons_to_be_delete
    wait_delete.sort()
    wait_delete.reverse()
    exon_start_end=copy.copy(exon_start_end_repeat_delete) #delete the within exons     
    for p in wait_delete:
        del exon_start_end[p]

    Length=exon_start_end[0][1]-exon_start_end[0][0]
    for numi in range(1,len(exon_start_end)):
        if exon_start_end[numi][0]<exon_start_end[numi-1][1] and exon_start_end[numi][1]>=exon_start_end[numi-1][1]:
             Length=Length+(exon_start_end[numi][1]-exon_start_end[numi-1][1])
        elif exon_start_end[numi][0]>=exon_start_end[numi-1][1]:
             Length=Length+(exon_start_end[numi][1]-exon_start_end[numi][0])
    exon_info.append([genename1]+[Length])

#############################################################################
######write to file

try:
    ws=open(outputFilename,'w')
except:
    print prog_base + ": error: cannot open file " + outputFilename
    sys.exit(1)

titles=['GeneName'+'\t'+'StrandDirection'+'\t'+'GeneMappableLength'+'\t'+'ExonsNumber'+'\t'+'(ExonStart,ExonEnd)'+'\n']
ws.writelines(titles)

for i in range(len(exon_number_include_within)):
    insertline=str(exon_number_include_within[i][0])+'\t'+str(exon_number_include_within[i][1])+'\t'+str(exon_info[i][1])+'\t'+str(exon_number_include_within[i][2])+'\t'
    for j in range(3,len(exon_number_include_within[i])):
        insertline=insertline+str(exon_number_include_within[i][j])
    ws.writelines(insertline+'\n')

ws.close()

print 'ExonCombination Done for Chromosome %s' % (chr)
