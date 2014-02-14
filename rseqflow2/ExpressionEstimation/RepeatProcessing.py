#!/usr/bin/env python

#merge the repeat region in transcripts
#and then transfer into location in Chromosome

import os
import sys
import optparse

# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-r", "--repeat-region", action = "store", type = "string", dest = "repeat_region",
                  help = "repeat region by chromosome")
parser.add_option("-e", "--exon-index", action = "store", type = "string", dest = "exon_index",
                  help = "exon index by chromosome")
parser.add_option("-d", "--detect", action = "store", type = "int", dest = "detect",
                  help = "repeat detect length")
parser.add_option("-t", "--output-transcript", action = "store", type = "string", dest = "output_transcript",
                  help = "output: repeat region in transcripts")
parser.add_option("-c", "--output-chromosome", action = "store", type = "string", dest = "output_chromosome",
                  help = "output: repeat region in chromosome")

(options, args) = parser.parse_args()

if (options.repeat_region is None or
    options.exon_index is None or
    options.detect is None or
    options.output_transcript is None or
    options.output_chromosome is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_repeat = options.repeat_region
fname_index = options.exon_index
repeat_detect_length = options.detect
output1 = options.output_transcript
output2 = options.output_chromosome

#fname_repeat='c:/Python26/test file/hg19_nonrepeat/Sample_RepeatRegionGenHg19_noMerge.txt'
#fname_index='C:/python26/test file/hg19_nonrepeat/GencodeExonIndex_chr1.txt'
#repeat_detect_length=75
#output1='RepeatRegionTrans.txt'
#output2='RepeatRegionChrs.txt'

#fname_repeat=sys.argv[1]
#fname_index=sys.argv[2]
#repeat_detect_length=sys.argv[3]
#output1=sys.argv[4]
#output2=sys.argv[5]

try:
    bk_repeat=open(fname_repeat)
except:
    print prog_base + ": error: cannot open file " + fname_repeat
    sys.exit(1)
sh_repeat=bk_repeat.readlines()
row_number_repeat=len(sh_repeat)
bk_repeat.close()


try:
    bk_index=open(fname_index)
except:
    print prog_base + ": error: cannot open file " + fname_index
    sys.exit(1)
sh_index=bk_index.readlines()
row_number_exon=len(sh_index)
bk_index.close()

GeneIsoform={}
for i in range(1,row_number_exon):
    temp=sh_index[i].split('\t')
    GeneIsoform[temp[0]]=temp[1]


isoform_info={}
isoform_mapping=[]
isoform_gene={}
gene_repeat_base={}
for i in sh_repeat:
    temp1=i[0:-1].split('_',2)
    temp2=temp1[2].split('=')
    isoform=temp2[0]
    if isoform not in GeneIsoform:
	continue
    gene_name=GeneIsoform[isoform]
    temp3=temp2[1].split(':')
    chr=temp3[0]
    mapping=int(temp3[2])
    isoform_info[isoform]=[]
    isoform_mapping.append([isoform,mapping])
    isoform_gene[isoform]=chr+':'+gene_name
    gene_repeat_base[gene_name]=[]


for i in range(0,len(isoform_mapping)):
    isoform=isoform_mapping[i][0]
    temp1=isoform_info[isoform]
    temp1.append(isoform_mapping[i][1])
    isoform_info[isoform]=temp1

for i in isoform_info.keys():
    temp1=isoform_info[i]
    temp1.sort()
    temp2=[]
    temp3=[[temp1[0],temp1[0]+int(repeat_detect_length)-1]]
    for j in range(1,len(temp1)):
        if temp1[j]<=temp1[j-1]+int(repeat_detect_length):
            temp3[len(temp3)-1][1]=temp1[j]+int(repeat_detect_length)-1
        else:
            temp3.append([temp1[j],temp1[j]+int(repeat_detect_length)-1])
    isoform_info[i]=temp3

try:
    file_output=open(output1, 'w')
except:
    print prog_base + ": error: cannot open file " + output1
    sys.exit(1)

titles='IsoformID'+'\t'+'RepetitionRange'+'\t'+'GeneName'+'\n'
file_output.writelines(titles)

for v in isoform_info.keys():
    insert_line=str(v)
    info=isoform_info[v]
    for j in range(0,len(info)):
        insert_line=insert_line+'\t'+str(info[j])
    insert_line=insert_line+'\t'+str(isoform_gene[v])+'\n'
    file_output.writelines(insert_line)

file_output.close( )

print 'RepeationProcessingForTranscriptsDone'

##merge the repeat region of transcripts
###################################################################################################################


#############################################################################################################
###transfer the repeat in transcripts into chromosome
exon_information={}
isoform_gene_pair={}
gene_strand={}

for w in range(1,len(sh_index)):###build dictory for exons for isoform
    temp1=sh_index[w][0:-1].split('\t')
    gene_name=temp1[1]
    strand=temp1[2]
    gene_strand[gene_name]=strand
    temp4=str(temp1[4])[1:-2].split('), ')
    temp5=[]
    exon_info=[]
    for vv in range(0,len(temp4)):
        temp4[vv]=temp4[vv][1:]
        temp6=temp4[vv].split(',')
        temp5=temp5+temp6
    exon_info=exon_info+temp5
    exon_information[temp1[0]]=exon_info
    isoform_gene_pair[temp1[0]]=gene_name

for v in isoform_info.keys(): #just look repeat region in isoform as reads,mapped to exon location
    gene_name=isoform_gene_pair.get(v,'none')
    exon_info=exon_information.get(v,'None')
    if exon_info=='None':
        continue
    repeat=isoform_info[v]
    for i in range(0,len(repeat)):
        for p in range(0,(len(exon_info)/4)):
            strand=gene_strand[gene_name]
            if strand=='+':
                if p!=(len(exon_info)/4-1): #not the last exon
                    if int(repeat[i][0])>=int(exon_info[4*p]) and int(repeat[i][1])<int(exon_info[4*(p+1)]):
                        temp22_start=int(exon_info[4*p+2])+1+int(repeat[i][0])-int(exon_info[4*p])
                        temp22_end=int(exon_info[4*p+2])+1+int(repeat[i][1])-int(exon_info[4*p])
                        temp11=gene_repeat_base[gene_name]
                        temp11.append(zip([temp22_start],[temp22_end]))
                        gene_repeat_base[gene_name]=temp11
                    elif int(repeat[i][0])>=int(exon_info[4*p]) and int(repeat[i][0])<int(exon_info[4*(p+1)]) and int(repeat[i][1])>=int(exon_info[4*(p+1)]):
                        temp22_start=int(exon_info[4*p+2])+1+int(repeat[i][0])-int(exon_info[4*p])
                        temp22_end=int(exon_info[4*p+3])
                        temp11=gene_repeat_base[gene_name]
                        temp11.append(zip([temp22_start],[temp22_end]))
                        gene_repeat_base[gene_name]=temp11
                    elif int(repeat[i][0])<=int(exon_info[4*p]) and int(repeat[i][1])>=int(exon_info[4*p]) and int(repeat[i][1])<int(exon_info[4*(p+1)]):
                        temp22_start=int(exon_info[4*p+2])+1
                        temp22_end=int(repeat[i][1])-int(exon_info[4*p])+int(exon_info[4*p+2])+1
                        temp11=gene_repeat_base[gene_name]
                        temp11.append(zip([temp22_start],[temp22_end]))
                        gene_repeat_base[gene_name]=temp11
                    elif int(repeat[i][0])<=int(exon_info[4*p]) and int(repeat[i][1])>=int(exon_info[4*(p+1)]):
                        temp22_start=int(exon_info[4*p+2])+1
                        temp22_end=int(exon_info[4*p+3])
                        temp11=gene_repeat_base[gene_name]
                        temp11.append(zip([temp22_start],[temp22_end]))
                        gene_repeat_base[gene_name]=temp11
                else:  #the last exon
                    if int(repeat[i][0])>=int(exon_info[4*p]):
                        temp22_start=int(exon_info[4*p+2])+1+int(repeat[i][0])-int(exon_info[4*p])
                        temp22_end=int(exon_info[4*p+2])+1+int(repeat[i][1])-int(exon_info[4*p])
                        temp11=gene_repeat_base[gene_name]
                        temp11.append(zip([temp22_start],[temp22_end]))
                        gene_repeat_base[gene_name]=temp11
                    elif int(repeat[i][0])<=int(exon_info[4*p]) and int(repeat[i][1])>=int(exon_info[4*p]):
                        temp22_start=int(exon_info[4*p+2])+1
                        temp22_end=int(repeat[i][1])-int(exon_info[4*p])+int(exon_info[4*p+2])+1
                        temp11=gene_repeat_base[gene_name]
                        temp11.append(zip([temp22_start],[temp22_end]))
                        gene_repeat_base[gene_name]=temp11
            elif strand=='-':
                if p!=(len(exon_info)/4-1): #not the last exon
                    if int(repeat[i][0])>=int(exon_info[4*p]) and int(repeat[i][1])<int(exon_info[4*(p+1)]):
                        temp22_end=int(exon_info[4*p+3])-(int(repeat[i][0])-int(exon_info[4*p]))
                        temp22_start=int(exon_info[4*p+3])-(int(repeat[i][1])-int(exon_info[4*p]))
                        temp11=gene_repeat_base[gene_name]
                        temp11.append(zip([temp22_start],[temp22_end]))
                        gene_repeat_base[gene_name]=temp11
                    elif int(repeat[i][0])>=int(exon_info[4*p]) and int(repeat[i][0])<int(exon_info[4*(p+1)]) and int(repeat[i][1])>=int(exon_info[4*(p+1)]):
                        temp22_start=int(exon_info[4*p+2])+1
                        temp22_end=int(exon_info[4*p+3])-(int(repeat[i][0])-int(exon_info[4*p]))
                        temp11=gene_repeat_base[gene_name]
                        temp11.append(zip([temp22_start],[temp22_end]))
                        gene_repeat_base[gene_name]=temp11

                    elif int(repeat[i][0])<=int(exon_info[4*p]) and int(repeat[i][1])>=int(exon_info[4*p]) and int(repeat[i][1])<int(exon_info[4*(p+1)]):
                        temp22_end=int(exon_info[4*p+3])
                        temp22_start=int(exon_info[4*p+3])-(int(repeat[i][1])-int(exon_info[4*p]))
                        temp11=gene_repeat_base[gene_name]
                        temp11.append(zip([temp22_start],[temp22_end]))
                        gene_repeat_base[gene_name]=temp11

                    elif int(repeat[i][0])<=int(exon_info[4*p]) and int(repeat[i][1])>=int(exon_info[4*(p+1)]):
                        temp22_start=int(exon_info[4*p+2])+1
                        temp22_end=int(exon_info[4*p+3])
                        temp11=gene_repeat_base[gene_name]
                        temp11.append(zip([temp22_start],[temp22_end]))
                        gene_repeat_base[gene_name]=temp11

                else:  #the last exon
                    if int(repeat[i][0])>=int(exon_info[4*p]):
                        temp22_end=int(exon_info[4*p+3])-(int(repeat[i][0])-int(exon_info[4*p]))
                        temp22_start=int(exon_info[4*p+3])-(int(repeat[i][1])-int(exon_info[4*p]))
                        temp11=gene_repeat_base[gene_name]
                        temp11.append(zip([temp22_start],[temp22_end]))
                        gene_repeat_base[gene_name]=temp11


                    elif int(repeat[i][0])<=int(exon_info[4*p]) and int(repeat[i][1])>=int(exon_info[4*p]):
                        temp22_end=int(exon_info[4*p+3])
                        temp22_start=int(exon_info[4*p+3])-(repeat[i][1]-int(exon_info[4*p]))
                        temp11=gene_repeat_base[gene_name]
                        temp11.append(zip([temp22_start],[temp22_end]))
                        gene_repeat_base[gene_name]=temp11
        
 
for i in gene_repeat_base.keys():
    combine=[]
    temp=gene_repeat_base[i]
    if temp!=[]:
        temp.sort()
        combine=[[temp[0][0][0],temp[0][0][1]]]

    for j in range(1,len(temp)):
        if temp[j][0][0]<=temp[j-1][0][1] and temp[j][0][1]>temp[j-1][0][1] and temp[j][0][1]>combine[len(combine)-1][1]:
            combine[len(combine)-1][1]=temp[j][0][1]
        elif temp[j][0][0]>temp[j-1][0][1] and temp[j][0][0]>combine[len(combine)-1][1]:
            combine.append([temp[j][0][0],temp[j][0][1]])
    gene_repeat_base[i]=combine

try:
    file_output=open(output2, 'w')
except:
    print prog_base + ": error: cannot open file " + output2
    sys.exit(1)

titles='GeneID'+'\t'+'RepeatRegion'+'\n'
file_output.writelines(titles)

for v in gene_repeat_base.keys():
    temp=gene_repeat_base[v]
    if temp!=[]:
        insert_line=str(v)
        for j in range(0,len(temp)):
            insert_line=insert_line+'\t'+str(temp[j])
        insert_line=insert_line+'\n'
        file_output.writelines(insert_line)

file_output.close( )

print 'RepeationProcessingForChromosomeDone'
