#!/usr/bin/env python

#modify the exon index file and modify the exon length

import os
import sys
import copy
import optparse

# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-r", "--repeat-region", action = "store", type = "string", dest = "repeat_region",
                  help = "repeat region in chromosome")
parser.add_option("-e", "--exon-index", action = "store", type = "string", dest = "exon_index",
                  help = "exon index by chromosome")
parser.add_option("-o", "--output-index", action = "store", type = "string", dest = "output_index",
                  help = "output: exon index no repeat region by chromosome")
parser.add_option("-l", "--output-length", action = "store", type = "string", dest = "output_length",
                  help = "output: exon length no repeat region by chromosome")

(options, args) = parser.parse_args()

if (options.repeat_region is None or
    options.exon_index is None or
    options.output_index is None or
    options.output_length is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_repeat = options.repeat_region
fname_index = options.exon_index
output1 = options.output_index
output2 = options.output_length

#fname_repeat='c:/Python26/test file/hg19_nonrepeat/RepeatRegionChrs.txt'
#fname_index='C:/python26/test file/hg19_nonrepeat/GencodeExonIndex_chr1.txt'
#output1='ExonIndexNonRepeatRegion_chr1.txt'
#output2='ExonLengthNonRepeatRegion_chr1.txt'

#fname_repeat=sys.argv[1]
#fname_index=sys.argv[2]
#output1=sys.argv[3]
#output2=sys.argv[4]

try:
    bk_repeat=open(fname_repeat)
except:
    print prog_base + ": error: cannot open file " + fname_repeat
    sys.exit(1)
sh_repeat=bk_repeat.readlines()
row_number_repeat=len(sh_repeat)
bk_repeat.close()

repeat_info={}
isoform_mapping=[]
isoform_gene={}
gene_repeat_base={}

for i in range(1,len(sh_repeat)):
    temp1=sh_repeat[i][0:-1].split('\t')
    gene_name=temp1[0]
    start=[]
    end=[]
    for j in range(1,len(temp1)):
        temp2=temp1[j][1:-1].split(',')
        start.append(int(temp2[0]))
        end.append(int(temp2[1]))
    repeat_info[gene_name]=zip(start,end)
    
try:
    bk_index=open(fname_index)
except:
    print prog_base + ": error: cannot open file " + fname_index
    sys.exit(1)
sh_index=bk_index.readlines()
row_number_exon=len(sh_index)
bk_index.close()

exon_information={}
isoform_gene_pair={}
gene_strand={}
exon_index_modi={}
exon_length={}

for w in range(1,len(sh_index)):###build dictory for exons for isoform
    temp1=sh_index[w][0:-1].split('\t')
    gene_name=temp1[1]
    strand=temp1[2]
    gene_strand[gene_name]=strand
    temp4=str(temp1[4])[1:-2].split('), ')
    temp5=[]
    for vv in range(0,len(temp4)):
        temp4[vv]=temp4[vv][1:]
        temp6=temp4[vv].split(',')
        for i in range(0,len(temp6)):
            temp6[i]=int(temp6[i])
        temp5.append(temp6)
        temp_str=gene_name+':'+str(temp6[1])+':'+temp1[0]# gene:exon_serial_number:isoform
        exon_length[temp_str]=temp6[3]-temp6[2]
    exon_information[temp1[0]]=temp5
    exon_index_modi[temp1[0]]=temp5
    isoform_gene_pair[temp1[0]]=gene_name

for v in exon_information.keys():
    gene=isoform_gene_pair[v]
    repeat=repeat_info.get(gene,'None')
    if repeat=='None':
        continue
    strand=gene_strand[gene]
    exon_info=exon_information[v]
    exon_modi=copy.copy(exon_info)
    for i in range(0,len(repeat)):
        for p in range(0,len(exon_info)):
            temp_str=gene+':'+str(exon_info[p][1])+':'+v
            exon_len=exon_length[temp_str]
            if strand=='+':
                if int(repeat[i][0])>int(exon_info[p][2])+1 and int(repeat[i][1])<int(exon_info[p][3]):
                    new1=[exon_info[p][0],exon_info[p][1],exon_info[p][2]+1,repeat[i][0]-1]
                    new2=[int(exon_info[p][0])+1+repeat[i][1]-int(exon_info[p][2])-1,exon_info[p][1],repeat[i][1]+1,exon_info[p][3]]  
                    exon_modi[p]=new1
                    exon_modi.append(new2)
                    exon_len=exon_len-(repeat[i][1]-repeat[i][0]+1)
                elif int(repeat[i][0])>int(exon_info[p][2])+1 and int(repeat[i][0])<int(exon_info[p][3]) and int(repeat[i][1])>=int(exon_info[p][3]):
                    new=[exon_info[p][0],exon_info[p][1],exon_info[p][2]+1,repeat[i][0]-1]
                    exon_modi[p]=new
                    exon_len=exon_len-(exon_info[p][3]-repeat[i][0]+1)
                elif int(repeat[i][0])<=int(exon_info[p][2])+1 and int(repeat[i][1])>=int(exon_info[p][2])+1 and int(repeat[i][1])<int(exon_info[p][3]):
                    new=[int(exon_info[p][0])+1+repeat[i][1]-int(exon_info[p][2])-1,exon_info[p][1],repeat[i][1]+1,exon_info[p][3]]
                    exon_modi[p]=new
                    exon_len=exon_len-(repeat[i][1]-exon_info[p][2])
                elif int(repeat[i][0])<=int(exon_info[p][2])+1 and int(repeat[i][1])>=int(exon_info[p][3]):
                    exon_modi[p]=[]
                    exon_len=0

            elif strand=='-':
                if int(repeat[i][0])>int(exon_info[p][2])+1 and int(repeat[i][1])<int(exon_info[p][3]):
                    new1=[exon_info[p][0],exon_info[p][1],repeat[i][1]+1,exon_info[p][3]]
                    new2=[int(exon_info[p][0])+int(exon_info[p][3])-repeat[i][0]+1,exon_info[p][1],exon_info[p][2]+1,repeat[i][0]]####
                    exon_modi[p]=new1
                    exon_modi.append(new2)
                    exon_len=exon_len-(repeat[i][1]-repeat[i][0]+1)
                elif int(repeat[i][0])>int(exon_info[p][2])+1 and int(repeat[i][0])<int(exon_info[p][3]) and int(repeat[i][1])>=int(exon_info[p][3]):
                    new=[int(exon_info[p][0])+int(exon_info[p][3])-repeat[i][0]+1,exon_info[p][1],exon_info[p][2]+1,repeat[i][0]-1]
                    exon_modi[p]=new
                    exon_len=exon_len-(exon_info[p][3]-repeat[i][0]+1)

                elif int(repeat[i][0])<=int(exon_info[p][2])+1 and int(repeat[i][1])>=int(exon_info[p][2])+1 and int(repeat[i][1])<int(exon_info[p][3]):
                    new=[exon_info[p][0],exon_info[p][1],repeat[i][1]+1,exon_info[p][3]]
                    exon_modi[p]=new
                    exon_len=exon_len-(repeat[i][1]-exon_info[p][2])

                elif int(repeat[i][0])<=int(exon_info[p][2])+1 and int(repeat[i][1])>=int(exon_info[p][3]):
                    exon_modi[p]=[]
                    exon_len=0

            exon_length[temp_str]=exon_len
        exon_modi=filter(lambda x:x !=[],exon_modi)    
        exon_info=copy.copy(exon_modi)
    exon_index_modi[v]=exon_modi

                        
try:
    file_output=open(output1, 'w')
except:
    print prog_base + ": error: cannot open file " + output1
    sys.exit(1)

titles=['IsoformID'+'\t'+'GeneID'+'\t'+'StrandDirection'+'\t'+'ExonNumber'+'\t'+'ExonLocation(StartPointInIsoform,ExonSerialNumber,StartPointInChromo,EndPointInChromo)'+'\n']
file_output.writelines(titles)

for v in exon_index_modi.keys():
    temp=exon_index_modi[v]
    if temp!=[]:
        temp.sort()
        insert_line=str(v)+'\t'+str(isoform_gene_pair[v])+'\t'+gene_strand[isoform_gene_pair[v]]+'\t'+str(len(temp))+'\t'+str(temp)+'\n'
        file_output.writelines(insert_line)

file_output.close()

print 'ExonIndexModifyForChromosomeDone'

try:
    file_output=open(output2, 'w')
except:
    print prog_base + ": error: cannot open file " + output2
    sys.exit(1)

titles=['GeneID:ExonSerialNumber:Isoform'+'\t'+'Length'+'\n']
file_output.writelines(titles)



for v in exon_length.keys():
    if exon_length[v]!=0:
        insert_line=str(v)+'\t'+str(exon_length[v])+'\n'
        file_output.writelines(insert_line)

file_output.close()

print 'ExonLengthModifyForChromosomeDone'
