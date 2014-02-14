#!/usr/bin/env python

#modify the junction index file and modify the junction length

import os
import sys
import copy
import optparse

# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-r", "--repeat-region", action = "store", type = "string", dest = "repeat_region",
		  help = "repeat region in transcripts")
parser.add_option("-j", "--junction-index", action = "store", type = "string", dest = "junction_index",
		  help = "junction index by chromosome")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output",
		  help = "output: junction index no repeat region by chromosome")

(options, args) = parser.parse_args()

if (options.repeat_region is None or
    options.junction_index is None or
    options.output is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_repeat = options.repeat_region
fname_index = options.junction_index
output1 = options.output

#fname_repeat='c:/Python26/test file/hg19_nonrepeat/RepeatRegionTrans.txt'
#fname_index='C:/python26/test file/hg19_nonrepeat/GencodeJunctionIndex_chr1.txt'
#output1='JunctionIndexNonRepeatRegion_chr1.txt'

#fname_repeat=sys.argv[1]
#fname_index=sys.argv[2]
#output1=sys.argv[3]

try:
    bk_repeat=open(fname_repeat)
except:
    print prog_base + ": error: cannot open file " + fname_repeat
    sys.exit(1)
sh_repeat=bk_repeat.readlines()
row_number_repeat=len(sh_repeat)
bk_repeat.close()

repeat_info={}

for i in range(1,len(sh_repeat)):
    temp1=sh_repeat[i][0:-1].split('\t')
    isoform_name=temp1[0]
    start=[]
    end=[]
    for j in range(1,len(temp1)-1):
        temp2=temp1[j][1:-1].split(',')
        start.append(int(temp2[0]))
        end.append(int(temp2[1]))
    repeat_info[isoform_name]=zip(start,end)
    
try:
    bk_index=open(fname_index)
except:
    print prog_base + ": error: cannot open file " + fname_index
    sys.exit(1)
sh_index=bk_index.readlines()
row_number_junction=len(sh_index)
bk_index.close()

junction_information={}
isoform_gene_pair={}
gene_strand={}
junction_index_modi={}

for w in range(1,len(sh_index)):###build dictory for junctions for isoform
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
        temp_str=gene_name+':'+str(temp6[1])+':'+temp1[0]# gene:junction_serial_number:isoform
    junction_information[temp1[0]]=temp5
    junction_index_modi[temp1[0]]=temp5
    isoform_gene_pair[temp1[0]]=gene_name

for v in junction_information.keys():
    repeat=repeat_info.get(v,'None')
    if repeat=='None':
        continue
    junction_info=junction_information[v]
    junction_modi=copy.copy(junction_info)
    for i in range(0,len(repeat)):
        for p in range(0,len(junction_info)):
            if int(repeat[i][0])<=int(junction_info[p][0]) and int(repeat[i][1])>int(junction_info[p][0]):
                junction_modi[p]=[]
        junction_modi=filter(lambda x:x !=[],junction_modi)    
        junction_info=copy.copy(junction_modi)
    junction_index_modi[v]=junction_modi

                        
try:
    file_output=open(output1, 'w')
except:
    print prog_base + ": error: cannot open file " + output1
    sys.exit(1)

titles=['IsoformID'+'\t'+'GeneID'+'\t'+'StrandDirection'+'\t'+'junctionNumber'+'\t'+'junctionLocation(StartPointInIsoform,junctionSerialNumber,StartPointInChromo,EndPointInChromo)'+'\n']
file_output.writelines(titles)

for v in junction_index_modi.keys():
    temp=junction_index_modi[v]
    if temp!=[]:
        insert_line=str(v)+'\t'+str(isoform_gene_pair[v])+'\t'+gene_strand[isoform_gene_pair[v]]+'\t'+str(len(temp))+'\t'+str(temp)+'\n'
        file_output.writelines(insert_line)

file_output.close( )

print 'junctionIndexModifyForChromosomeDone'

