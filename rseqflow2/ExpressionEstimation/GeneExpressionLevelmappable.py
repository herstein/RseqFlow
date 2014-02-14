#!/usr/bin/env python

# this program will count the bases mapped to exons, and then count RPKM 
# to get the approximate exon's expression level 

import os
import sys
import optparse

# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-s", "--uniq-map", action = "store", type = "string", dest = "uniq_map",
                  help = "uniq map sam file")
parser.add_option("-m", "--mappable-exon-index", action = "store", type = "string", dest = "mappable_exon_index",
                  help = "mappable exon index no repeat regions")
parser.add_option("-r", "--repeat-region", action = "store", type = "string", dest = "repeat_region",
                  help = "repetition region in chromosome")
parser.add_option("-e", "--exon-combination", action = "store", type = "string", dest = "exon_combination",
                  help = "Exon combination produced before hand")
parser.add_option("-c", "--mapping-count", action = "store", type = "string", dest = "mapping_count",
                  help = "mapping count")
parser.add_option("-o", "--output1", action = "store", type = "string", dest = "output1",
                  help = "output: gene expression level file")
#parser.add_option("-n", "--output2", action = "store", type = "string", dest = "output2",
#                  help = "output: reads number falling into gene")



(options, args) = parser.parse_args()

if (options.uniq_map is None or
    options.mappable_exon_index is None or
    options.repeat_region is None or
    options.exon_combination is None or
    options.mapping_count is None or
    options.output1 is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_mapping_result = options.uniq_map
fname_exon = options.mappable_exon_index
fname_repeat = options.repeat_region
fname_length = options.exon_combination
fname_mapping_log = options.mapping_count
output_filename = options.output1
#output_filename2=options.output2

#fname_mapping_result='c:/Python26/test file/hg19_nonrepeat/Mapping_test.txt'
#fname_exon='C:/python26/test file/hg19_nonrepeat/ExonIndexNonRepeatRegion_chr1.txt'
#fname_repeat='C:/Python26/test file/hg19_nonrepeat/RepeatRegionChrs.txt'
#fname_length='C:/Python26/test file/hg19_nonrepeat/GencodeExonCombination_chr1.txt'
#fname_mapping_log='c:/Python26/test file/hg19_nonrepeat/75A3gencodev3c.log'
#output_filename='GeoncodeExonExpressionLevel_chr1.txt'

#fname_mapping_result=sys.argv[1]
#fname_exon=sys.argv[2]
#fname_repeat=sys.argv[3]
#fname_length=sys.argv[4]
#fname_mapping_log=sys.argv[5]
#output_filename=sys.argv[6]

try:
    bk_exon=open(fname_exon)
except:
    print prog_base + ": error: cannot open " + fname_exon
    sys.exit(1)
sh_exon=bk_exon.readlines()
row_number_exon=len(sh_exon)
bk_exon.close()

gene_express={}
exon_information={}
isoform_gene_pair={}
gene_strand={}
gene_read={}
gene_readnumber={}
for w in range(1,len(sh_exon)):###build dictory for exons for isoform
    temp1=sh_exon[w][0:-1].split('\t')
    gene_name=temp1[1]
    strand=temp1[2]
    gene_strand[gene_name]=strand
    gene_read[gene_name]=0
    gene_express[gene_name]=0
    gene_readnumber[gene_name]=0
    if temp1[4]!='[]':
        temp4=str(temp1[4])[2:-2].split('], [')
        temp5=[]
        exon_info=[]
        for vv in range(0,len(temp4)):
            temp6=temp4[vv].split(',')
            temp5=temp5+temp6
        exon_info=exon_info+temp5
        exon_information[temp1[0]]=exon_info
    isoform_gene_pair[temp1[0]]=gene_name

try:
    bk_repeat=open(fname_repeat)
except:
    print prog_base + ": error: cannot open " + fname_repeat
    sys.exit(1)
sh_repeat=bk_repeat.readlines()
row_number_repeat=len(sh_repeat)
bk_repeat.close()

repeat_info={}

for i in range(1,len(sh_repeat)):
    temp1=sh_repeat[i][0:-1].split('\t')
    gene_name=temp1[0]
    start=[]
    end=[]
    for j in range(1,len(temp1)):
        temp2=temp1[j][1:-1].split(',')
        start.append(int(temp2[0]))
        end.append(int(temp2[1]))
    RepeatLength=end[0]-start[0]+1
    for k in range(1,len(start)):
        if start[k]<end[k-1] and end[k]>=end[k-1]:
             RepeatLength=RepeatLength+end[k]-end[k-1]+1
        elif start[k]>end[k-1]:
             RepeatLength=RepeatLength+end[k]-start[k]+1

    repeat_info[gene_name]=RepeatLength
    
gene_length={}

try:
    bk_length=open(fname_length)
except:
    print prog_base + ": error: cannot open " + fname_length
    sys.exit(1)
sh_length=bk_length.readlines()
bk_length.close()

for i in range(1,len(sh_length)):
    temp1=sh_length[i][0:-1].split('\t')
    gene_old_length=int(temp1[2])
    genename=temp1[0]
    lengthchange=repeat_info.get(genename,'None')
    if lengthchange=='None':
        gene_length[genename]=gene_old_length
    else:
        gene_length[genename]=gene_old_length-lengthchange

try:
    bk_mapping_result=open(fname_mapping_result)
except:
    print prog_base + ": error: cannot open " + fname_mapping_result
    sys.exit(1)
sh_mapping_result=bk_mapping_result.readlines()
row_number_mapping_result=len(sh_mapping_result)
bk_mapping_result.close()
if row_number_mapping_result==0:
    print "No reads mapped to only one genen for current chromosome, so no expression level can be caculated base on repeat region library."
    sys.exit(0)

try:
    bk_mapping_log=open(fname_mapping_log)
except:
    print prog_base + ": error: cannot open " + fname_mapping_log
    sys.exit(1)
sh_mapping_log=bk_mapping_log.readlines()
row_number_mapping_log=len(sh_mapping_log)
bk_mapping_log.close()

###look for read_length and #mappedreads

temp_position=sh_mapping_log[1].find(':')
temp3=sh_mapping_log[0].find(':')
if temp3!=-1 and temp_position!=-1:
   mapped_reads=int(sh_mapping_log[0][temp3+1:-1])
   read_length=int(sh_mapping_log[1][temp_position+1:-1])


##############################################################################
####### the mapping result
chromosome = ''
for v in range(0,row_number_mapping_result): #
    mapping_info=sh_mapping_result[v][0:-1].split('\t')
    temp2=mapping_info[2].split('_',2)
    temp3=temp2[2].split(':',1)
    temp4=temp3[0].split('=',1)
    chromosome=temp4[1]
    mapping_info[2]=temp4[0]
    gene_name=isoform_gene_pair.get(mapping_info[2],'none')
    if gene_name=='none':
        continue
    gene_read[gene_name]=1
    if exon_information.get(mapping_info[2],'none')=='none':
        continue
    exon_info=exon_information[mapping_info[2]]
    for p in range(0,(len(exon_info)/4)):
        if p>=0: 
            if int(mapping_info[3])>=int(exon_info[4*p]) and int(mapping_info[3])+read_length-1<=(int(exon_info[4*p])+int(exon_info[4*p+3])-int(exon_info[4*p+2])):
                temp11=gene_express[gene_name]
                gene_express[gene_name]=temp11+read_length
            elif int(mapping_info[3])>=int(exon_info[4*p]) and int(mapping_info[3])<=(int(exon_info[4*p])+int(exon_info[4*p+3])-int(exon_info[4*p+2])) and int(mapping_info[3])+read_length-1>=(int(exon_info[4*p])+int(exon_info[4*p+3])-int(exon_info[4*p+2])):
                temp11=gene_express[gene_name]
                gene_express[gene_name]=temp11+(int(exon_info[4*p])+int(exon_info[4*p+3])-int(exon_info[4*p+2]))-int(mapping_info[3])+1
            elif int(mapping_info[3])<=int(exon_info[4*p]) and int(mapping_info[3])+read_length-1>=int(exon_info[4*p]) and int(mapping_info[3])+read_length-1<=(int(exon_info[4*p])+int(exon_info[4*p+3])-int(exon_info[4*p+2])):
                temp11=gene_express[gene_name]
                gene_express[gene_name]=temp11+int(mapping_info[3])+read_length-int(exon_info[4*p])
            elif int(mapping_info[3])<=int(exon_info[4*p]) and int(mapping_info[3])+read_length-1>=(int(exon_info[4*p])+int(exon_info[4*p+3])-int(exon_info[4*p+2])):
                temp11=gene_express[gene_name]
                gene_express[gene_name]=temp11+int(exon_info[4*p+3])-int(exon_info[4*p+2])+1
       
        
for v in gene_express.keys():
    if gene_express[v]!=0 and gene_length[v]!=0:
        rpkm=(float(gene_express[v])/read_length/gene_length[v]/int(mapped_reads))*1000000000
        generead=(float(gene_express[v]))/read_length
        gene_express[v]=rpkm
        gene_readnumber[v]=generead
    elif gene_express[v]==0 and gene_read[v]==1:
        gene_express[v]=0
    elif gene_express[v]==0:
        gene_readnumber[v]=0
 
try:
    file_output=open(output_filename, 'w')
except:
    print prog_base + ": error: cannot open " + output_filename
    sys.exit(1)


#note='None(RR) means there are some reads mapped to this gene, but after removing repetitive region, there is no reads left'+'\n'
method_manual='# There are 4 methods to deal with multi-mapped reads:'+'\n'
method_manual+='# 1.Unique: remove the multi-mapped reads directly.'+'\n'
method_manual+='# 2.Random: assign the multi-mapped reads randomly.'+'\n'
method_manual+='# 3.Proportion: assign the multi-mapped reads according to the proportion of gene expression level.'+'\n'
method_manual+='# 4.Mappable Library: only count the reads falling into mappable regions, and it is based on the unmappable repetitive region library, which can be built automatically.'+'\n'
method_manual+='# You chose the 4th one. For this method, only counts the reads falling into mappable regions, and it is based on the unmappable repetitive region library, which can be built automatically.'+'\n'
#file_output.writelines(method_manual)
titles='#GeneID'+'\t'+'Chromosome'+'\t'+'Strand'+'\t'+'Reads_Number'+'\t'+'RPKM'+'\n'
#file_output.writelines(note)
file_output.writelines(titles)
for v in gene_express.keys():
    insert_line=v+'\t'+chromosome+'\t'+gene_strand[v]+'\t'+str(gene_readnumber[v])+'\t'+str(gene_express[v])+'\n'
    file_output.writelines(insert_line)
file_output.close()
print 'GeneExpressionLevel Done for Chromosome %s' % (chromosome)


#try:
#    file_output2=open(output_filename2, 'w')
#except:
#    print prog_base + ": error: cannot open " + output_filename2
#    sys.exit(1)


#for v in gene_readnumber.keys():
#    insert_line=v+'\t'+chromosome+'\t'+gene_strand[v]+'\t'+str(gene_readnumber[v])+'\n'
#    file_output2.writelines(insert_line)
#file_output2.close()
#format(gene_read_number[v],'.2f')               
