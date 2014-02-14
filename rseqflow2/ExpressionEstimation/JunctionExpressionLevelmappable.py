#!/usr/bin/env python

import os
import sys
import optparse

# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-s", "--uniq-map", action = "store", type = "string", dest = "uniq_map",
		  help = "uniq mapped sam file for each chromosome")
parser.add_option("-j", "--mappable-junction-index", action = "store", type = "string", dest = "mappable_junction_index",
		  help = "junction index no repeat regions")
parser.add_option("-m", "--mapping-count", action = "store", type = "string", dest = "mapping_count",
		  help = "mapping count")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output",
		  help = "output: junction expression level file")

(options, args) = parser.parse_args()

if (options.uniq_map is None or
    options.mappable_junction_index is None or
    options.mapping_count is None or
    options.output is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_mapping_result = options.uniq_map
fname_junction = options.mappable_junction_index
fname_mapping_log = options.mapping_count
output_filename = options.output



try:
    bk_junc=open(fname_junction)
except:
    print prog_base + ": error: cannot open file " + fname_junction
    sys.exit(1)
sh_junc=bk_junc.readlines()
row_number_junc=len(sh_junc)
bk_junc.close()

junc_express={}
junc_information={}
isoform_gene_pair={}
gene_strand={}

for w in range(1,len(sh_junc)):###build dictory for junction for isoform
    temp1=sh_junc[w][0:-1].split('\t')
    gene_name=temp1[1]
    strand=temp1[2]
    gene_strand[gene_name]=strand
    temp4=str(temp1[4])[1:-2].split('], ')
    temp5=[]
    junc_info=[]
    for vv in range(0,len(temp4)):
        temp4[vv]=temp4[vv][1:]
        temp6=temp4[vv].split(',')
        junc_express[str(gene_name)+':'+str(int(temp6[2]))+':'+str(int(temp6[3]))]=0
        temp5=temp5+temp6
    junc_info=junc_info+temp5
    junc_information[temp1[0]]=junc_info
    isoform_gene_pair[temp1[0]]=gene_name
    
try:
    bk_mapping_log=open(fname_mapping_log)
except:
    print prog_base + ": error: cannot open file " + fname_mapping_log
    sys.exit(1)
sh_mapping_log=bk_mapping_log.readlines()
row_number_mapping_log=len(sh_mapping_log)
bk_mapping_log.close()

temp_position=sh_mapping_log[1].find(':')
temp3=sh_mapping_log[0].find(':')
if temp3!=-1 and temp_position!=-1:
   mapped_reads=int(sh_mapping_log[0][temp3+1:-1])
   read_length=int(sh_mapping_log[1][temp_position+1:-1])

  
try:
    bk_mapping_result=open(fname_mapping_result)
except:
    print prog_base + ": error: cannot open file " + fname_mapping_result
    sys.exit(1)
sh_mapping_result=bk_mapping_result.readlines()
row_number_mapping_result=len(sh_mapping_result)
bk_mapping_result.close()
if row_number_mapping_result==0:
    print "No reads mapped to only one genen for current chromosome, so no expression level can be caculated base on repeat region library."
    sys.exit(0)

isoform_name=[]
chromosome = ''
 
##############################################################################
####### the mapping result
for v in range(0,row_number_mapping_result): #
    mapping_info=sh_mapping_result[v][0:-1].split('\t')
    temp2=mapping_info[2].split('_',2)
    temp3=temp2[2].split(':',1)
    temp4=temp3[0].split('=',1)
    chromosome=temp4[1]
    mapping_info[2]=temp4[0]                #####take out the isoform name
    gene_name=isoform_gene_pair.get(mapping_info[2],'none')
    isoform_name.append(mapping_info[2])
    if gene_name=='none':
        continue
    
    junc_info=[]
    junc_info=junc_information[isoform_name[v]]
    if junc_info!=[]:   #if there are junctions in the isoform
        for p in range(0,(len(junc_info)/4)):
            if int(mapping_info[3])>=(int(junc_info[4*p])-read_length+2) and int(mapping_info[3])<=(int(junc_info[4*p])):
                temp22=gene_name+':'+str(int(junc_info[4*p+2]))+':'+str(int(junc_info[4*p+3]))
                temp11=junc_express[temp22]
                junc_express[temp22]=temp11+1
junc_ratio={}
for v in junc_express.keys():
    read_ratio=(float(junc_express[v])/mapped_reads)*1000000
    junc_ratio[v]=read_ratio
junc_express_list=junc_express.items()                        

try:
    file_output=open(output_filename, 'w')
except:
    print prog_base + ": error: cannot open file " + output_filename
    sys.exit(1)
method_manual='# There are 4 methods to deal with multi-mapped reads:'+'\n'
method_manual+='# 1.Unique: remove the multi-mapped reads directly.'+'\n'
method_manual+='# 2.Random: assign the multi-mapped reads randomly.'+'\n'
method_manual+='# 3.Proportion: assign the multi-mapped reads according to the proportion of gene expression level.'+'\n'
method_manual+='# 4.Mappable Library: only count the reads falling into mappable regions, and it is based on the unmappable repetitive region library, which can be built automatically.'+'\n'
method_manual+='# You chose the 4th one. For this method, only counts the reads falling into mappable regions, and it is based on the unmappable repetitive region library, which can be built automatically.'+'\n'
#file_output.writelines(method_manual)
titles='#GeneID:JunctionStartLocation:JunctionEndLocation'+'\t'+'Chromosome'+'\t'+'Strand'+'\t'+'Reads_Number'+'\t'+'#Reads/#Mapped Reads(Million)'+'\n'
file_output.writelines(titles)

for l in range(0,len(junc_express_list)):
    temp=junc_express_list[l][0].split(':')
    geneName=":".join(temp[0:-2])
    junc_id=junc_express_list[l][0]
    insert_line=str(junc_express_list[l][0])+'\t'+chromosome+'\t'+gene_strand[geneName]+'\t'+str(junc_ratio[junc_id])+'\t'+str(junc_express_list[l][1])+'\n'
    file_output.writelines(insert_line)
file_output.close( )

print 'JunctionExpressionLevel Done for Chromosome %s' % (chromosome)     
