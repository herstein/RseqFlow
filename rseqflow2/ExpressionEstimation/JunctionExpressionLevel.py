#!/usr/bin/env python

import os
import sys
import optparse

#fname_mapping_result='/home/cmb-03/tc1/yw_003/ExtentionRseqFlow/ExpressionEstimation/UniqueMap//test_input/test_uniq.sam'
#fname_junction='/home/cmb-03/tc1/yw_003/ExtentionRseqFlow/ExpressionEstimation/UniqueMap//test_input/JunctionIndex_human_gencodeV12_chr3.txt'
#fname_mapping_log='/home/cmb-03/tc1/yw_003/ExtentionRseqFlow/ExpressionEstimation/UniqueMap//test_input/uniq_mapping_info.log'
#output_filename='/home/cmb-03/tc1/yw_003/ExtentionRseqFlow/ExpressionEstimation/UniqueMap//test_output/geoncodeV12JunctionExpressionLevel_chr3.txt'

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-u", "--uniqueMapped-result", action = "store", type = "string", dest = "uniqueMapped_sam_file",
                  help = "Transcript list of uniquely-mapped reads")
parser.add_option("-i", "--junction-index", action = "store", type = "string", dest = "junction_index_file",
                  help = "txt file of junction index")
parser.add_option("-l", "--mapping-log", action = "store", type = "string", dest = "mapping_log_file",
                  help = "log file which containts information about mapped reads and read length")
parser.add_option("-o", "--junction-expression", action = "store", type = "string", dest = "junction_expression_level",
                  help = "output file of junction expression level")


(options, args) = parser.parse_args()

if (options.uniqueMapped_sam_file is None or
    options.junction_index_file is None or
    options.mapping_log_file is None or
    options.junction_expression_level is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_mapping_result=options.uniqueMapped_sam_file
fname_junction=options.junction_index_file
fname_mapping_log=options.mapping_log_file
output_filename=options.junction_expression_level

bk_junc=open(fname_junction)
sh_junc=bk_junc.readlines()
row_number_junc=len(sh_junc)
bk_junc.close()

#JSH 2013-08-26
### create empty file if chr not present (required for PEGASUS workflow) ###
if(row_number_junc==0):
    print "JunctionExpressionLevel.py: No chr data for " + fname_junction
    try:
        open(output_filename,'w').close()
    except:
        print prog_base + ": error: cannot create file " + output_filename
        sys.exit(1)

    sys.exit(0)


junc_express={}
junc_information={}
isoform_gene_pair={}
gene_strand={}
for w in range(1,len(sh_junc)):###build dictory for junction for isoform
    temp1=sh_junc[w][0:-1].split('\t')
    gene_name=temp1[1]
    strand=temp1[2]
    gene_strand[gene_name]=strand
    temp4=str(temp1[4])[1:-2].split('), ')
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

bk_mapping_log=open(fname_mapping_log)
sh_mapping_log=bk_mapping_log.readlines()
row_number_mapping_log=len(sh_mapping_log)
bk_mapping_log.close()
for w in range(0,len(sh_mapping_log)):###look for read_length and #mappedreads
    temp1=sh_mapping_log[w].find('Read Length:')
    if temp1!=-1:
        read_length=int(sh_mapping_log[w][temp1+12:-1])
    temp2=sh_mapping_log[w].find('Number of mapped reads:')
    if temp2!=-1:
        mapped_reads=int(sh_mapping_log[w][temp2+23:-1])
#############check####################
try:
    read_length+=0
except:    
    print "Error:Can't read read length from " + str(fname_mapping_log)
    sys.exit(1)
try:
    read_length+=0
except:    
    print "Error:Can't read number of mapped reads from " + str(fname_mapping_log)
    sys.exit(1)
  
bk_mapping_result=open(fname_mapping_result)
#sh_mapping_result=bk_mapping_result.readlines()
#row_number_mapping_result=len(sh_mapping_result)
#bk_mapping_result.close()


isoform_name=[]
chromosome = ''
 
##############################################################################
####### the mapping result
count=0
for v in bk_mapping_result.xreadlines(): #
    temp=v.split('\t')
    temp1=temp[2].split('_',2)
    temp2=temp1[2].split('=')
    temp4=temp2[1].split(':')
    chromosome=temp4[0]
    mapping_info=[]
    mapping_info.append(temp2[0])                #####take out the isoform name
    mapping_info.append(temp[3])
    isoform_name.append(mapping_info[0])
    try:
        gene_name=isoform_gene_pair[mapping_info[0]]
    except:
        print "Warning:Transcript:%s is not found!" % (mapping_info[0])
        continue
    
    junc_info=[]
    junc_info=junc_information[mapping_info[0]]
    if junc_info!=[]:   #if there are junctions in the isoform
        for p in range(0,(len(junc_info)/4)):
            if int(mapping_info[1])>=(int(junc_info[4*p])-read_length+2) and int(mapping_info[1])<=(int(junc_info[4*p])):
                temp22=gene_name+':'+str(int(junc_info[4*p+2]))+':'+str(int(junc_info[4*p+3]))
                temp11=junc_express[temp22]
                junc_express[temp22]=temp11+1
bk_mapping_result.close()
junc_read_number={}
for v in junc_express.keys():
    junc_read_number[v]=junc_express[v]
    read_ratio=(float(junc_read_number[v])/mapped_reads)*1000000
    junc_express[v]=read_ratio


junc_express_list=junc_express.items()                        
file_output=open(output_filename, 'w')
method_manual='# There are 4 methods to deal with multi-mapped reads:'+'\n'
method_manual+='# 1.Unique: remove the multi-mapped reads directly.'+'\n'
method_manual+='# 2.Random: assign the multi-mapped reads randomly.'+'\n'
method_manual+='# 3.Proportion: assign the multi-mapped reads according to the proportion of gene expression level.'+'\n'
method_manual+='# 4.Mappable Library: only count the reads falling into mappable regions, and it is based on the unmappable repetitive region library, which can be built automatically.'+'\n'
method_manual+='# You chose the 1st one. This method removes the multi-mapped reads directly.'+'\n'
#file_output.writelines(method_manual)
titles='#GeneID:JunctionStartLocation:JunctionEndLocation'+'\t'+'Chromosome'+'\t'+'Strand'+'\t'+'Reads_Number'+'\t'+'#Reads/#Mapped Reads(Million)'+'\n'
file_output.writelines(titles)
for v in junc_express:
    if junc_express[v]!=0:
        temp=v.split(':')
	geneName=":".join(temp[0:-2])
        insert_line=str(v)+'\t'+chromosome+'\t'+gene_strand[geneName]+'\t'+str(junc_read_number[v])+'\t'+str(junc_express[v])+'\n'
        file_output.writelines(insert_line)
file_output.close( )
print 'JunctionExpressionLevel (Uniquely Mapped Reads) Done for Chromosome %s' % (chromosome)
