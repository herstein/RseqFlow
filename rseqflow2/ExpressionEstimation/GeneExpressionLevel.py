#!/usr/bin/python
'''
Created on 2010-2-7

@author: yingwang
'''

import os
import sys
import optparse
import re
#fname_mapped='/home/cmb-03/tc1/yw_003/ExtentionRseqFlow/ExpressionEstimation/UniqueMap/test_input/test_uniq.sam'
#fname_gene_isoform_pair='/home/cmb-03/tc1/yw_003/ExtentionRseqFlow/ExpressionEstimation/UniqueMap/test_input/ExonIndex_human_gencodeV12_chr3.txt'
#fname_gene_length='/home/cmb-03/tc1/yw_003/ExtentionRseqFlow/ExpressionEstimation/UniqueMap/test_input/ExonCombination_human_gencodeV12_chr3.txt'
#fname_mapping_log='/home/cmb-03/tc1/yw_003/ExtentionRseqFlow/ExpressionEstimation/UniqueMap/test_input/uniq_mapping_info.log'
#output_filename='/home/cmb-03/tc1/yw_003/ExtentionRseqFlow/ExpressionEstimation/UniqueMap/test_output/geoncodeV12GeneExpressioniLevel_chr3.txt'

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-u", "--uniqueMapped-result", action = "store", type = "string", dest = "uniqueMapped_sam_file",
		  help = "Transcript list of multiply-mapped reads")
parser.add_option("-i", "--exon-index", action = "store", type = "string", dest = "exon_index_file",
		  help = "txt file of exon index")		  
parser.add_option("-c", "--exon-combination", action = "store", type = "string", dest = "exon_combination_file",
		  help = "txt file of exon combination")
parser.add_option("-l", "--mapping-log", action = "store", type = "string", dest = "mapping_log_file",
		  help = "log file which containts information about mapped reads and read length")		  
parser.add_option("-o", "--gene-expression", action = "store", type = "string", dest = "gene_expression_level",
		  help = "output file of gene expression level")	    	  
		  

(options, args) = parser.parse_args()

if (options.uniqueMapped_sam_file is None or
    options.exon_index_file is None or
    options.exon_combination_file is None or
    options.mapping_log_file is None or 
    options.gene_expression_level is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_mapped=options.uniqueMapped_sam_file
fname_gene_isoform_pair=options.exon_index_file
fname_gene_length=options.exon_combination_file
fname_mapping_log=options.mapping_log_file
output_filename=options.gene_expression_level


bk_mapped=open(fname_mapped)
#sh_mapped=bk_mapped.readlines()
#row_number_mapped=len(sh_mapped)
#bk_mapped.close()

bk_pair=open(fname_gene_isoform_pair)
sh_pair=bk_pair.readlines()
row_number_pair=len(sh_pair)
bk_pair.close()

#JSH 2013-08-26
### create empty file if chr not present (required for PEGASUS workflow) ###
if(row_number_pair==0):
    print "GeneExpressionLevel.py: No chr data for " + fname_gene_isoform_pair
    try:
        open(output_filename,'w').close()
    except:
        print prog_base + ": error: cannot create file " + output_filename
        sys.exit(1)

    sys.exit(0)


bk_length=open(fname_gene_length)
sh_length=bk_length.readlines()
row_number_length=len(sh_length)
bk_length.close()

#JSH 2013-08-26
### create empty file if chr not present (required for PEGASUS workflow) ###
if(row_number_length==0):
    print "GeneExpressionLevel.py: No chr data for " + fname_gene_length
    try:
        open(output_filename,'w').close()
    except:
        print prog_base + ": error: cannot create file " + output_filename
        sys.exit(1)

    sys.exit(0)


gene_length={}
gene_strand={}
for w in range(1,row_number_length):
    temp1=sh_length[w].split('\t')
    gene_length[temp1[0]]=int(temp1[2])
    gene_strand[temp1[0]]=temp1[1]

bk_mapping_log=open(fname_mapping_log)
sh_mapping_log=bk_mapping_log.readlines()
row_number_mapping_log=len(sh_mapping_log)
bk_mapping_log.close()
#read_length=-1
#mapped_reads=-1
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
    mapped_reads+=0
except:
    print "Error:Can't read number of mapped reads from " + str(fname_mapping_log)
    sys.exit(1)


isoform_gene_pair={}
gene_read_number={}
gene_rpkm={}
for v in range(1,row_number_pair):
    temp=sh_pair[v].split('\t')
    transcript_id=temp[0]
    gene_id=temp[1]
    isoform_gene_pair[transcript_id]=gene_id
    #gene_read_number[gene_id]=0

mapped_info=[]
chromosome = ''
for v in bk_mapped.xreadlines():
    temp=v.split('\t')
    temp1=temp[2].split('_',2)
    #JSH
    # Check for chr
    if re.search(r'=chr',temp1[2]):
        temp2=temp1[2].split('=')
        temp4=temp2[1].split(':')
        chromosome=temp4[0]
        mapped_info.append(temp2[0])
    else:
 	continue
bk_mapped.close()    
    
for v in mapped_info:
    try:
        gene_name=isoform_gene_pair[v]
    except:
        print "Warning:Transcript:%s is not found!" % (v)
        continue
    try:
        gene_read_number[gene_name]+=1
    except:
        gene_read_number[gene_name]=1        
        
for v in gene_read_number.keys():
    rpkm=(float(gene_read_number[v])/gene_length[v]/mapped_reads)*1000000000
    gene_rpkm[v]=rpkm

    #JSH
    #print "GeneName:%s  gene_read_number:%s gene_length:%s  mapped_reads:%s rpkm:%s" % (v,gene_read_number[v],gene_length[v],mapped_reads,rpkm)

    

file_output=open(output_filename, 'w')
method_manual='# There are 4 methods to deal with multi-mapped reads:'+'\n'
method_manual+='# 1.Unique: remove the multi-mapped reads directly.'+'\n'
method_manual+='# 2.Random: assign the multi-mapped reads randomly.'+'\n'
method_manual+='# 3.Proportion: assign the multi-mapped reads according to the proportion of gene expression level.'+'\n'
method_manual+='# 4.Mappable Library: only count the reads falling into mappable regions, and it is based on the unmappable repetitive region library, which can be built automatically.'+'\n'
method_manual+='# You chose the 1st one. This method removes the multi-mapped reads directly.'+'\n'
#file_output.writelines(method_manual)
titles='#GeneID'+'\t'+'chromosome'+'\t'+'Strand'+'\t'+'Reads_Number'+'\t'+'RPKM'+'\n'
file_output.writelines(titles)
for v in gene_read_number.keys():
    insert_line=str(v)+'\t'+chromosome+'\t'+gene_strand[v]+'\t'+str(gene_read_number[v])+'\t'+str(gene_rpkm[v])+'\n'
    file_output.writelines(insert_line)
file_output.close()
print 'GeneExpressionLevel (Uniquely Mapped Reads) Done for Chromosome %s' % (chromosome)
