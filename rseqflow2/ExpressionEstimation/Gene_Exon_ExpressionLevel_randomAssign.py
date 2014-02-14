#!/usr/bin/env python

'''
Created on 2012-08-15

@author: linliu
'''
import os
import sys
import optparse
import random
import re

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-m", "--reads-multiIsoformList", action = "store", type = "string", dest = "reads_multiIsoform_list",
		  help = "Transcript list of multiply-mapped reads")
parser.add_option("-i", "--exon-index", action = "store", type = "string", dest = "exon_index_file",
		  help = "txt file of exon index")		  
parser.add_option("-c", "--exon-combination", action = "store", type = "string", dest = "exon_combination_file",
		  help = "txt file of exon combination")		  
parser.add_option("-l", "--mapping-log", action = "store", type = "string", dest = "mapping_log_file",
		  help = "log file which containts information about mapped reads and read length")		  
parser.add_option("-g", "--gene-expression", action = "store", type = "string", dest = "gene_expression_level",
		  help = "output file of gene expression level")
parser.add_option("-e", "--exon-expression", action = "store", type = "string", dest = "exon_expression_level",
		  help = "output file of exon expression level")		    	  
		  

(options, args) = parser.parse_args()

if (options.reads_multiIsoform_list is None or
    options.exon_index_file is None or
    options.exon_combination_file is None or
    options.mapping_log_file is None or 
    options.gene_expression_level is None or
    options.exon_expression_level is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_read_isoformlist=options.reads_multiIsoform_list
fname_exon_index=options.exon_index_file
fname_exon_combination=options.exon_combination_file
fname_mapping_log=options.mapping_log_file
fname_output=options.gene_expression_level
fname_output2=options.exon_expression_level

bk_read_isoformlist=open(fname_read_isoformlist)
sh_read_isoformlist=bk_read_isoformlist.readlines()
row_number_read_isoformlist=len(sh_read_isoformlist)
bk_read_isoformlist.close()

bk_exon_index=open(fname_exon_index)
sh_exon_index=bk_exon_index.readlines()
row_number_exon_index=len(sh_exon_index)
bk_exon_index.close()

bk_exon_combination=open(fname_exon_combination)
sh_exon_combination=bk_exon_combination.readlines()
row_number_exon_combination=len(sh_exon_combination)
bk_exon_combination.close()

bk_mapping_log=open(fname_mapping_log)
sh_mapping_log=bk_mapping_log.readlines()
row_number_mapping_log=len(sh_mapping_log)
bk_mapping_log.close()
############get isoform_gene pair##############
isoform_gene_pair={}
isoform_exonInformation={}
for v in range(1, row_number_exon_index):
    temp=sh_exon_index[v][0:-1].split('\t')
#JSH
    # check if the line is a comment
#    if ('#' in temp[0]):
    if re.match('^#',temp[0]):
        continue

    transcript_id=temp[0]
    gene_id=temp[1]
    isoform_gene_pair[transcript_id]=gene_id
    isoform_exonInformation[transcript_id]=[]
    exon_list=temp[4][2:-2].split('), (')
    for vv in range(0, len(exon_list)):
        temp1=exon_list[vv].split(', ')
        isoform_exonInformation[transcript_id].append(temp1)
############get length of each gene##############   
gene_length={}
gene_strand={}
for v in range(1,row_number_exon_combination):
    temp1=sh_exon_combination[v].split('\t')
    gene_length[temp1[0]]=int(temp1[2])
    gene_strand[temp1[0]]=temp1[1]
############get number of mapped reads##############
for v in range(0,row_number_mapping_log):###look for read_length and #mappedreads
    temp1=sh_mapping_log[v].find('Read Length:')
    if temp1!=-1:
        read_length=int(sh_mapping_log[v][temp1+12:-1])
    temp2=sh_mapping_log[v].find('Number of mapped reads:')
    if temp2!=-1:
        mapped_reads=int(sh_mapping_log[v][temp2+23:-1])
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
############Assign reads randomly to a gene##############
gene_read_number={}
exon_express={}  
chromosome=''
for v in range(0, row_number_read_isoformlist):
    temp=sh_read_isoformlist[v][0:-1].split('\t')
    readID=temp[0]
    chromosome=temp[1]
    isoformlist=temp[2][0:-1].split(',')
    length=len(isoformlist)
    randomly=random.randint(0,length-1)
    temp2=isoformlist[randomly].split(':') 
    del randomly
    isoform_name=temp2[0]
    map_start=int(temp2[1])
    map_end=map_start+read_length-1
    #---------------gene expression level---------------#
    try:
        gene_name=isoform_gene_pair[isoform_name]
    except:
        print "Error:Transcript:%s is not found!" % (isoform_name)
        sys.exit(1)
    try:
        gene_read_number[gene_name]+=1
    except:
        gene_read_number[gene_name]=1
    #---------------exon expression level---------------#    
    exon_list=isoform_exonInformation[isoform_name]
    for vv in range(0, len(exon_list)):
        exon_location=exon_list[vv]
        exon_start=int(exon_location[0])
        exon_end=exon_start+int(exon_location[3])-int(exon_location[2])
        exon_name=gene_name+':'+str(exon_location[2])+':'+str(exon_location[3])
        if map_start>=exon_start and map_end<exon_end:
            try:
                exon_express[exon_name]+=read_length
            except:   
                exon_express[exon_name]=read_length
        elif map_start>=exon_start and map_start<exon_end and map_end>=exon_end:
            try:
                exon_express[exon_name]+=exon_end-map_stat
            except:
                exon_express[exon_name]=exon_end-map_start    
        elif map_start<=exon_start and map_end>=exon_start and map_end<exon_end:
            try:
                exon_express[exon_name]+=map_end-exon_start+1
            except:
                exon_express[exon_name]=map_end-exon_start+1    
        elif map_start<=exon_start and map_end>=exon_end:
            try:
                exon_express[exon_name]+=exon_end-exon_start
            except:
                exon_express[exon_name]=exon_end-exon_start   
gene_rpkm={}
for v in gene_read_number:
    rpkm=(float(gene_read_number[v])/gene_length[v]/mapped_reads)*1000000000
    gene_rpkm[v]=rpkm
        
exon_rpkm={}   
exon_read_number={}
for v in exon_express:
    temp=v.split(':')
    exon_length=int(temp[-1])-int(temp[-2])+1
    exon_read_number[v]=float(exon_express[v])/float(read_length)
    rpkm=(exon_read_number[v]/exon_length/mapped_reads)*1000000000
    exon_rpkm[v]=rpkm

file_output=open(fname_output, 'w')
method_manual='# There are 4 methods to deal with multi-mapped reads:'+'\n'
method_manual+='# 1.Unique: remove the multi-mapped reads directly.'+'\n'
method_manual+='# 2.Random: assign the multi-mapped reads randomly.'+'\n'
method_manual+='# 3.Proportion: assign the multi-mapped reads according to the proportion of gene expression level.'+'\n'
method_manual+='# 4.Mappable Library: only count the reads falling into mappable regions, and it is based on the unmappable repetitive region library, which can be built automatically.'+'\n'
method_manual+='# You chose the 2nd one. This method assigns the multi-mapped reads randomly.'+'\n'
#file_output.writelines(method_manual)
titles='#GeneID'+'\t'+'Chromosome'+'\t'+'Strand'+'\t'+'Reads_Number'+'\t'+'RPKM'+'\n'
file_output.writelines(titles)
for v in gene_read_number:
    insert_line=str(v)+'\t'+chromosome+'\t'+gene_strand[v]+'\t'+str(gene_read_number[v])+'\t'+str(gene_rpkm[v])+'\n'
    file_output.writelines(insert_line)
file_output.close()

file_output2=open(fname_output2, 'w')
method_manual='# There are 4 methods to deal with multi-mapped reads:'+'\n'
method_manual+='# 1.Unique: remove the multi-mapped reads directly.'+'\n'
method_manual+='# 2.Random: assign the multi-mapped reads randomly.'+'\n'
method_manual+='# 3.Proportion: assign the multi-mapped reads according to the proportion of gene expression level.'+'\n'
method_manual+='# 4.Mappable Library: only count the reads falling into mappable regions, and it is based on the unmappable repetitive region library, which can be built automatically.'+'\n'
method_manual+='# You chose the 2nd one. This method assigns the multi-mapped reads randomly.'+'\n'
#file_output2.writelines(method_manual)
titles='#GeneID:ExonStartLocation:ExonEndLocation'+'\t'+'Chromosome'+'\t'+'Strand'+'\t'+'Reads_Number'+'\t'+'RPKM'+'\n'
file_output2.writelines(titles)
for v in exon_express:
    temp=v.split(':')
    geneName=":".join(temp[0:-2])
    insert_line=str(v)+'\t'+chromosome+'\t'+gene_strand[geneName]+'\t'+format(exon_read_number[v],'.2f')+'\t'+str(exon_rpkm[v])+'\n'
    file_output2.writelines(insert_line)
file_output2.close( )
print 'GeneExpressionLevel (Assign Reads Randomly) Done for Chromosome %s' % (chromosome)
print 'ExonExpressionLevel (Assign Reads Randomly) Done for Chromosome %s' % (chromosome)
