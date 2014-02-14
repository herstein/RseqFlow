#!/usr/bin/env python

'''
Created on 2012-08-16

@author: linliu
'''
import os
import sys
import optparse

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-m", "--reads-multiIsoformList", action = "store", type = "string", dest = "reads_multiIsoform_list",
		  help = "Transcript list of multiply-mapped reads")
parser.add_option("-i", "--exon-index", action = "store", type = "string", dest = "exon_index_file",
		  help = "txt file of exon index")
parser.add_option("-c", "--exon-combination", action = "store", type = "string", dest = "exon_combination_file",
		  help = "txt file of exon combination")		  
parser.add_option("-u", "--unique-geneExpression", action = "store", type = "string", dest = "unique_gene_expressionLevel",
		  help = "gene expression level from uniquely-gene mapped reads data")  	  
parser.add_option("-l", "--mapping-log", action = "store", type = "string", dest = "mapping_log_file",
		  help = "log file which containts information about mapped reads and read length")
parser.add_option("-o", "--gene-expression", action = "store", type = "string", dest = "gene_expression_level",
		  help = "output file of gene expression level")
		  

(options, args) = parser.parse_args()

if (options.reads_multiIsoform_list is None or
    options.exon_index_file is None or
    options.exon_combination_file is None or
    options.unique_gene_expressionLevel is None or
    options.mapping_log_file is None or
    options.gene_expression_level is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_read_isoformlist=options.reads_multiIsoform_list
fname_exon_index=options.exon_index_file
fname_exon_combination=options.exon_combination_file
fname_uniq_geneExpression=options.unique_gene_expressionLevel
fname_mapping_log=options.mapping_log_file
fname_output=options.gene_expression_level

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

bk_uniq_geneExpression=open(fname_uniq_geneExpression)
sh_uniq_geneExpression=bk_uniq_geneExpression.readlines()
row_number_uniq_geneExpression=len(sh_uniq_geneExpression)
bk_uniq_geneExpression.close()

bk_mapping_log=open(fname_mapping_log)
sh_mapping_log=bk_mapping_log.readlines()
row_number_mapping_log=len(sh_mapping_log)
bk_mapping_log.close()
############get isoform_gene pair##############
isoform_gene_pair={}
isoform_exonInformation={}
for v in range(1, row_number_exon_index):
    temp=sh_exon_index[v][0:-1].split('\t')
    transcript_id=temp[0]
    gene_id=temp[1]
    isoform_gene_pair[transcript_id]=gene_id
############get length of each gene##############   
gene_length={}
gene_strand={}
for v in range(1,row_number_exon_combination):
    temp1=sh_exon_combination[v].split('\t')
    gene_length[temp1[0]]=int(temp1[2])
    gene_strand[temp1[0]]=temp1[1]    
############get rpkm of each gene##############   
gene_uniq_rpkm={}
for v in range(1,row_number_uniq_geneExpression):
    temp1=sh_uniq_geneExpression[v][0:-1].split('\t')
#JSH
# check if the line is a comment
    if (not '#' in temp1[0]):
       gene_uniq_rpkm[temp1[0]]=float(temp1[4])

#    temp1=sh_uniq_geneExpression[v][0:-1].split('\t')
#    gene_uniq_rpkm[temp1[0]]=float(temp1[4])
############get information of mapping##############
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
############Assign reads base on proportion of genes##############
gene_read_number={} 
chromosome=''
for v in range(0, row_number_read_isoformlist):
    temp=sh_read_isoformlist[v][0:-1].split('\t')
    readID=temp[0]
    chromosome=temp[1]
    isoformlist=temp[2][0:-1].split(',')
    length=len(isoformlist)
    gene_proportion={}
    total=0.0
    for vv in range(0, length):
        temp2=isoformlist[vv].split(':') 
        isoform_name=temp2[0]
        try:
            gene_name=isoform_gene_pair[isoform_name]
        except:
            print "Error:Transcript:%s is not found!" % (isoform_name)
            sys.exit(1)
        try:    
            temp_rpkm=gene_uniq_rpkm[gene_name]
        except:
            temp_rpkm=0
        gene_proportion[gene_name]=temp_rpkm    
        total+=temp_rpkm
    if total==0:
        total+=len(gene_proportion.keys()) #let total = float(length)
        for g in gene_proportion:
            proportion=1/total
            try:
                gene_read_number[g]+=proportion
            except:    
                gene_read_number[g]=proportion
    else:            
        for g in gene_proportion:
            proportion=gene_proportion[g]/total
            try:
                gene_read_number[g]+=proportion
            except:    
                gene_read_number[g]=proportion

gene_multi_rpkm={}
for v in gene_read_number:
    rpkm=(float(gene_read_number[v])/gene_length[v]/mapped_reads)*1000000000
    gene_multi_rpkm[v]=rpkm

file_output=open(fname_output, 'w')
method_manual='# There are 4 methods to deal with multi-mapped reads:'+'\n'
method_manual+='# 1.Unique: remove the multi-mapped reads directly.'+'\n'
method_manual+='# 2.Random: assign the multi-mapped reads randomly.'+'\n'
method_manual+='# 3.Proportion: assign the multi-mapped reads according to the proportion of gene expression level.'+'\n'
method_manual+='# 4.Mappable Library: only count the reads falling into mappable regions, and it is based on the unmappable repetitive region library, which can be built automatically.'+'\n'
method_manual+='# You chose the 3rd one. This method assigns the multi-mapped reads according to the proportion of gene expression level.'+'\n'
#file_output.writelines(method_manual)
titles='#GeneID'+'\t'+'Chromosome'+'\t'+'Strand'+'\t'+'Reads_Number'+'\t'+'RPKM'+'\n'
file_output.writelines(titles)
for v in gene_read_number:
    insert_line=str(v)+'\t'+chromosome+'\t'+gene_strand[v]+'\t'+format(gene_read_number[v],'.2f')+'\t'+str(gene_multi_rpkm[v])+'\n'
    file_output.writelines(insert_line)
file_output.close()
print 'GeneExpressionLevel (Assign Reads Based on Proportion) Done for Chromosome %s' % (chromosome)
