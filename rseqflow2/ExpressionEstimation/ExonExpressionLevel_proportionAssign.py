#!/usr/bin/env python
'''
Created on 2012-08-16

@author: linliu
'''
import os
import sys
import optparse
import re

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-m", "--reads-multiIsoformList", action = "store", type = "string", dest = "reads_multiIsoform_list",
		  help = "Transcript list of multiply-mapped reads")
parser.add_option("-i", "--exon-index", action = "store", type = "string", dest = "exon_index_file",
		  help = "txt file of exon index")
#parser.add_option("-c", "--exon-combination", action = "store", type = "string", dest = "exon_combination_file",
#		  help = "txt file of exon combination")		  
parser.add_option("-u", "--unique-exonExpression", action = "store", type = "string", dest = "unique_exon_expressionLevel",
		  help = "exon expression level from uniquely-gene mapped reads data")
parser.add_option("-l", "--mapping-log", action = "store", type = "string", dest = "mapping_log_file",
		  help = "log file which containts information about mapped reads and read length")
parser.add_option("-o", "--exon-expression", action = "store", type = "string", dest = "exon_expression_level",
		  help = "output file of exon expression level")		    	  
		  

(options, args) = parser.parse_args()

if (options.reads_multiIsoform_list is None or
    options.exon_index_file is None or
    options.unique_exon_expressionLevel is None or
    options.mapping_log_file is None or
    options.exon_expression_level is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_read_isoformlist=options.reads_multiIsoform_list
fname_exon_index=options.exon_index_file
#fname_exon_combination=options.exon_combination_file
fname_uniq_exonExpression=options.unique_exon_expressionLevel
fname_mapping_log=options.mapping_log_file
fname_output=options.exon_expression_level

bk_read_isoformlist=open(fname_read_isoformlist)
sh_read_isoformlist=bk_read_isoformlist.readlines()
row_number_read_isoformlist=len(sh_read_isoformlist)
bk_read_isoformlist.close()

bk_exon_index=open(fname_exon_index)
sh_exon_index=bk_exon_index.readlines()
row_number_exon_index=len(sh_exon_index)
bk_exon_index.close()

#bk_exon_combination=open(fname_exon_combination)
#sh_exon_combination=bk_exon_combination.readlines()
#row_number_exon_combination=len(sh_exon_combination)
#bk_exon_combination.close()

bk_uniq_exonExpression=open(fname_uniq_exonExpression)
sh_uniq_exonExpression=bk_uniq_exonExpression.readlines()
row_number_uniq_exonExpression=len(sh_uniq_exonExpression)
bk_uniq_exonExpression.close()

bk_mapping_log=open(fname_mapping_log)
sh_mapping_log=bk_mapping_log.readlines()
row_number_mapping_log=len(sh_mapping_log)
bk_mapping_log.close()
############get isoform_gene pair and exon index##############
isoform_gene_pair={}
isoform_exonInformation={}
gene_strand={}
chromosome=''
for v in range(1, row_number_exon_index):
    temp=sh_exon_index[v][0:-1].split('\t')
#JSH
# check if the line is a comment
#    if ('#' in temp[0]):
    if re.match('^#',temp[0]):
        continue

    transcript_id=temp[0]
    gene_id=temp[1]
    gene_strand[gene_id]=temp[2]
    isoform_gene_pair[transcript_id]=gene_id
    isoform_exonInformation[transcript_id]=[]
    exon_list=temp[4][2:-2].split('), (')
    for vv in range(0, len(exon_list)):
        temp1=exon_list[vv].split(', ')
        isoform_exonInformation[transcript_id].append(temp1)
############get length of each gene##############   
#gene_length={}
#gene_strand={}
#for v in range(1,row_number_exon_combination):
#    temp1=sh_exon_combination[v].split('\t')
#    gene_length[temp1[0]]=int(temp1[2])
#    gene_strand[temp1[0]]=temp1[1]    
############get rpkm of each gene##############   
exon_uniq_rpkm={}
for v in range(1,row_number_uniq_exonExpression):
    temp1=sh_uniq_exonExpression[v][0:-1].split('\t')
#JSH
# check if the line is a comment
    if re.match('^#',temp1[0]):
        continue
    exon_uniq_rpkm[temp1[0]]=float(temp1[3])
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
    sys.exit(0)
try:
    mapped_reads+=0
except:   
    print "Error:Can't read number of mapped reads from " + str(fname_mapping_log)
    sys.exit(0)          
############Assign reads base on proportion of exons##############
exon_express={} 
for v in range(0, row_number_read_isoformlist):
    temp=sh_read_isoformlist[v][0:-1].split('\t')
   
    readID=temp[0]
    chromosome=temp[1]
    isoformlist=temp[2][0:-1].split(',')
    length=len(isoformlist)
    exon_proportion={}
    exon_covered_bp={}
    total=0.0
    for vv in range(0, length):
        temp2=isoformlist[vv].split(':') 
        isoform_name=temp2[0]
        map_start=int(temp2[1])
        map_end=map_start+read_length-1
        exon_list=isoform_exonInformation[isoform_name]
        try:
            gene_name=isoform_gene_pair[isoform_name]
        except:
            print "Error:Transcript:%s is not found!" % (isoform_name)
            sys.exit(1)
        for w in range(0, len(exon_list)):
            exon_location=exon_list[w]
            exon_start=int(exon_location[0])
            exon_end=exon_start+int(exon_location[3])-int(exon_location[2])
            exon_name=gene_name+':'+str(exon_location[2])+':'+str(exon_location[3])
            found=False
            if map_start>=exon_start and map_end<exon_end:
                found=True
                exon_covered_bp[exon_name]=read_length   
            elif map_start>=exon_start and map_start<exon_end and map_end>=exon_end:
                found=True
                exon_covered_bp[exon_name]=exon_end-map_start    
            elif map_start<=exon_start and map_end>=exon_start and map_end<exon_end:
                found=True
                exon_covered_bp[exon_name]=map_end-exon_start+1    
            elif map_start<=exon_start and map_end>=exon_end:
                found=True
                exon_covered_bp[exon_name]=exon_end-exon_start
            if found==True:
                try:
                    temp_rpkm=exon_uniq_rpkm[exon_name]
                except:
                    temp_rpkm=0
                exon_proportion[exon_name]=temp_rpkm    
                total+=temp_rpkm
    if total==0:
        total+=len(exon_proportion.keys())
        for e in exon_proportion:
            proportion=exon_covered_bp[e]*1/total
            try:
                exon_express[e]+=proportion
            except:
                exon_express[e]=proportion
    else:        
        for e in exon_proportion:
            proportion=exon_covered_bp[e]*exon_proportion[e]/total
            try:
                exon_express[e]+=proportion
            except:
                exon_express[e]=proportion

exon_multi_rpkm={}
exon_read_number={}
for v in exon_express:
    temp=v.split(':')
    exon_length=int(temp[-1])-int(temp[-2])+1
    exon_read_number[v]=float(exon_express[v])/float(read_length)
    rpkm=(exon_read_number[v]/exon_length/mapped_reads)*1000000000
    exon_multi_rpkm[v]=rpkm

file_output=open(fname_output, 'w')
method_manual='# There are 4 methods to deal with multi-mapped reads:'+'\n'
method_manual+='# 1.Unique: remove the multi-mapped reads directly.'+'\n'
method_manual+='# 2.Random: assign the multi-mapped reads randomly.'+'\n'
method_manual+='# 3.Proportion: assign the multi-mapped reads according to the proportion of gene expression level.'+'\n'
method_manual+='# 4.Mappable Library: only count the reads falling into mappable regions, and it is based on the unmappable repetitive region library, which can be built automatically.'+'\n'
method_manual+='# You chose the 3rd one. This method assigns the multi-mapped reads according to the proportion of gene expression level.'+'\n'
#file_output.writelines(method_manual)
titles='#GeneID:ExonStartLocation:ExonEndLocation'+'\t'+'Chromosome'+'\t'+'Strand'+'\t'+'Reads_Number'+'\t'+'RPKM'+'\n'
file_output.writelines(titles)
for v in exon_express:
    temp=v.split(':')
    geneName=":".join(temp[0:-2])
    insert_line=str(v)+'\t'+chromosome+'\t'+gene_strand[geneName]+'\t'+format(exon_read_number[v],'.2f')+'\t'+str(exon_multi_rpkm[v])+'\n'
    file_output.writelines(insert_line)
file_output.close( )
print 'ExonExpressionLevel (Assign Reads Based on Proportion) Done for Chromosome %s' % (chromosome)
