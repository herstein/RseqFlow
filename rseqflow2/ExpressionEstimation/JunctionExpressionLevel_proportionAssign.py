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
parser.add_option("-i", "--junc-index", action = "store", type = "string", dest = "junc_index_file",
		  help = "txt file of junc index")
#parser.add_option("-c", "--junc-combination", action = "store", type = "string", dest = "junc_combination_file",
#		  help = "txt file of junc combination")		  
parser.add_option("-u", "--unique-juncExpression", action = "store", type = "string", dest = "unique_junc_expressionLevel",
		  help = "junc expression level from uniquely-gene mapped reads data")
parser.add_option("-l", "--mapping-log", action = "store", type = "string", dest = "mapping_log_file",
		  help = "log file which containts information about mapped reads and read length")
parser.add_option("-o", "--junc-expression", action = "store", type = "string", dest = "junc_expression_level",
		  help = "output file of junc expression level")		    	  
		  

(options, args) = parser.parse_args()

if (options.reads_multiIsoform_list is None or
    options.junc_index_file is None or
    options.unique_junc_expressionLevel is None or
    options.mapping_log_file is None or
    options.junc_expression_level is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_read_isoformlist=options.reads_multiIsoform_list
fname_junc_index=options.junc_index_file
#fname_junc_combination=options.junc_combination_file
fname_uniq_juncExpression=options.unique_junc_expressionLevel
fname_mapping_log=options.mapping_log_file
fname_output=options.junc_expression_level

bk_read_isoformlist=open(fname_read_isoformlist)
sh_read_isoformlist=bk_read_isoformlist.readlines()
row_number_read_isoformlist=len(sh_read_isoformlist)
bk_read_isoformlist.close()

bk_junc_index=open(fname_junc_index)
sh_junc_index=bk_junc_index.readlines()
row_number_junc_index=len(sh_junc_index)
bk_junc_index.close()

#bk_junc_combination=open(fname_junc_combination)
#sh_junc_combination=bk_junc_combination.readlines()
#row_number_junc_combination=len(sh_junc_combination)
#bk_junc_combination.close()

bk_uniq_juncExpression=open(fname_uniq_juncExpression)
sh_uniq_juncExpression=bk_uniq_juncExpression.readlines()
row_number_uniq_juncExpression=len(sh_uniq_juncExpression)
bk_uniq_juncExpression.close()

bk_mapping_log=open(fname_mapping_log)
sh_mapping_log=bk_mapping_log.readlines()
row_number_mapping_log=len(sh_mapping_log)
bk_mapping_log.close()
############get isoform_gene pair and junc index##############
isoform_gene_pair={}
isoform_juncInformation={}
gene_strand={}
for v in range(1, row_number_junc_index):
    temp=sh_junc_index[v][0:-1].split('\t')
    transcript_id=temp[0]
    gene_id=temp[1]
    gene_strand[gene_id]=temp[2]
    isoform_gene_pair[transcript_id]=gene_id
    isoform_juncInformation[transcript_id]=[]
    junc_list=temp[4][2:-2].split('), (')
    for vv in range(0, len(junc_list)):
        temp1=junc_list[vv].split(', ')
        isoform_juncInformation[transcript_id].append(temp1)
############get length of each gene##############   
#gene_length={}
#gene_strand={}
#for v in range(1,row_number_junc_combination):
#    temp1=sh_junc_combination[v].split('\t')
#    gene_length[temp1[0]]=int(temp1[2])
#    gene_strand[temp1[0]]=temp1[1]
############get rpkm of each gene##############   
junc_uniq_rpkm={}
for v in range(1,row_number_uniq_juncExpression):
    temp1=sh_uniq_juncExpression[v][0:-1].split('\t')
    #JSH
    # check if the line is a comment
    if (not '#' in temp1[0]):
       junc_uniq_rpkm[temp1[0]]=float(temp1[3])
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
############Assign reads base on proportion of juncs##############
junc_express={} 
chromosome=''
for v in range(0, row_number_read_isoformlist):
    temp=sh_read_isoformlist[v][0:-1].split('\t')
    readID=temp[0]
    chromosome=temp[1]
    isoformlist=temp[2][0:-1].split(',')
    length=len(isoformlist)
    junc_proportion={}
    total=0.0
    for vv in range(0, length):
        temp2=isoformlist[vv].split(':') 
        isoform_name=temp2[0]
        map_start=int(temp2[1])
        map_end=map_start+read_length-1
        try:
            gene_name=isoform_gene_pair[isoform_name]
        except:
            #print "Warning:Transcript:%s is not found!" % (isoform_name)
            continue
        junc_list=isoform_juncInformation[isoform_name]
        for w in range(0, len(junc_list)):
            junc_location=junc_list[w]
            junc_start=int(junc_location[0])
            junc_end=junc_start+int(junc_location[3])-int(junc_location[2])
            junc_name=gene_name+':'+str(junc_location[2])+':'+str(junc_location[3])
            if map_end>=(junc_start+1) and map_start<=junc_start:
                try:
                    temp_rpkm=junc_uniq_rpkm[junc_name]
                except:
                    temp_rpkm=0
                junc_proportion[junc_name]=temp_rpkm    
                total+=temp_rpkm
    if total==0:
        total+=len(junc_proportion.keys())
        for e in junc_proportion:
            proportion=1/total
            try:
                junc_express[e]+=proportion
            except:
                junc_express[e]=proportion
    else:        
        for e in junc_proportion:
            proportion=junc_proportion[e]/total
            try:
                junc_express[e]+=proportion
            except:
                junc_express[e]=proportion

junc_read_ration={}    
for v in junc_express:
    read_ratio=(float(junc_express[v])/mapped_reads)*1000000
    junc_read_ration[v]=read_ratio

file_output=open(fname_output, 'w')
method_manual='# There are 4 methods to deal with multi-mapped reads:'+'\n'
method_manual+='# 1.Unique: remove the multi-mapped reads directly.'+'\n'
method_manual+='# 2.Random: assign the multi-mapped reads randomly.'+'\n'
method_manual+='# 3.Proportion: assign the multi-mapped reads according to the proportion of gene expression level.'+'\n'
method_manual+='# 4.Mappable Library: only count the reads falling into mappable regions, and it is based on the unmappable repetitive region library, which can be built automatically.'+'\n'
method_manual+='# You chose the 3rd one. This method assigns the multi-mapped reads according to the proportion of gene expression level.'+'\n'
#file_output.writelines(method_manual)
titles='#GeneID:JunctionStartLocation:JunctionEndLocation'+'\t'+'Chromosome'+'\t'+'Strand'+'\t'+'Reads_Number'+'\t'+'#Reads/#Mapped Reads(Million)'+'\n'
file_output.writelines(titles)
for v in junc_express:
    temp=v.split(':')
    geneName=":".join(temp[0:-2])
    insert_line=str(v)+'\t'+chromosome+'\t'+gene_strand[geneName]+'\t'+format(junc_express[v],'.2f')+'\t'+str(junc_read_ration[v])+'\n'
    file_output.writelines(insert_line)
file_output.close()
print 'JunctionExpressionLevel (Assign Reads Based on Proportion) Done for Chromosome %s' % (chromosome)
