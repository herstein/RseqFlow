#!/usr/bin/env python

'''
Created on 2012-08-15

@author: linliu
'''
import os
import sys
import optparse
import random

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-m", "--reads-multiIsoformList", action = "store", type = "string", dest = "reads_multiIsoform_list",
		  help = "Transcript list of multiply-mapped reads")
parser.add_option("-i", "--junction-index", action = "store", type = "string", dest = "junction_index_file",
		  help = "txt file of junction index")		  
#parser.add_option("-c", "--junction-combination", action = "store", type = "string", dest = "junction_combination_file",
#		  help = "txt file of junction combination")		  
parser.add_option("-l", "--mapping-log", action = "store", type = "string", dest = "mapping_log_file",
		  help = "log file which containts information about mapped reads and read length")		  
parser.add_option("-j", "--gene-expression", action = "store", type = "string", dest = "junction_expression_level",
		  help = "output file of junction expression level")	    	  
		  

(options, args) = parser.parse_args()

if (options.reads_multiIsoform_list is None or
    options.junction_index_file is None or
    options.mapping_log_file is None or
    options.junction_expression_level is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_read_isoformlist=options.reads_multiIsoform_list
fname_junction_index=options.junction_index_file
#fname_junction_combination=options.junction_combination_file
fname_mapping_log=options.mapping_log_file
fname_output=options.junction_expression_level

bk_read_isoformlist=open(fname_read_isoformlist)
sh_read_isoformlist=bk_read_isoformlist.readlines()
row_number_read_isoformlist=len(sh_read_isoformlist)
bk_read_isoformlist.close()

bk_junction_index=open(fname_junction_index)
sh_junction_index=bk_junction_index.readlines()
row_number_junction_index=len(sh_junction_index)
bk_junction_index.close()

#bk_junction_combination=open(fname_junction_combination)
#sh_junction_combination=bk_junction_combination.readlines()
#row_number_junction_combination=len(sh_junction_combination)
#bk_junction_combination.close()

bk_mapping_log=open(fname_mapping_log)
sh_mapping_log=bk_mapping_log.readlines()
row_number_mapping_log=len(sh_mapping_log)
bk_mapping_log.close()
############get isoform_gene pair##############
isoform_gene_pair={}
isoform_jcInformation={}
gene_strand={}
for v in range(1, row_number_junction_index):
    temp=sh_junction_index[v][0:-1].split('\t')
    transcript_id=temp[0]
    gene_id=temp[1]
    gene_strand[gene_id]=temp[2]
    isoform_gene_pair[transcript_id]=gene_id
    isoform_jcInformation[transcript_id]=[]
    junction_list=temp[4][2:-2].split('), (')
    for vv in range(0, len(junction_list)):
        temp1=junction_list[vv].split(', ')
        isoform_jcInformation[transcript_id].append(temp1)
############get length of each gene##############   
#gene_length={}
#gene_strand={}
#for v in range(1,row_number_junction_combination):
#    temp1=sh_junction_combination[v].split('\t')
#    gene_length[temp1[0]]=int(temp1[2])
#    gene_strand[temp1[0]]=temp1[1]
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
junc_express={}  
chromosome=''
for v in range(0, row_number_read_isoformlist):
    temp=sh_read_isoformlist[v][0:-1].split('\t')
    readID=temp[0]
    chromosome=temp[1]
    isoformlist=temp[2][0:-1].split(',')
    length=len(isoformlist)
    junc_cordinate=[]
    for t in range(0,length):
        temp2=isoformlist[t].split(':')
        isoform_name=temp2[0]
        map_start=int(temp2[1])
        map_end=map_start+read_length-1
        try:
            gene_name=isoform_gene_pair[isoform_name]
        except:
            #print "Warning:Transcript:%s is not found!" % (isoform_name)
            #print "        You can consider to ignore this warning when calculating junction expression level."
            #print "        Maybe transcript:%s has no junction, just with single exon." % (isoform_name)
            continue
        junction_list=isoform_jcInformation[isoform_name]    
        for vv in range(0,len(junction_list)):
            junction_location=junction_list[vv]
            junction_start=int(junction_location[0])
            junction_name=gene_name+':'+str(junction_location[2])+':'+str(junction_location[3])
            if map_end>=(junction_start+1) and map_start<=junction_start:
                junc_cordinate.append(junction_name)
    junc_length=len(junc_cordinate)
    if junc_length==0:
        continue
    randomly=random.randint(0,junc_length-1)
    try:
        junc_express[junc_cordinate[randomly]]+=1
    except:
        junc_express[junc_cordinate[randomly]]=1         
    del randomly
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
method_manual+='# You chose the 2nd one. This method assigns the multi-mapped reads randomly..'+'\n'
#file_output.writelines(method_manual)
titles='#GeneID:JunctionStartLocation:JunctionEndLocation'+'\t'+'Chromosome'+'\t'+'Strand'+'\t'+'Reads_Number'+'\t'+'#Reads/#Mapped Reads(Million)'+'\n'
file_output.writelines(titles)
for v in junc_express:
    temp=v.split(':')
    geneName=":".join(temp[0:-2])
    insert_line=str(v)+'\t'+chromosome+'\t'+gene_strand[geneName]+'\t'+str(junc_express[v])+'\t'+str(junc_read_ration[v])+'\n'
    file_output.writelines(insert_line)
file_output.close()
print 'JunctionExpressionLevel (Assign Reads Randomly) Done for Chromosome %s' % (chromosome)
