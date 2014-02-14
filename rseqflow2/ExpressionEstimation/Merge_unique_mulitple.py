#!/usr/bin/env python
'''
Created on 2012-08-17

@author: linliu
Modified by J.Herstein 2013-04-05
'''
import os
import sys
import optparse
import decimal
import re

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-u", "--unique-expressionLevel", action = "store", type = "string", dest = "unique_expression_result",
		  help = "expression level result base on uniquely mapped reads in txt format(in gene, exon or junction level)")
parser.add_option("-m", "--multiple-expressionLevel", action = "store", type = "string", dest = "multiple_expression_result",
		  help = "expression level result base on multiple mapped reads in txt format(in gene, exon or junction level)")
parser.add_option("-t", "--file-type", action = "store", type = "string", dest = "type_of_inputfile",
                  help = "specify the type of input file:gene, exon or junction(the argument is gene, exon or junction)")
parser.add_option("-o", "--output-mergeResult", action = "store", type = "string", dest = "output_merge_result",
		  help = "merge unique result and multiple result into one txt file(in gene, exon or junction level)")	    	  
		  

(options, args) = parser.parse_args()

if (options.unique_expression_result is None or
    options.multiple_expression_result is None or
    options.type_of_inputfile is None or
    options.output_merge_result is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_unique_expression=options.unique_expression_result
fname_multiple_expression=options.multiple_expression_result
fname_output=options.output_merge_result
file_type=options.type_of_inputfile

if (file_type != 'gene') and (file_type != 'exon') and (file_type != 'junction'):
    print "Invalid type! The file type should be gene, exon or junction."
    sys.exit(1)

try:
    bk_unique_expression=open(fname_unique_expression)
except:
    print >>std.err, "File "+fname_unique_expression+" doesn't exist!"
    sys.exit(1)
sh_unique_expression=bk_unique_expression.readlines()
row_number_unique_expression=len(sh_unique_expression)
bk_unique_expression.close()

try:
    bk_multiple_expression=open(fname_multiple_expression)
except:
    print "File "+fname_multiple_expression+" doesn't exist!"
    sys.exit(1)
sh_multiple_expression=bk_multiple_expression.readlines()
row_number_multiple_expression=len(sh_multiple_expression)
bk_multiple_expression.close()

############get unique expression level##############
ID_chromosome={}
ID_strand={}
ID_read_number={}
ID_rpkm={}
for v in range(1, row_number_unique_expression):
    temp=sh_unique_expression[v][0:-1].split('\t')
    ID=temp[0]

    #JSH
    # Skip commented lines
    test=re.match('^#',ID)
    if test:
        continue

    chromosome=temp[1]
    strand=temp[2]
    read_number=float(temp[3])
    rpkm=float(temp[4])
    ID_chromosome[ID]=chromosome
    ID_strand[ID]=strand
    ID_read_number[ID]=read_number
    ID_rpkm[ID]=rpkm
###########merge unique and multiple expression level##########
for v in range(1, row_number_multiple_expression):
    temp=sh_multiple_expression[v][0:-1].split('\t')
    ID=temp[0]

    #JSH
    # Skip commented lines
    test=re.match('^#',ID)
    if test:
        continue

    chromosome=temp[1]
    strand=temp[2]
    read_number=float(temp[3])
    rpkm=float(temp[4])
    try:
        ID_read_number[ID]+=read_number
        ID_rpkm[ID]+=rpkm
        if chromosome!=ID_chromosome[ID] or strand!=ID_strand[ID]:
           print "Warning:chromosome or strand are not identical."
    except:
        ID_read_number[ID]=read_number
        ID_rpkm[ID]=rpkm
        ID_chromosome[ID]=chromosome
        ID_strand[ID]=strand    
################write to file#################
file_output=open(fname_output, 'w')
if file_type=='gene':
    titles='#GeneID'+'\t'+'Chromosome'+'\t'+'Strand'+'\t'+'Reads_Number'+'\t'+'RPKM'+'\n'
elif file_type=='exon':
    titles='#GeneID:ExonStartLocation:ExonEndLocation'+'\t'+'Chromosome'+'\t'+'Strand'+'\t'+'Reads_Number'+'\t'+'RPKM'+'\n'
elif file_type=='junction':
    titles='#GeneID:JunctionStartLocation:JunctionEndLocation'+'\t'+'Chromosome'+'\t'+'Strand'+'\t'+'Reads_Number'+'\t'+'#Reads/#Mapped Reads(Million)'+'\n'

file_output.writelines(titles)
for ID in ID_read_number:
    insert_line=str(ID)+'\t'+ID_chromosome[ID]+'\t'+ID_strand[ID]+'\t'+str(ID_read_number[ID])+'\t'+str(ID_rpkm[ID])+'\n'
    file_output.writelines(insert_line)
file_output.close()

print 'Merge Done for Input Files'
