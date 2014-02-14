#!/usr/bin/env python

'''
Created on 2012-07-26


@author: Liu Lin
'''

import os
import sys
import optparse
import copy
# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-g", "--gene-annotation", action = "store", type = "string", dest = "gene_annotation",
		  help = "gene annotation from UCSC")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output",
		  help = "output: exon combination by chromosome")

(options, args) = parser.parse_args()

if (options.gene_annotation is None or
    options.output is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname = options.gene_annotation
outputFilename = options.output

#fname=sys.argv[1]
#outputFilename=sys.argv[2]

#fname='GencodeAnnotationManual_chrY.gtf'
#outputFilename='singleGene_chrY.gtf'

try:
    bk=open(fname)
except:
    print prog_base + ": error: cannot open file " + fname
    sys.exit(1)
sh=bk.readlines()
row_number=len(sh)
bk.close()

try:
    ws=open(outputFilename,'w')
except:
    print prog_base + ": error: cannot open file " + outputFilename
    sys.exit(1)

genename_repeat=[]
genename=[]
whole_info=[]

##############################################################################
#-------------------get gene list in chromosome------------------------------#

#************************initial*******************#
gtfLine=[]
gene_appear_count=0
pre_transcript_id="none"
temp=sh[0][0:-1].split('\t')
if len(temp) != 9:
    print "warning:gtf format should have 9 column!"
    sys.exit(1)
temp1					=temp[8]
attributes		=temp1.split(';')
temp2					=attributes[0].split('"')
pre_gene_id		=temp2[1]

for v in range(0,row_number):
    temp=sh[v][0:-1].split('\t')
    if len(temp) != 9:
        print "warning:gtf format should have 9 column!"
        sys.exit(1)
    chr			=temp[0]
    type		=temp[2]
    start		=int(temp[3])-1
    end			=int(temp[4])
    strand	=temp[6]
    temp1		=temp[8]
    attributes		=temp1.split(';')
    temp2					=attributes[0].split('"')
    gene_id				=temp2[1]
    temp3					=attributes[1].split('"')
    transcript_id	=temp3[1]
    if gene_id!=pre_gene_id and transcript_id!=pre_transcript_id:
        if gene_appear_count==1:
            for vv in range(0,len(gtfLine)):
                ws.writelines(gtfLine[vv])
        gene_appear_count = 1
        gtfLine=[]
        gtfLine.append(sh[v])
    if gene_id==pre_gene_id and transcript_id!=pre_transcript_id:
        gene_appear_count += 1
        gtfLine=[]
        gtfLine.append(sh[v])
    if gene_id==pre_gene_id and transcript_id==pre_transcript_id:
        gtfLine.append(sh[v])
    pre_gene_id=gene_id
    pre_transcript_id=transcript_id     
#**************************save the last transcript information**************    
if gene_appear_count==1:
    for vv in range(0,len(gtfLine)):
        ws.writelines(gtfLine[vv])
ws.close()
###################################################################################        
