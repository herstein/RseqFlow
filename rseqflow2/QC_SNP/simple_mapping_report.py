#!/usr/bin/env python
'''
Created on 2012-07-30

@author: Liu Lin 
Modified by J.Herstein 2013-06-04 to allow fastq.gz format
Modified by R. Mayani  2013-08-30 to print title at beginning
'''

import os
import sys
import optparse
import copy

#JSH
import gzip
import subprocess
import re

# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-s", "--sam-alignment", action = "store", type = "string", dest = "alignment_sam_file",
		  help = "sam file of alignment")
parser.add_option("-q", "--fq-reads", action = "store", type = "string", dest = "reads_fq_file",
		  help = "fq file of reads")
parser.add_option("-t", "--title", action = "store", type = "string", dest = "title",
		  default='', help = "Title to add to the report")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output",
		  help = "output: mapping report")

(options, args) = parser.parse_args()

if (options.alignment_sam_file is None or
    options.reads_fq_file is None or
    options.output is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_sam = options.alignment_sam_file
fname_fq  = options.reads_fq_file
outputFilename = options.output+".simple_mapping_report.txt"
title = options.title
#fname=sys.argv[1]
#outputFilename=sys.argv[2]

#fname='c:/Python26/test file/GencodeAnnotationManual_chrY.txt'
#outputFilename='GencodeExonCombination_chrY.txt'
###############################get total reads###################################
try:
    #JSH
    #bk_fq=open(fname_fq)
    if fname_fq.endswith(".gz"):
          bk=subprocess.Popen(["zcat", fname_fq], stdout=subprocess.PIPE,close_fds=True)
          bk_fq=bk.stdout
    else:
          bk_fq=open(fname_fq)

except:
    print prog_base + ": error: cannot open file " + fname_fq
    sys.exit(1)

#JSH
#sh_fq=bk_fq.readlines()
#row_number_fq=len(sh_fq)
row_number_fq=0
for line in bk_fq.xreadlines():
     row_number_fq=row_number_fq+1

#JSH
if not fname_fq.endswith(".gz"):
    bk_fq.close()

total_reads=row_number_fq/4
###############################get mapped reads###################################
try:
    bk_sam=open(fname_sam)
except:
    print prog_base + ": error: cannot open file " + fname_sam
    sys.exit(1)
#JSH
#sh_sam=bk_sam.readlines()
#row_number_sam=len(sh_sam)
#bk_sam.close()

#readID_repeat_list=[]
#readID_list=[]
read_count=0

#JSH
for line in bk_sam.xreadlines():
       fields=line[0:-1].split('\t')
       chrom=fields[2]
       # Count as read if column 3 is not '*'
       test=re.match('^\*$',chrom)
       if test:
	  continue
       read_count=read_count + 1
bk_sam.close()
mapped_reads=read_count

#JSH
#for v in range(0,row_number_sam):
#    if sh_sam[v][0]=='@':
#    	  continue
#    temp=sh_sam[v][0:-1].split('\t')
#    readID=temp[0]
#    readID_repeat_list.append(readID)
#print "done"
#readID_list=list(set(readID_repeat_list))
#mapped_reads=len(readID_list)
#############read log file##################
OUT=open(outputFilename, 'w')
if title:
       print >>OUT, "Title: %s" % title

print >>OUT, "#=================================================="
print >>OUT, "%s\t%s" % ('Number of Reads:', format(total_reads,','))
print >>OUT, "%s\t%s" % ('Number of Mapped Reads:', format(mapped_reads,','))
print >>OUT, "%s\t\t%.2f%%" % ('Mapping Rate:', mapped_reads*100.00/total_reads)
print >>OUT, "#=================================================="
OUT.close()
