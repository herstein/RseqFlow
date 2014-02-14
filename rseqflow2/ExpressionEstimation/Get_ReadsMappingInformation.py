#!/usr/bin/env python

#if a read mapped to multiple genes, only keep the reads mapped to single gene
#and then split the multi-mapped reads into the several genes evenly
#or assign the multi-mapped reads to one of the genes randomly

import os
import sys
import random
import optparse
import re

#fname_sam='c:/Python26/sample.txt'
#fname_table='c:/Python26/isoform_gene_table_ref19.txt'
#output='MappingSingleGene.txt'
#output2='SplitbyGenelength.txt'
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-s", "--sam-mappingResult", action = "store", type = "string", dest = "mappingResult_sam_file",
		  help = "mapping result in sam format")
parser.add_option("-l", "--output-log", action = "store", type = "string", dest = "mappingInfo_log_output",
		  help = "mapping infomation in log file")	    	  
		  

(options, args) = parser.parse_args()

if (options.mappingResult_sam_file is None or
    options.mappingInfo_log_output is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_sam=options.mappingResult_sam_file # the original mapping file
output=options.mappingInfo_log_output#
###########calculate mapped reads########
bk_sam=open(fname_sam)
#sh_sam=bk_sam.readlines()
#bk_sam.close()

read_mappingCount={}
samline=0
read_count=0
read_length=0
for line in bk_sam.xreadlines():
  samline+=1
  if line[0]=='@':
    #print "This is header of sam file."
    continue
  temp1=line[0:-1].split('\t')
  readID=temp1[0]
  #JSH 
  # Check that reference field is not "*"
  if re.match('^\*',temp1[2]):
      continue
  temp2=temp1[2].split('_',2)
  temp3=temp2[2].split('=')
  isoform=temp3[0]
  temp4=temp3[1].split(':')
  chromosome=temp4[0]
  #JSH
  # Check for valid read length
  if re.search(r'M$',temp1[5]):
    temp6=temp1[5].split('M')
    read_length=int(temp6[0])
  if readID not in read_mappingCount:
    read_mappingCount[readID]=0
    read_count+=1
bk_sam.close()
############get read length and chromosome##########
#temp=sh_sam[-1].split('\t')
#temp1=temp[5].split('M')
#read_length=int(temp1[0])
#temp2=temp[2].split('=')
#temp3=temp2[1].split(':')
#chromosome=temp3[0]
		
file_output=open(output, 'w')

print "Number of mapped reads:%d" % (read_count)
print "Read Length:%d" % (read_length)
print >> file_output, "Number of mapped reads:%d" % (read_count)
print >>file_output,  "Read Length:%d" % (read_length)
file_output.close()

print "Uniquely mapped result and multiple mapped result Done for Current Chromosome"
