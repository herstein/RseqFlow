#!/usr/bin/env python
'''
Created on 2012-07-30

@author: Liu Lin 
Modified by J.Herstein 2013-06-04 to allow fastq.gz format
'''

import os
import sys
import optparse
import copy

#JSH
import gzip
import subprocess

# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-s", "--sam-alignment", action = "store", type = "string", dest = "alignment_sam_file",
		  help = "sam file of alignment")
parser.add_option("-q", "--fq-reads", action = "store", type = "string", dest = "reads_fq_file",
		  help = "fq file of reads")
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
outputFilename = options.output+".mapping_report.txt"

#fname=sys.argv[1]
#outputFilename=sys.argv[2]

#fname='c:/Python26/test file/GencodeAnnotationManual_chrY.txt'
#outputFilename='GencodeExonCombination_chrY.txt'

try:
    #JSH 2013-06-04
    #bk_fq=open(fname_fq)
    if fname_fq.endswith(".gz"):
         bk=subprocess.Popen(["zcat", fname_fq], stdout=subprocess.PIPE,close_fds=True)
         bk_fq=bk.stdout
    else:
         bk=open(fname_fq)
         bk_fq=bk.xreadlines()
except:
    print prog_base + ": error: cannot open file " + fname_fq
    sys.exit(1)
#sh_fq=bk_fq.readlines()
#row_number_fq=len(sh_fq)
#bk_fq.close()
row_number_fq=0
#JSH
#for v in bk_fq.xreadlines():
for v in bk_fq:
	row_number_fq+=1
#JSH
if not fname_fq.endswith(".gz"):
    bk_fq.close()
total_reads=row_number_fq/4

try:
    bk_sam=open(fname_sam)
except:
    print prog_base + ": error: cannot open file " + fname_sam
    sys.exit(1)
#sh_sam=bk_sam.readlines()
#row_number_sam=len(sh_sam)
#bk_sam.close()

##########################################################################################
#########################SingleEnd or PairEnd?############################################
temp=bk_sam.readline()
temp1=temp[0:-1].split('\t')
flag=int(temp1[1])
bk_sam.seek(0)
if (flag%2)==0:
    type="singleEnd"
else:
    type="pairEnd"    
########initial########
total_map_reads=0
unique_reads=0
multiple_reads=0
#######################start to make the report##########################################
##############################################################################
#-------------------get reads information from sam file------------------------------#
####------------------------------singleEnd------------------------------------######
if type=="singleEnd":
    readname_repeat=[]
    read_count=0
    #read_list=[]
    for v in bk_sam.xreadlines():
        if v[0]!='@':
            temp=v[0:-1].split('\t')
            readID=temp[0]
            readname_repeat.append(readID)

    readname_repeat.sort()
    pre_readID='none'
    repeat_number=len(readname_repeat)

    for v in range(0,repeat_number):
        readID=readname_repeat[v]
        if readID!=pre_readID:
            if read_count==1:
                unique_reads+=1
            elif read_count>1:
                multiple_reads+=1
            read_count=1
            total_map_reads+=1
            pre_readID=readID
        else:
            read_count+=1

    if read_count==1:
        unique_reads+=1
    elif read_count>1:
        multiple_reads+=1
####-----------------------------pairEnd------------------------------------######        
else:
    readname_repeat=[]
    read_count=0
    for v in bk_sam.xreadlines():
        if v[0]!='@':
            temp=v[0:-1].split('\t')
            readID=temp[0]
            mapping_type=temp[6]
            if mapping_type=='=':
                readname_repeat.append(readID)
            #if readID not in readID_twoType_count:
            #    readID_twoType_count[readID]=[0,0]
            #if mapping_type == '=':
            #    readID_twoType_count[readID][0]+=1
            #elif mapping_type == '*':
            #    readID_twoType_count[readID][1]+=1   
    readname_repeat.sort()
    pre_readID='none'
    repeat_number=len(readname_repeat)

    for v in range(0,repeat_number):
        readID=readname_repeat[v]
        if readID!=pre_readID:
            if read_count==2:
                unique_reads+=1
            elif read_count>2:
                multiple_reads+=1
            read_count=1
            total_map_reads+=1
            pre_readID=readID
        else:
            read_count+=1

    if read_count==2:
        unique_reads+=1
    elif read_count>2:
        multiple_reads+=1      

bk_sam.close()       
    #for readID in readID_twoType_count:
    #    read_count=readID_twoType_count[readID][0]/2+readID_twoType_count[readID][1]
    #    total_map_reads+=1
    #    if read_count==1:
    #        unique_reads+=1
    #    elif read_count>1:
    #        multiple_reads+=1        
##############################################################################
       
##############################################################################

######write to file

try:
    OUT=open(outputFilename,'w')
except:
    print prog_base + ": error: cannot open file " + outputFilename
    sys.exit(1)
if type=="singleEnd":
    print >>OUT,"\nThis is Single End Data."
else:
    print >>OUT,"\nThis is Paired End Data.(If reads cannot be pair mapped, the reads are considered to be unmapped.)"    
print >>OUT,"#=================================================="
print >>OUT,"%-30s\t%s" % ("Number of Total Reads:",format(total_reads,','))
#print >>OUT,"%-30s\t%s\t%.2f%%" % ("Number of Mapped Reads:",format(total_map_reads,','),total_map_reads*100.00/total_reads)
print >>OUT,"%-30s\t%s\t%.2f%%" % ("Number of Mapped Reads:",format(total_map_reads,','),total_map_reads*100.00/total_reads)
print >>OUT,"%-30s%s\t%.2f%%" % ("Number of Multiple-mapped reads:",format(multiple_reads,','),multiple_reads*100.00/total_reads)
print >>OUT,"%-30s\t%s\t%.2f%%" % ("Number of Unique-mapped reads:",format(unique_reads,','),unique_reads*100.00/total_reads)
print >>OUT,"#=================================================="
OUT.close()
