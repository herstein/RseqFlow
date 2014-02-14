#!/usr/bin/env python
'''
Created on 2012-08-21

@author: linliu
'''

import os,sys
import re
import string
import optparse
import warnings
import string
import collections
import math
import sets
from time import strftime

#from bx.bitset import *
#from bx.bitset_builders import *
from bx.intervals import *

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-s", "--sam-file", action = "store", type = "string", dest = "sam_alignment_file",
		  help = "Input alignment file in SAM format")
parser.add_option("-a", "--annotation-file", action = "store", type = "string", dest = "annotation_bed_file",
		  help = "Reference gene model in bed fomat.")
parser.add_option("-n", "--number-ofSample", action = "store", type = "string", dest = "number_of_sample", default=200000,
		  help = "Number of reads sampled from SAM/BAM file. default=%default")
parser.add_option("-o", "--output-prefix", action = "store", type = "string", dest = "output_prefix",
		  help = "refix of output file(s). [required]")	  

(options, args) = parser.parse_args()

if (options.sam_alignment_file is None or
    options.annotation_bed_file is None or
    options.output_prefix is None ):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)
    
fname_sam=options.sam_alignment_file
fname_annotation=options.annotation_bed_file
sample_size=options.number_of_sample


if sample_size<1000:
    print "Warning: Sample Size too small to give a accurate estimation"
    
bk_sam=open(fname_sam)
#sh_sam=bk_sam.readlines()
#row_number_sam=len(sh_sam)
#bk_sam.close()

bk_annotation=open(fname_annotation)
sh_annotation=bk_annotation.readlines()
row_number_annotation=len(sh_annotation)
bk_annotation.close()

##########################################read annotation##########################################
p_strandness={}
s_strandness={}
chrRange_strand={}
p_strandness=collections.defaultdict(int)
s_strandness=collections.defaultdict(int)
print "Reading annotation...."

for v in range(0, row_number_annotation):
    fields=sh_annotation[v].split('\t')
    if len(fields)!=12:
        print >>sys.stderr, "\n File: %s , line %d " % (fname_annotation, v+1)
        print >>sys.stderr, "  Bed file should be 12-column, so skip this line."
        continue
    chromosome= fields[0]
    tx_start  = int(fields[1])+1
    tx_end    = int(fields[2])
    geneName  = fields[3]
    strand    = fields[5]
		
    if chromosome not in chrRange_strand:
        chrRange_strand[chromosome]=Intersecter()
    chrRange_strand[chromosome].insert(tx_start, tx_end, strand) 
##########################################read sam file##########################################
read_paired=1
read_mapped=2
read_unmapped=4
mate_unmapped=8
read_reverse=16
mate_reverse=32
is_read1=64
is_read2=128
secondary=256
qc_fail=512
duplicate=1024
count=0
print "Reading sam file..."
#########get read length########
for v in bk_sam.xreadlines():
    fields=v[0:-1].split('\t')
    cigar=fields[5]
    temp=cigar.split('M')
    if len(temp)==2:
        read_length=int(temp[0])
        break    
bk_sam.seek(0)        
#JSH
lineNum = 0
badLineCount = 0
for v in bk_sam.xreadlines():
    #JSH
    lineNum=lineNum+1
    if count >=sample_size:
        break
    fields=v[0:-1].split('\t')
    flag=int(fields[1])
    if (flag&read_unmapped or 
        flag&secondary or
        flag&qc_fail or
        flag&duplicate):
	badLineCount=badLineCount+1
	#print >>sys.stderr, "\n File: %s: line number %d" % (fname_sam, lineNum)
        #print >>sys.stderr, "  The read fails to pass the criterion, so skip this line."
        continue
       
    chromosome=fields[2]
    if flag&read_paired:
        if flag&is_read1:
            read_id='1'
        elif flag&is_read2:
            read_id='2'
        else:
            print >>sys.stderr, "\n File: %s , line %d" % (fname_sam, v+1)
            print >>sys.stderr, "  No flag tells if the read is 1 or 2, so skip this line."
            continue
        if flag&read_reverse:
	    map_strand='-'
        else:
	    map_strand='+'
        map_start=int(fields[3])
	map_end=map_start+read_length-1
	if chromosome in chrRange_strand:
	    if len(set(chrRange_strand[chromosome].find(map_start,map_end)))==1:
		gene_strand=set(chrRange_strand[chromosome].find(map_start,map_end)).pop()
                s_strandness[read_id+map_strand + gene_strand]+=1  
		#JSH
		count += 1
    else:
	if flag&read_reverse:
	    map_strand='-'
	else:
            map_strand='+'
	map_start=int(fields[3])
	map_end=map_start+read_length-1
	if chromosome in chrRange_strand:
	    if len(set(chrRange_strand[chromosome].find(map_start,map_end)))==1:
		gene_strand=set(chrRange_strand[chromosome].find(map_start,map_end)).pop()
                s_strandness[map_strand + gene_strand]+=1
		count += 1

bk_sam.close()
##########################################caculate strand specificity##########################################							
print "Total number reads: %d" % lineNum
print "Total " + str(count) + " usable reads were sampled"
print "Number reads that failed to pass criteria: %d" % badLineCount
protocol="unknown"
strandness=None
spec1_plus=0.0
spec1_redu=0.0
spec2_plus=0.0
spec2_redu=0.0
other=0.0
if len(p_strandness) >0 and len(s_strandness) ==0 :
			protocol="PairEnd"
			#for k,v in p_strandness.items():
			#	print >>sys.stderr, k + '\t' + str(v)
			spec1_plus = (p_strandness['1++'] + p_strandness['2+-'])/float(sum(p_strandness.values()))
			spec1_redu = (p_strandness['1--'] + p_strandness['2-+'])/float(sum(p_strandness.values()))
			spec2_plus = (p_strandness['1+-'] + p_strandness['2++'])/float(sum(p_strandness.values()))
			spec2_redu = (p_strandness['1-+'] + p_strandness['2--'])/float(sum(p_strandness.values()))
			#spec1= (p_strandness['1++'] + p_strandness['1--'] + p_strandness['2+-'] + p_strandness['2-+'])/float(sum(p_strandness.values()))
			#spec2= (p_strandness['1+-'] + p_strandness['1-+'] + p_strandness['2++'] + p_strandness['2--'])/float(sum(p_strandness.values()))
			other = 1-spec1_plus-spec1_redu-spec2_plus-spec2_redu
			
elif len(s_strandness) >0 and len(p_strandness) ==0 :
			protocol="SingleEnd"
			#for k,v in s_strandness.items():
			#	print  >>sys.stderr, k + '\t' + str(v)
			spec1_plus = s_strandness['++']/float(sum(s_strandness.values()))
			spec1_redu = s_strandness['--']/float(sum(s_strandness.values()))
			spec2_plus = s_strandness['+-']/float(sum(s_strandness.values()))
			spec2_redu = s_strandness['-+']/float(sum(s_strandness.values()))
			#spec1 = (s_strandness['++'] + s_strandness['--'])/float(sum(s_strandness.values()))
			#spec2 = (s_strandness['+-'] + s_strandness['-+'])/float(sum(s_strandness.values()))
			other = 1-spec1_plus-spec1_redu-spec2_plus-spec2_redu
else:
			protocol="Mixture"
			spec1 = "NA"
			spec2 = "NA"
			other = "NA"			
##########################################write to ouput file##########################################			
if other <0: other=0.0
OUT=open(options.output_prefix+".strand_stat.txt",'w')
print >>OUT, "#=========================================="
if protocol == "PairEnd":
		print >>OUT,"This is PairEnd Data"
		print >>OUT,"Fraction of reads explained by \"1++,2+-\": %.4f" % spec1_plus
		print >>OUT,"Fraction of reads explained by \"1--,2-+\": %.4f" % spec1_redu
		print >>OUT,"Fraction of reads explained by \"1+-,2++\": %.4f" % spec2_plus
		print >>OUT,"Fraction of reads explained by \"1-+,2--\": %.4f" % spec2_redu
		#print "Fraction of reads explained by \"1++,1--,2+-,2-+\": %.4f" % sp1
		#print "Fraction of reads explained by \"1+-,1-+,2++,2--\": %.4f" % sp2
	        #print >>OUT,"Fraction of reads explained by other combinations: %.4f" % other"
		print >>OUT, "#=========================================="
		print >>OUT, "Note:"
		print >>OUT, "1++,2+-:read1 mapped to '+' strand indicates parental gene on '+' strand"
		print >>OUT, "        read2 mapped to '+' strand indicates parental gene on '-' strand"
		print >>OUT, "1--,2-+:read1 mapped to '-' strand indicates parental gene on '-' strand"
		print >>OUT, "        read2 mapped to '-' strand indicates parental gene on '+' strand"
		print >>OUT, "1+-,2++:read1 mapped to '+' strand indicates parental gene on '-' strand"
		print >>OUT, "        read2 mapped to '+' strand indicates parental gene on '+' strand"
		print >>OUT, "1-+,2--:read1 mapped to '-' strand indicates parental gene on '+' strand"
		print >>OUT, "        read2 mapped to '-' strand indicates parental gene on '-' strand"
elif protocol == "SingleEnd":
		print >>OUT,"This is SingleEnd Data"
		print >>OUT,"Fraction of reads explained by \"++\": %.4f" % spec1_plus
		print >>OUT,"Fraction of reads explained by \"--\": %.4f" % spec1_redu
		print >>OUT,"Fraction of reads explained by \"+-\": %.4f" % spec2_plus
		print >>OUT,"Fraction of reads explained by \"-+\": %.4f" % spec2_redu
		#print "Fraction of reads explained by \"++,--\": %.4f" % spec1
		#print "Fraction of reads explained by \"+-,-+\": %.4f" % sp2
		#print >>OUT,"Fraction of reads explained by other combinations: %.4f" % other"
		print >>OUT, "#=========================================="
		print >>OUT, "Note:"
		print >>OUT, "++:read mapped to '+' strand indicates parental gene on '+' strand"
		print >>OUT, "--:read mapped to '-' strand indicates parental gene on '-' strand"
		print >>OUT, "+-:read mapped to '+' strand indicates parental gene on '-' strand"
		print >>OUT, "-+:read mapped to '-' strand indicates parental gene on '+' strand"
else:
		print >>OUT,"Unknown Data type"
		print >>OUT, "#=========================================="
OUT.close()		
