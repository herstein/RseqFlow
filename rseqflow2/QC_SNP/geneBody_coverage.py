#!/usr/bin/env python
'''
Created on 2012-08-22

@author: linliu
'''

import os,sys
import re
import string
import optparse
import warnings
import collections
import math
import sets
import subprocess
#from time import strftime

from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

def percentile_list(N):
    """
    Find the percentile of a list of values.
    @parameter N - is a list of values. Note N MUST BE already sorted.
    @return - the list of percentile of the values
    """
    if not N:return None
    per_list=[]
    for i in range(0,101):
	k = (len(N)-1) * i/100.0
	f = math.floor(k)
	c = math.ceil(k)
	if f == c:
	    per_list.append( int(N[int(k)])  )
	else:
	    d0 = N[int(f)] * (c-k)
	    d1 = N[int(c)] * (k-f)
	    per_list.append(int(round(d0+d1)))	
    return per_list

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-s", "--sam-file", action = "store", type = "string", dest = "sam_alignment_file",
		  help = "Input alignment file in SAM format")
parser.add_option("-a", "--annotation-file", action = "store", type = "string", dest = "annotation_bed_file",
		  help = "Reference gene model in bed fomat.")
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

bk_sam=open(fname_sam)
#sh_sam=bk_sam.readlines()
#row_number_sam=len(sh_sam)
#bk_sam.close()

bk_annotation=open(fname_annotation)
#sh_annotation=bk_annotation.readlines()
#row_number_annotation=len(sh_annotation)
#bk_annotation.close()

outfile1 = options.output_prefix + ".geneBodyCoverage_plot.r"
outfile2 = options.output_prefix + ".geneBodyCoverage.txt"
OUT1 = open(outfile1,'w')
OUT2 = open(outfile2,'w')

ranges={}
totalReads=0
fragment_num=0		#splice reads will counted twice
rpkm={}
		
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
_cigar_split=re.compile(r'(\d+)[M|N]')
print "reading "+ fname_sam + '...',
for v in bk_sam.xreadlines():
    if v.startswith("@"):continue
    fields=v[0:-1].split('\t')
    flag=int(fields[1])
    if flag&read_unmapped: continue		#skip unmap reads
    totalReads +=1
			
    chrom = fields[2].upper()
    chromStart = string.atoi(fields[3])-1
    comb=[int(i) for i in _cigar_split.findall(fields[5])]	#"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
    fragment_num += (len(comb) +1)/2
    blockStart=[]
    blockSize=[]
			
    for i in range(0,len(comb),2):
	blockStart.append(chromStart + sum(comb[:i]) )
				
    for i in range(0,len(comb),2):
	blockSize.append(comb[i])
			
    for st,size in zip(blockStart,blockSize):
	if chrom not in ranges:ranges[chrom] = Intersecter()
	else:ranges[chrom].add_interval( Interval( st, st+size ) )  
bk_sam.close()	
##########################################calculating coverage over gene body##########################################
print "calculating coverage over gene body ..."
coverage=collections.defaultdict(int)
flag=0
for v in bk_annotation.xreadlines():
    try:
	if v.startswith(('#','track','browser')):continue  
        # Parse fields from gene tabls
	fields 	= v.split()
	chrom   = fields[0].upper()
	tx_start  = int( fields[1] )
	tx_end    = int( fields[2] )
	geneName  = fields[3]
	strand    = fields[5]
				
	exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
	exon_starts = map((lambda x: x + tx_start ), exon_starts)
	exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
	exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
    except:
	print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
	continue
    gene_all_base=[]
    percentile_base=[]
    mRNA_len =0
    flag=0
    for st,end in zip(exon_starts,exon_ends):
	gene_all_base.extend(range(st+1,end+1))		#0-based coordinates on genome
	mRNA_len = len(gene_all_base)
	if mRNA_len <100:
	    flag=1
	    break
    if flag==1: continue
    if strand == '-':gene_all_base.sort(reverse=True)			#deal with gene on minus stand
    else:gene_all_base.sort(reverse=False)
    percentile_base = percentile_list (gene_all_base)	#get 101 points from each gene's coordinates
			
    for i in range(0,len(percentile_base)):
	if chrom in ranges:
	    coverage[i] += len(ranges[chrom].find(percentile_base[i], percentile_base[i]+1))
bk_annotation.close()
x_coord=[]
y_coord=[]
print >>OUT2, "Total reads: " + str(totalReads)
print >>OUT2, "Fragment number: " + str(fragment_num)
print >>OUT2, "percentile\tcount"
for i in coverage:
    x_coord.append(str(i))
    y_coord.append(str(coverage[i]))
    print >>OUT2, str(i) + '\t' + str(coverage[i])
print >>OUT1, 'pdf(\"%s\")' % (options.output_prefix+".geneBody_coverage.pdf")
print >>OUT1, "x=0:100"
print >>OUT1, "y=c(" + ','.join(y_coord) + ')'
print >>OUT1, "plot(x,y,xlab=\"percentile of gene body (5'->3')\",ylab='read number',type='s')"
print >>OUT1, "dev.off()"

OUT1.close()
OUT2.close()
try:
    subprocess.call("Rscript " + options.output_prefix + '.geneBodyCoverage_plot.r',shell=True)
except:
    print >>sys.stderr, "Cannot generate pdf file from " + options.output_prefix + '.geneBodyCoverage_plot.r'
    pass
