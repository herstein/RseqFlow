#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
Calculating Phred Quality Score for each position on read. Note that each read should have 
the fixed (same) length

Modified by J.Herstein 2013-06-04 to allow for fastq.gz files
-------------------------------------------------------------------------------------------------'''

#import built-in modules
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
import copy
import subprocess

#JSH
import gzip
import subprocess


#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *


#changing history to this module



__copyright__ = "Copyright 2010, Wei Li's Lab"
__credits__ = []
__license__ = "GPL"
__version__="2.3"
__maintainer__ = "Liguo Wang"
__email__ = "liguow@bcm.edu"
__status__ = "Development" #Prototype or Production
prog_base = os.path.split(sys.argv[0])[1]
parser = optparse.OptionParser()

def printlog (mesg):
    '''print progress into stderr and log file'''
    mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
    LOG=open('class.log','a')
    print >>sys.stderr,mesg
    print >>LOG,mesg

def readsQual_boxplot(inputfile,outfile,shrink=1000):
    '''calculate phred quality score for each base in read (5->3)'''

    output = outfile + ".read_qual.r"
    FO=open(output,'w')
	

    quality = collections.defaultdict(dict)	#read_pos=>quality score=>count
    q_max = -1
    q_min = 10000
    q_list=[]
    i_box={}	#key is read postion,value is 
    ##################read input fq file##########
    #JSH
    if inputfile.endswith(".gz"):
        print "Load reads file(fastq.gz format)..."
    else:
	print "Load reads file(fastq format)..."
    try:
	#JSH 2013-06-04
    	#bk=open(inputfile)
	if inputfile.endswith(".gz"):
            bk=subprocess.Popen(["zcat", inputfile], stdout=subprocess.PIPE,close_fds=True)
            bk_reads=bk.stdout
        else:
            bk=open(inputfile)
            bk_reads=bk.xreadlines()

    except:
    	print prog_base + ": error: cannot open file " + inputfile
        sys.exit(1)
    #sh=bk.readlines()
    #row_number=len(sh)
    #bk.close()
    #################get reads quality############
    print "Reading "+inputfile+"... "
    count=0
    #JSH
    #for v in bk.xreadlines():
    for v in bk_reads:
	count+=1
	if (count%4)!=0:
	    continue
	temp=v[0:-1]
	qual_str = temp
	read_len = len(temp)
	for i,j in enumerate(qual_str):
	    q=ord(j)-33
	    if q > q_max: q_max = q
	    if q < q_min: q_min = q
	    try:
		quality[i][q] += 1
	    except:
	        quality[i][q] = 1		
    #JSH
    if not inputfile.endswith(".gz"):
       bk.close()			

    for p in range(0,read_len):
	#print str(p) + ':',
	val=[]
	occurrence=[]
        for q in range(q_min,q_max+1):
	    if quality.has_key(p) and quality[p].has_key(q):
	        val.append(str(q))				
	        occurrence.append(str(quality[p][q]))	
		q_list.append(str(quality[p][q]))
	    else:
		q_list.append(str(0))
	i_box[p] = 'rep(c(' + ','.join(val) + '),times=c(' + ','.join(occurrence) + ')/' + str(shrink)+ ')'
		
  	
    #generate R script for boxplot
	
    print >>FO, "pdf(\'%s\')" % (outfile + ".read_qual.boxplot.pdf")
    for i in sorted(i_box):
	print >>FO,'p'+str(i) + '<-' + i_box[i]
    print >>FO, 'boxplot(' + ','.join(['p'+str(i) for i in i_box]) + ',xlab=\"Position of Read(5\'->3\')\",ylab=\"Phred Quality Score\",outline=F' + ')'
    print "start"
    print >>FO,"dev.off()"
		
		
    #generate R script for heatmap
    print >>FO, '\n'
    print >>FO, "pdf(\'%s\')" % (outfile + ".read_qual.heatmap.pdf")
    print >>FO, "qual=c(" + ','.join(q_list)  + ')'
    print >>FO, "mat=matrix(qual,ncol=%s,byrow=F)" % (read_len)
    print >>FO, 'Lab.palette <- colorRampPalette(c("blue", "orange", "red3","red2","red1","red"), space = "rgb",interpolate=c(\'spline\'))'
    print >>FO, "heatmap(mat,Rowv=NA,Colv=NA,xlab=\"Position of Read\",ylab=\"Phred Quality Score\",labRow=seq(from=%s,to=%s),col = Lab.palette(256),scale=\"none\" )" % (q_min,q_max)
    print >>FO, 'dev.off()'

def main():
    #parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="reads file in fastq(fq) format. [required]")
    parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s). [required]")
    parser.add_option("-r","--reduce",action="store",type="int",dest="reduce_fold",default=1000,help="To avoid making huge vector in R, nucleotide with particular phred score represented less than this number will be ignored. Increase this number save more memory while reduce precision. This option only applies to the 'boxplot'. default=%default")
    (options,args)=parser.parse_args()

    if not (options.output_prefix and options.input_file):
	parser.print_help()
	sys.exit(0)
    inputfile=options.input_file
    outfile=options.output_prefix
    if os.path.exists(options.input_file):
	readsQual_boxplot(inputfile,outfile)
        try:
	    subprocess.call("Rscript " + options.output_prefix + ".read_qual.r",shell=True)
        except:
	    pass
    else:
	print >>sys.stderr, '\n\n' + options.input_file + " does NOT exists" + '\n'
	sys.exit(0)




if __name__ == '__main__':
	main()

