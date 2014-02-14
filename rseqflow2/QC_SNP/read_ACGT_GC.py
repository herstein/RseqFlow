#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
Calculate distribution of reads' GC content
modify by liulin
modify date:2012-9-11
modified by J.Herstein 2013-06-04 to allow fastq.gz format
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
import subprocess
import decimal
import gzip
import subprocess


#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

#import my own modules
#from qcmodule import SAM
#changes to the paths

#changing history to this module

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Reads file in fastq(fq) format.")
parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s).")

(options,args)=parser.parse_args()

if (options.input_file is None or
    options.output_prefix is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)
    
inputfile=options.input_file
outputPrefix=options.output_prefix


'''GC content distribution of reads'''
outfile1 = outputPrefix + ".GC_hist.xls"
outfile2 = outputPrefix + ".GC_hist_plot.r"
FO=open(outfile1,'w')
RS1=open(outfile2,'w')

Output_ACGT_file = outputPrefix + ".ACGT_Percentage.txt"
Output_GC_content_file = outputPrefix + ".GC_Percentage.txt"
Rfile = outputPrefix + ".ATCG_GC_Percentage_plot.r"
		
############open input file###################
if inputfile.endswith(".gz"):
    print "Load reads file(fastq.gz format)..."
else:
    print "Load reads file(fastq format)..."

try:
    #JSH 2013-06-04
    if inputfile.endswith(".gz"):
          bk_r=subprocess.Popen(["zcat", inputfile], stdout=subprocess.PIPE,close_fds=True)
          bk=bk_r.stdout
          temp=bk.readline()
          temp=bk.readline()
          read_length=len(temp[0:-1])
          # Clunky, but I'm not sure how to "rewind" to beginning 
          bk_r=subprocess.Popen(["zcat", inputfile], stdout=subprocess.PIPE,close_fds=True)
          bk_reads=bk_r.stdout
    else:
          bk=open(inputfile)
          temp=bk.readline()
          temp=bk.readline()
          read_length=len(temp[0:-1])
          bk.seek(0)
          bk_reads=bk.xreadlines()

except:
    print prog_base + ": error: cannot open file " + inputfile
    sys.exit(1)
#sh=bk.readlines()
#row_number=len(sh)
#bk.close()
###########read sequence into a list##########
gc_hist=collections.defaultdict(int)	#key is GC percent, value is count of reads
#temp=bk.readline()
#temp=bk.readline()
#read_length=len(temp[0:-1])
#bk.seek(0)
initial_value=[0,0,0,0,0]
A_num=[]
C_num=[]
G_num=[]
T_num=[]
N_num=[]
for v in range(0, read_length):
    A_num.append(0)
    C_num.append(0)
    G_num.append(0)
    T_num.append(0)
    N_num.append(0)
total_reads=0
count=0
#JSH
#for v in bk.xreadlines():
for v in bk_reads:
    count+=1
    if (count%4)!=2:
	continue
    total_reads+=1
    temp=v[0:-1]
    aligned_read = temp
    RNA_read = aligned_read.upper()
    gc_percent = "%4.2f" % ((RNA_read.count('C') + RNA_read.count('G'))/(read_length+0.0)*100)
    gc_hist[gc_percent] += 1
    for vv in range(0,read_length):
        if RNA_read[vv]=='A':
	    A_num[vv]+=1
	elif RNA_read[vv]=='C':
            C_num[vv]+=1
        elif RNA_read[vv]=='G':
            G_num[vv]+=1 
        elif RNA_read[vv]=='T':
            T_num[vv]+=1
	elif RNA_read[vv]=='N':
            N_num[vv]+=1	 
#JSH 
if not inputfile.endswith(".gz"):
   bk.close()
####caculating GC content###		                                       	
print >>sys.stdout, "writing GC content ..."	
print >>FO, "GC%\tread_count"
gc_sort_keys=[]
for v in gc_hist:
    gc_sort_keys.append(decimal.Decimal(v))
gc_sort_keys.sort()
for i in gc_sort_keys:
    i_str=str(i)+'%'
    print >>FO, i_str + '\t' + format(gc_hist[str(i)],',')			
print >>sys.stdout, "writing R script ..."
print >>RS1, "pdf(\"%s\")" % (outputPrefix +  ".GC_hist_plot.pdf")
print >>RS1, 'gc=rep(c(' + ','.join([i for i in gc_hist.keys()]) + '),' + 'times=c(' + ','.join([str(i) for i in gc_hist.values()]) + '))'
print >>RS1, 'hist(gc,probability=T,breaks=%d,xlab="GC content (%%)",ylab="Density of Reads",border="blue",main="")' % 100
print >>RS1,"dev.off()"
RS1.close()
###caculating A,C,G,T,N percentage#########
read_position=[]
totalCount_str=[]
A_str=[]
C_str=[]
G_str=[]
T_str=[]
N_str=[]
OUT_ACGT=open(Output_ACGT_file,'w')
OUT_GC=open(Output_GC_content_file,'w')
print >>OUT_ACGT, "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s" % ('Position In Read','#Reads','A%','C%','G%','T%','N%')
print >>OUT_GC, "%10s\t%10s\t%10s" % ('Position In Read','#Reads','GC content')
for v in range(0,read_length):
    position=v+1
    ap=A_num[v]*100/float(total_reads)
    cp=C_num[v]*100/float(total_reads)
    gp=G_num[v]*100/float(total_reads)
    tp=T_num[v]*100/float(total_reads)
    np=N_num[v]*100/float(total_reads)
    print >>OUT_ACGT, "%16d\t%10s\t%9.2f%%\t%9.2f%%\t%9.2f%%\t%9.2f%%\t%9.2f%%" % (int(position),format(total_reads,','),ap,cp,gp,tp,np)
    gc_p=gp+cp
    print >>OUT_GC, "%16d\t%10s\t%9.2f%%" % (int(position),format(total_reads,','),gc_p)
    A_str.append(str(A_num[v]))
    C_str.append(str(C_num[v]))
    G_str.append(str(G_num[v]))
    T_str.append(str(T_num[v]))
    N_str.append(str(N_num[v]))
    read_position.append(str(position))
    totalCount_str.append(str(total_reads))
OUT_ACGT.close()
OUT_GC.close()

RS2=open(Rfile,'w')
print >>sys.stdout, "generating R script  ..."
print >>RS2, "position=c(" + ','.join(read_position) + ')'
print >>RS2, "A_count=c(" + ','.join(A_str) + ')'
print >>RS2, "C_count=c(" + ','.join(C_str) + ')'
print >>RS2, "G_count=c(" + ','.join(G_str) + ')'
print >>RS2, "T_count=c(" + ','.join(T_str) + ')'
print >>RS2, "total=c("+','.join(totalCount_str)+')'

print >>RS2, "ym=max(A_count/total,C_count/total,G_count/total,T_count/total) + 0.05"
print >>RS2, "yn=min(A_count/total,C_count/total,G_count/total,T_count/total)"
print >>RS2, "xm=max(position) + 2"
print >>RS2, "xn=0" 
		
print >>RS2, 'pdf(\"%s\")' % (outputPrefix +".ACGT_Percentage_plot.pdf")
print >>RS2, 'plot(position,A_count/total,type="o",pch=20,xlim=c(xn,xm),ylim=c(yn,ym),col="green",xlab="Position in read(bp)",ylab="Nucleotide Frequency",main="Sequence content across all bases")'
print >>RS2, 'lines(position,T_count/total,type="o",pch=20,col="red")'
print >>RS2, 'lines(position,G_count/total,type="o",pch=20,col="black")'
print >>RS2, 'lines(position,C_count/total,type="o",pch=20,col="blue")'
print >>RS2, 'legend('+ str(len(read_position)-12) + ',ym,legend=c("%A","%T","%G","%C"),col=c("green","red","black","blue"),lwd=2,pch=20,text.col=c("green","red","black","blue"))'

print >>RS2, "ym=max((C_count+G_count)/total) + 0.05"
print >>RS2, "yn=min((C_count+G_count)/total) - 0.05"
print >>RS2, 'pdf(\"%s\")' % (outputPrefix +".GC_Percentage_plot.pdf")
print >>RS2, 'plot(position,(C_count+G_count)/total,type="o",pch=20,xlim=c(xn,xm),ylim=c(yn,ym),col="red",xlab="Position in read(bp)",ylab="Nucleotide Frequency",main="GC content across all bases")'
print >>RS2, 'legend('+ str(len(read_position)-12) + ',ym,legend=c("%GC"),col=c("red"),lwd=2,pch=20,text.col=c("red"))'
		
RS2.close()

try:
    subprocess.call("Rscript " + outfile2, shell=True)
    subprocess.call("Rscript " + Rfile, shell=True)
except:
    pass
print "Caculating ATCG percentage done."
print "GC content across the read length and hist of GC content done."
