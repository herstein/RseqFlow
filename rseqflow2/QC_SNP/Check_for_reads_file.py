#!/usr/bin/env python
'''
Created on 2012-08-21

@author: linliu
Modified by J.Herstein 2013-06-04: allow fastq.gz files in addition to fastq.
'''
import os
import sys
import optparse
import string
import gzip
import subprocess

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-r", "--reads-fastq", action = "store", type = "string", dest = "reads_input_File",
		  help = "reads file (format is fastq, fq, or fastq.gz)")		  

(options, args) = parser.parse_args()

if (options.reads_input_File is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_reads=options.reads_input_File

#JSH 2013-06-04
if fname_reads.endswith(".gz"):
    bk=subprocess.Popen(["zcat", fname_reads], stdout=subprocess.PIPE,close_fds=True)
    bk_reads=bk.stdout
else:
    bk=open(fname_reads)
    bk_reads=bk.xreadlines() 
   

#bk_reads=open(fname_reads)
#sh_reads=bk_reads.readlines()
#row_number_reads=len(sh_reads)
#bk_reads.close()
####################################################check reads fq file####################################################
sequence_set=set('acgtnACGTN')
v=0
#JSH
#for fqline in bk_reads.xreadlines():
for fqline in bk_reads:
    ####-------------------------first line-------------------------####
    temp=fqline[0:-1]
    if v%4==0:
        if temp[0]!= '@':
            print "Error in file: line "+str(v+1)
            print " please check your file format(fastq file):"
            print "  Line 1 should begin with a '@' character and is followed by a sequence identifier."
            sys.exit(1)
    ####-------------------------second line-------------------------####        
    #if v+1==row_number_reads:
    #    print "Error in file: line "+str(v+1)
    #    print " Missing some information line:"
    #    print "  The information of a read shoubld be 4 lines." 
    #    sys.exit(1)
    if v%4==1:
    	sLength=len(temp)
    	temp_set=set(temp)
    	if not temp_set<=sequence_set:
            print "Error in file: line "+str(v+2)
            print " please check your file format(fastq file):"
            print "  Line 2 should be the raw sequence letters." 
            print "  And sequence shouldn't contain any characters not in 'ACGTN' or 'acgtn'."
            sys.exit(1)
    ####-------------------------third line-------------------------####
    #if v+2==row_number_reads:
    #    print "Error in file: line "+str(v+1)
    #    print " Missing some information line:"
    #    print "  The information of a read shoubld be 4 lines." 
    #    sys.exit(1)
    if v%4==2:
    	if temp[0] != '+':
            print "Error in file: line "+str(v+3)
            print " please check your file format(fastq file):"
            print "  Line 3 should begin with a '+' character."
            sys.exit(1)
####-------------------------fourth line-------------------------####
    if v%4==3:
    #f v+3==row_number_reads:
    #   print "Error in file: line "+str(v+1)
    #   print " Missing some information line:"
    #   print "  The information of a read shoubld be 4 lines." 
    #   sys.exit(1)
    	qLength=len(temp)
    #if temp[0] == '@':
    #    print "Error in file: line "+str(v+4)
    #    print " please check your file format(fastq file):"
    #    print "  Line 4 shouldn't begin with a '@' or '+' character."
    #    print "  Line 4 should encode the quality values for the sequence in Line 2."
    #    sys.exit(1)
    	if qLength!=sLength:
            print "Error in file: line "+str(v+1)
            print " please check your file format(fastq file):"
            print "  Line 4 must contain the same number of symbols as letters in sequence."
            sys.exit(1)
    v+=1
#JSH
if not fname_reads.endswith(".gz"):
    bk.close()
#bk_reads.close()
print "Reads File passed the check."
