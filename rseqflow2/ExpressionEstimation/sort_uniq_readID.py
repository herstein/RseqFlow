#!/usr/bin/env python

'''
Created on 2009-10-30
Modified on 2012-07-18

@author: Administrator
@modify by Liu Lin
'''

import os
import sys
import optparse
import copy
# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-s", "--sam-alignment", action = "store", type = "string", dest = "alignment_sam_file",
		  help = "sam file of alignment")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output",
		  help = "output: readID list after sorting")

(options, args) = parser.parse_args()

if (options.alignment_sam_file is None or
    options.output is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname = options.alignment_sam_file
outputFilename = options.output

#fname=sys.argv[1]
#outputFilename=sys.argv[2]

try:
    bk=open(fname)
except:
    print prog_base + ": error: cannot open file " + fname
    sys.exit(1)
#sh=bk.readlines()
#row_number=len(sh)
#bk.close()

readID_repeat=[]
readID_list=[]

##############################################################################
#-------------------get reads information from sam file------------------------------#
for v in bk.xreadlines():
    temp=v[0:-1].split('\t')
    if temp[0][0] != '@':
    		readID=temp[0]
    		readID_repeat.append(readID)
bk.close()
readID_list=list(set(sorted(readID_repeat)))    		
##############################################################################            
#############################################################################
######write to file

try:
    OUT=open(outputFilename,'w')
except:
    print prog_base + ": error: cannot open file " + outputFilename
    sys.exit(1)
for v in range(0,len(readID_list)):
		print >>OUT, "%s" % (str(readID_list[v]))    
OUT.close()
