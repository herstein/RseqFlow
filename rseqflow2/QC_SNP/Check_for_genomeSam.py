#!/usr/bin/env python
'''
Created on 2012-08-19

@author: linliu
'''
import os
import sys
import optparse
import string
import re

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-s", "--sam-alnTranto", action = "store", type = "string", dest = "aln_genome_sam",
		  help = "alignments to genome (format is sam)")

(options, args) = parser.parse_args()

if (options.aln_genome_sam is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(0)

fname_genomeSam=options.aln_genome_sam

bk_genomeSam=open(fname_genomeSam)
#sh_genomeSam=bk_genomeSam.readlines()
#row_number_genomeSam=len(sh_genomeSam)
#bk_genomeSam.close()
fileline=0
############################################check genome alignments#############################################################     
for v in bk_genomeSam.xreadlines():
    fileline+=1    
    fields=v.split('\t')
    if len(fields)<12:
        print "\n File: %s , line %d" % (fname_genomeSam, fileline)
        print "  Error:SAM file should have at least 12 columns."
        sys.exit(1)
    chromosome=fields[2]

    #JSH
    # Skip entry if column 3 is '*'
    test=re.match('^\*$',chromosome)
    if test:
        continue

    if chromosome[0:3]!='chr':
	print "\n File: %s , line %d" % (fname_genomeSam, fileline)
	print "  Error:chromosome name should begin with 'chr'."
	print " please check the genomescripts name in reference sequence, line %d" % (fileline)
	print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
	sys.exit(1)
    if not chromosome[3:].isalnum():
        print "\n File: %s , line %d" % (fname_genomeSam, fileline)
        print "  Error:chromosome name should consist of digits and uppercase/lowercase characters."
        print " please check the chromosome name. Ex:chr1-cds is incorrect format. In line %d" % (fileline)
        sys.exit(1)
##############################################check genome reference#####################################################
print fname_genomeSam+" passed successfully."
