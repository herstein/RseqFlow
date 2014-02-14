#!/usr/bin/env python
'''
Created on 2012-08-19

@author: linliu
Modified by J.Herstein 2013-04-05
'''
import os
import sys
import optparse
import string
import re

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-s", "--sam-alnTranto", action = "store", type = "string", dest = "aln_transcriptome_sam",
		  help = "alignments to transcriptome (format is sam)")

(options, args) = parser.parse_args()

if (options.aln_transcriptome_sam is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(0)

fname_transcriptome=options.aln_transcriptome_sam

bk_tranSam=open(fname_transcriptome)
#sh_tranSam=bk_tranSam.readlines()
#row_number_tranSam=len(sh_tranSam)
#bk_tranSam.close()
fileline=0
############################################check tran alignments#############################################################     
for v in bk_tranSam.xreadlines():
    fileline+=1    
    fields=v.split('\t')
    if len(fields)<12:
        print "\n File: %s , line %d" % (fname_transcriptome, fileline)
        print "  Error:SAM file should have at least 12 columns."
        sys.exit(1)
    temp1=fields[2]
   
    #JSH
    # Skip entry if column 3 is '*'
    test=re.match('^\*$',temp1)
    if test:
        continue

    temp2=temp1.split('_',2)
    if len(temp2)!=3:
        print "\n File: %s , line %d" % (fname_transcriptome, fileline)
        print "  Error:some information in transcripts title is missing."
        print " please check the transcripts name in reference sequence."
        print "  The correct format should be '$tranName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
        sys.exit(1)
    if temp2[0]=='':
        print "\n File: %s , line %d" % (fname_transcriptome, fileline)
        print "  Error:tranName is missing."
        print " please check the transcripts name in reference sequence(fa file or fa file). In line %d" % (fileline)
        print "  The correct format should be '$tranName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
        sys.exit(1)
    tranName=temp2[0]
    #geneList_tran_reference.add(tranName)
    if temp2[1]=='':
        print "\n File: %s , line %d" % (fname_transcriptome, fileline)
        print "  Error:Annotation Source is missing."
        print " please check the transcripts name in reference sequence. In line %d" % (fileline)
        print "  The correct format should be '$tranName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
        sys.exit(1)
    temp3=temp2[2].split('=') 
    if len(temp3)!=2:
   	print "\n File: %s , line %d" % (fname_transcriptome, fileline)
	print "  Error:transcriptID or chromosome information is missing."
	print " please check the transcripts name in reference sequence. In line %d" % (fileline)
	print "  The correct format should be '$tranName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
	sys.exit(1)
    if temp3[0]=='':
	print "\n File: %s , line %d" % (fname_transcriptome, fileline)
	print "  Error:transcriptID information is missing."
	print " please check the transcripts name in reference sequence. In line %d" % (fileline)
	print "  The correct format should be '$tranName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
	sys.exit(1)
    transcript_id=temp3[0]  
    temp4=temp3[1].split(':')
    if len(temp4)!=2:
	print "\n File: %s , line %d" % (fname_transcriptome, fileline)
	print "  Error:chromosome or Start-End information is missing."
        print " please check the transcripts name in reference sequence. In line %d" % (fileline)
        print "  The correct format should be '$tranName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
        sys.exit(1)
    chromosome=temp4[0]
    if chromosome[0:3]!='chr':
	print "\n File: %s , line %d" % (fname_transcriptome, fileline)
	print "  Error:chromosome name should begin with 'chr'."
	print " please check the transcripts name in reference sequence. In line %d" % (fileline)
	print "  The correct format should be '$tranName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
	sys.exit(1)
    if not chromosome[3:].isalnum():
        print "\n File: %s , line %d" % (fname_transcriptome, fileline)
        print "  Error:chromosome name should consist of digits and uppercase/lowercase characters."
        print " please check the chromosome name. Ex:chr1-cds is incorrect format. In line %d" % (fileline)
        sys.exit(1)
    
    temp5=temp4[1].split('-')
    Start=temp5[0]
    End=temp5[1]
    if not Start.isdigit():
        print "\n File: %s , line %d" % (fname_transcriptome, fileline)
        print "  Error:Start should consist of only digits."
        sys.exit(1)
    if not End.isdigit():
        print "\n File: %s , line %d" % (fname_transcriptome, fileline)
        print "  Error:Start should consist of only digits."
        sys.exit(1)
bk_tranSam.close()
##############################################check transcriptome reference#####################################################
print fname_transcriptome+" passed successfully."
