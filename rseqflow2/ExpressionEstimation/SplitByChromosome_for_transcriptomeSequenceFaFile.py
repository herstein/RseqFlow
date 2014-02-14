#!/usr/bin/env python

'''
Created on 2012-08-19

@author: linliu
'''
import os
import sys
import optparse
import string

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-i", "--input-file", action = "store", type = "string", dest = "input_fa_file",
		  help = "input file to be split (format is fa)")
parser.add_option("-p", "--prefix-output", action = "store", type = "string", dest = "output_prefix",
                  help = "prefix of output")
parser.add_option("-o", "--output-chrList", action = "store", type = "string", dest = "output_chromosome_list",
		  help = "output file of chromosome list.")	    	  
		  

(options, args) = parser.parse_args()

if (options.input_fa_file is None or
    options.output_prefix is None or
    options.output_chromosome_list is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_input=options.input_fa_file
fname_output=options.output_chromosome_list
prefix_output=options.output_prefix

bk_input=open(fname_input)
sh_input=bk_input.readlines()
row_number_input=len(sh_input)
bk_input.close()

############split input file##############
chr_fileLine={}
faline=0
if sh_input[0][0]!='>':
    print "\n File: %s , line %d" % (fname_input, 1)
    print "  Error:sequence name should begin with the character '>'."
    sys.exit(1)
for v in range(0, row_number_input):
    faline+=1
    if sh_input[v][0]!='>':
        chr_fileLine[chromosome].append(sh_input[v])
        continue
    temp=sh_input[v][0:-1].split(' ')
    temp1=temp[0]
    first_underLine=temp1.find('_')
    if first_underLine==-1:
        print "\n File: %s , line %d" % (fname_input, faline)
        print "  Error:some information in transcripts title is missing."
        print " please check the transcripts name in reference sequence(fa file or fa file)."
        print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
        sys.exit(1)
    if first_underLine==0:
        print "\n File: %s , line %d" % (fname_input, faline)
        print "  Error:genomeName is missing."
        print " please check the transcripts name in reference sequence(fa file or fa file). In line %d" % (faline)
        print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
        sys.exit(1)
    temp2=temp1[first_underLine+1:]
    second_underLine=temp2.find('_')
    if second_underLine==-1:
        print "\n File: %s , line %d" % (fname_input, faline)
        print "  Error:some information in transcripts title is missing."
        print " please check the transcripts name in reference sequence(fa file or fa file)."
        print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
        sys.exit(1)
    if second_underLine==0:
        print "\n File: %s , line %d" % (fname_input, faline)
        print "  Error:Annotation Source is missing."
        print " please check the transcripts name in reference sequence(fa file or fa file). In line %d" % (faline)
        print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
        sys.exit(1)
    temp3=temp2[second_underLine+1:].split('=') 
    if len(temp3)!=2:
   	print "\n File: %s , line %d" % (fname_input, faline)
	print "  Error:transcriptID or chromosome information is missing."
	print " please check the transcripts name in reference sequence(fa file or fa file). In line %d" % (faline)
	print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
	sys.exit(1)
    temp4=temp3[1].split(':')
    if len(temp4)!=2:
	print "\n File: %s , line %d" % (fname_input, faline)
	print "  Error:chromosome or Start-End information is missing."
        print " please check the transcripts name in reference sequence(fa file or fa file). In line %d" % (faline)
        print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
        sys.exit(1)
    chromosome=temp4[0]
    if chromosome[0:3]!='chr':
	print "\n File: %s , line %d" % (fname_input, faline)
	print "  Error:chromosome name should begin with 'chr'."
	print " please check the transcripts name in reference sequence(fa file or fa file). In line %d" % (faline)
	print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
	sys.exit(1)
    if not chromosome[3:].isalnum():
        print "\n File: %s , line %d" % (fname_input, faline)
        print "  Error:chromosome name should consist of digits, uppercase and lowercase characters."
        print " please check the chromosome name. Ex:chr1-cds is incorrect format. In line %d" % (faline)
        sys.exit(1)
    temp5=temp4[1].split('-')
    Start=temp5[0]
    End=temp5[1]
    if not Start.isdigit():
        print "\n File: %s , line %d" % (fname_input, faline)
        print "  Error:Start should consist of digits only."
        sys.exit(1)
    if not End.isdigit():
        print "\n File: %s , line %d" % (fname_input, faline)
        print "  Error:Start should consist of digits only."
        sys.exit(1)
    try:
        chr_fileLine[chromosome].append(sh_input[v])
    except:
        chr_fileLine[chromosome]=[]
        chr_fileLine[chromosome].append(sh_input[v])

##########write to files#######################
chrlist_output=open(fname_output, 'w')
for chr in chr_fileLine:
    chr_out=open(prefix_output+'_'+chr+'_sequence.fa', 'w')
    for v in range(0, len(chr_fileLine[chr])):
        insertline=chr_fileLine[chr][v]
        chr_out.writelines(insertline)
    chr_out.close()
    chrlist_output.writelines(chr+'\n')
chrlist_output.close()    

print 'Split Done for Input Fasta Files'

