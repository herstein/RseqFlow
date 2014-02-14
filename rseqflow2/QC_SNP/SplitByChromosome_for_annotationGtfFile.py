#!/usr/bin/env python
'''
Created on 2012-08-19

@author: linliu
Modified by J.Herstein 2013-04-12
'''
import os
import sys
import optparse
import string
import re

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-i", "--input-file", action = "store", type = "string", dest = "input_gtf_file",
		  help = "input file to be split (format is gtf)")
parser.add_option("-p", "--prefix-output", action = "store", type = "string", dest = "output_prefix",
                  help = "prefix of output")
parser.add_option("-o", "--output-chrList", action = "store", type = "string", dest = "output_chromosome_list",
		  help = "output file of chromosome list.")	    	  
		  

(options, args) = parser.parse_args()

if (options.input_gtf_file is None or
    options.output_prefix is None or
    options.output_chromosome_list is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_input=options.input_gtf_file
fname_output=options.output_chromosome_list
prefix_output=options.output_prefix


bk_input=open(fname_input)
sh_input=bk_input.readlines()
row_number_input=len(sh_input)
bk_input.close()

### JSH 2013-08-26
############ create initial empty files (required for PEGASUS workflow) #######
complete_chr_list=[]
# range endpoint is never included so use 23 to get 22 chromosomes. 
for v in range(1,23):
      complete_chr_list.append("chr"+str(v))
# don't know if gtf will have chrM or chrMT, so including both
complete_chr_list.append("chrX")
complete_chr_list.append("chrY")
complete_chr_list.append("chrM")
complete_chr_list.append("chrMT")
for v in range(0,len(complete_chr_list)):
      open(prefix_output+"_"+complete_chr_list[v]+"_annotation.gtf", 'w').close()

############split input file##############
chr_fileLine={}
chr_transcriptID={}
gtfline=0
pre_transcript_id=''
for v in range(0, row_number_input):
    gtfline+=1
    temp=sh_input[v][0:-1].split('\t')
    chromosome=temp[0]
    #JSH
    # skip comment lines
    #    if ('#' in temp[0]):
    if re.match('^#',chromosome):
        continue

    if chromosome[0:3]!='chr':
	      print "\n File: %s , line %d" % (fname_input, gtfline)
	      print "  Error:chromosome name should begin with 'chr'."
	      sys.exit(1)
    if not chromosome[3:].isalnum():
        print "\n File: %s , line %d" % (fname_input, gtfline)
        print "  Error:chromosome name should consist of digits and uppercase/lowercase characters."
        print " please check the chromosome name. Ex:chr1-cds is incorrect format. In line %d" % (gtfline)
        sys.exit(1)
    type=temp[2]
    attribute=temp[8][0:-1].split('; ')
    gene_content=attribute[0]
    transcript_content=attribute[1]
    gene_temp=gene_content.split('"')
    transcript_temp=transcript_content.split('"')
    if gene_temp[0]!='gene_id ':
        print "\n File: %s , line %d" % (fname_input, gtfline)
        print "  Error:the first attribute should be gene id and it should begin with 'gene_id '."
        print "   The correct format is gene_id \"value\"."
        sys.exit(1)
    if len(gene_temp)!=3 or gene_temp[1]=='':
        print "\n File: %s , line %d" % (fname_input, gtfline)
        print "  Error:value of gene_id is null or the format is incorrect, please check it."
        print "   The correct format is gene_id \"value\"."
        sys.exit(1)
    if transcript_temp[0]!='transcript_id ':
        print "\n File: %s , line %d" % (fname_input, gtfline)
        print "  Error:the second attribute should be transcript id and it should begin with 'transcript_id '."
        print "   The correct format is transcript_id \"value\"."
        sys.exit(1)
    if len(transcript_temp)!=3 or transcript_temp[1]=='':
        print "\n File: %s , line %d" % (fname_input, gtfline)
        print "  Error:value of transcript_id is null or the format is incorrect, please check it."
        print "   The correct format is transcript_id \"value\"."
        sys.exit(1)
    gene_id=gene_temp[1]
    transcript_id=transcript_temp[1]
    #if type!='gene' and gene_id==transcript_id:
    #    print "\n File: %s , line %d" % (fname_input, gtfline)
    #    print "  Error:the type is not gene, but gene_id and transcript_id are identical."
    #    sys.exit(1)
    if type=='gene':
        continue
    if transcript_id!=pre_transcript_id:
        pre_transcript_id=transcript_id
        if chromosome in chr_transcriptID and transcript_id in chr_transcriptID[chromosome]:
            print "\n File: %s , line %d" % (fname_input, gtfline)
            print "  Error:gtf file is not sorted. Transcript:%s has been used twice." % (transcript_id)
            sys.exit(1)
        try:
            chr_transcriptID[chromosome].append(transcript_id)    
        except:
            chr_transcriptID[chromosome]=[]
            chr_transcriptID[chromosome].append(transcript_id)
    try:
        chr_fileLine[chromosome].append(sh_input[v])
    except:
        chr_fileLine[chromosome]=[]
        chr_fileLine[chromosome].append(sh_input[v])

##########write to files#######################
chrlist_output=open(fname_output, 'w')
for chr in chr_fileLine:
    chr_out=open(options.output_prefix+'_'+chr+'_annotation.gtf', 'w')
    for v in range(0, len(chr_fileLine[chr])):
        insertline=chr_fileLine[chr][v]
        chr_out.writelines(insertline)
    chr_out.close()
    chrlist_output.writelines(chr+'\n')
chrlist_output.close()    

print 'Split Done for Input GTF Files'

