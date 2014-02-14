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
parser.add_option("-t", "--transcriptome-reference", action = "store", type = "string", dest = "transcriptome_reference_faFile",
		  help = "transcriptome reference (format is fa)")
parser.add_option("-a", "--annotation-input", action = "store", type = "string", dest = "annotation_input_gtfFile",
		  help = "transcriptome annotation (format is gtf)") 
		  

(options, args) = parser.parse_args()

if (options.transcriptome_reference_faFile is None or
    options.annotation_input_gtfFile is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(0)

fname_transcriptome=options.transcriptome_reference_faFile
fname_annotation=options.annotation_input_gtfFile

bk_transcriptome=open(fname_transcriptome)
#sh_transcriptome=bk_transcriptome.readlines()
#row_number_transcriptome=len(sh_transcriptome)
#bk_transcriptome.close()

bk_annotation=open(fname_annotation)
sh_annotation=bk_annotation.readlines()
row_number_annotation=len(sh_annotation)
bk_annotation.close()

##############################################check transcriptome reference#####################################################
chrList_tran_reference=set([])
#geneList_tran_reference=set([])
transcript_tran_reference=set([])
sequence_set=set("ACGTNacgtn")
fileline=0
temp_transcriptome=bk_transcriptome.readline()
if temp_transcriptome[0][0]!='>':
    print "\n File: %s , line %d" % (fname_transcriptome, 1)
    print "  Error:sequence name should begin with the character '>'."
    sys.exit(1)
bk_transcriptome.seek(0)
for v in bk_transcriptome.xreadlines():
    fileline+=1
    ############check sequence###############
    if v[0]!='>':
        temp_set=set(v[0:-1])
        if not temp_set<=sequence_set:
            print "\n File: %s , line %d" % (fname_transcriptome, fileline)
            print "  Error:sequence shouldn't contain any character(s) not in 'ACGTN' or 'acgtn'."
            sys.exit(1)
        continue
    ###########check sequece name############    
    temp=v[1:-1].split(' ')
    temp1=temp[0]
    temp2=temp1.split('_',2)
    if len(temp2)!=3:
        print "\n File: %s , line %d" % (fname_transcriptome, fileline)
        print "  Error:some information in transcripts title is missing."
        print " please check the transcripts name in reference sequence."
        print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
        sys.exit(1)
    if temp2[0]=='':
        print "\n File: %s , line %d" % (fname_transcriptome, fileline)
        print "  Error:genomeName is missing."
        print " please check the transcripts name in reference sequence(fa file or fa file). In line %d" % (fileline)
        print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
        sys.exit(1)
    genomeName=temp2[0]
    #geneList_tran_reference.add(genomeName)
    if temp2[1]=='':
        print "\n File: %s , line %d" % (fname_transcriptome, fileline)
        print "  Error:Annotation Source is missing."
        print " please check the transcripts name in reference sequence. In line %d" % (fileline)
        print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
        sys.exit(1)      
    temp3=temp2[2].split('=') 
    if len(temp3)!=2:
   	print "\n File: %s , line %d" % (fname_transcriptome, fileline)
	print "  Error:transcriptID or chromosome information is missing."
	print " please check the transcripts name in reference sequence. In line %d" % (fileline)
	print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
	sys.exit(1)
    if temp3[0]=='':
	print "\n File: %s , line %d" % (fname_transcriptome, fileline)
	print "  Error:transcriptID information is missing."
	print " please check the transcripts name in reference sequence. In line %d" % (fileline)
	print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
	sys.exit(1)
    transcript_id=temp3[0]
    transcript_tran_reference.add(transcript_id)
	    
    temp4=temp3[1].split(':')
    if len(temp4)!=2:
	print "\n File: %s , line %d" % (fname_transcriptome, fileline)
	print "  Error:chromosome or Start-End information is missing."
        print " please check the transcripts name in reference sequence. In line %d" % (fileline)
        print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
        sys.exit(1)
    chromosome=temp4[0]
    if chromosome[0:3]!='chr':
	print "\n File: %s , line %d" % (fname_transcriptome, fileline)
	print "  Error:chromosome name should begin with 'chr'."
	print " please check the transcripts name in reference sequence. In line %d" % (fileline)
	print "  The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
	sys.exit(1)
    if not chromosome[3:].isalnum():
        print "\n File: %s , line %d" % (fname_transcriptome, fileline)
        print "  Error:chromosome name should consist of digits and uppercase/lowercase characters."
        print " please check the chromosome name. Ex:chr1-cds is incorrect format. In line %d" % (fileline)
        sys.exit(1)
    chrList_tran_reference.add(chromosome)
    
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
##############################################check annotation reference#####################################################
chrList_tran_annotation=set([])
#geneList_tran_annotation=set([])
transcript_tran_annotation=set([])
fileline=0
#pre_gene_id=''
pre_transcript_id=''
for v in range(0, row_number_annotation):
    fileline+=1
    temp=sh_annotation[v][0:-1].split('\t')

    # JSH
    # Skip commented lines
    comment=re.match('^#',temp[0])
    if comment:
	continue

    if len(temp)!=9:
        print "\n File: %s , line %d" % (fname_annotation, fileline)
        print "  Error:GTF file should be 9-column."
        sys.exit(1)
    chromosome=temp[0]
    if chromosome[0:3]!='chr':
	print "\n File: %s , line %d" % (fname_annotation, fileline)
	print "  Error:chromosome name should begin with 'chr'."
	sys.exit(1)
    if not chromosome[3:].isalnum():
        print "\n File: %s , line %d" % (fname_annotation, fileline)
        print "  Error:chromosome name should consist of digits and uppercase/lowercase characters."
        print " please check the chromosome name. Ex:chr1-cds is incorrect format. In line %d" % (fileline)
        sys.exit(1)
    chrList_tran_annotation.add(chromosome)    
        
    type=temp[2]
    attribute=temp[8][0:-1].split('; ')
    gene_content=attribute[0]
    transcript_content=attribute[1]
    gene_temp=gene_content.split('"')
    transcript_temp=transcript_content.split('"')
    if gene_temp[0]!='gene_id ':
        print "\n File: %s , line %d" % (fname_annotation, fileline)
        print "  Error:the first attribute should be gene id and it should begin with 'gene_id '."
        print "   The correct format is gene_id \"value\"."
        sys.exit(1)
    if len(gene_temp)!=3 or gene_temp[1]=='':
        print "\n File: %s , line %d" % (fname_annotation, fileline)
        print "  Error:value of gene_id is null or the format is incorrect, please check it."
        print "   The correct format is gene_id \"value\"."
        sys.exit(1)
    if transcript_temp[0]!='transcript_id ':
        print "\n File: %s , line %d" % (fname_annotation, fileline)
        print "  Error:the second attribute should be transcript id and it should begin with 'transcript_id '."
        print "   The correct format is transcript_id \"value\"."
        sys.exit(1)
    if len(transcript_temp)!=3 or transcript_temp[1]=='':
        print "\n File: %s , line %d" % (fname_annotation, fileline)
        print "  Error:value of transcript_id is null or the format is incorrect, please check it."
        print "   The correct format is transcript_id \"value\"."
        sys.exit(1)
    gene_id=gene_temp[1]
    transcript_id=transcript_temp[1]
    #if type!='gene' and gene_id==transcript_id:
    #    print "\n File: %s , line %d" % (fname_annotation, fileline)
    #    print "  Error:the type is not gene, but gene_id and transcript_id are identical."
    #    sys.exit(1)
    if type=='gene':
        continue
    if transcript_id!=pre_transcript_id:
        pre_transcript_id=transcript_id
        if transcript_id in transcript_tran_annotation:
            print "\n File: %s , line %d" % (fname_annotation, fileline)
            print "  Error:gtf file is not sorted. Transcript:%s has been used twice." % (transcript_id)
            sys.exit(1)
        transcript_tran_annotation.add(transcript_id)
    #if gene_id!=pre_gene_id:
    #    pre_gene_id=gene_id
    #    gene_tran_annotation.add(gene_id)    

#####################################Check consistency of all input files##################################################
    
if not chrList_tran_reference <= chrList_tran_annotation:
    print " Error:chromosomes in transcriptome reference and annotation are not consistent!"
    sys.exit(1)
#if not geneList_tran_reference == geneList_tran_annotation:
#    print " Error:genes in transcriptome reference and annotation are not consistent!"
#    sys.exit(1)
if not transcript_tran_reference <= transcript_tran_annotation:
    print " Error:transcripts in transcriptome reference and annotation are not consistent!"
    sys.exit(1)
