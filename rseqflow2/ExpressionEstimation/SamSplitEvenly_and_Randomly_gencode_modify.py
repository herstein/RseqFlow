#!/usr/bin/python

#JSH 2013-10-31
#if a read mapped to multiple genes, only keep the reads mapped to single gene
# Bowtie2 default is to keep only the best valid alignment
# It uses the SAM optional field XS: to report the alignment score for the 
#      2nd best alignment for a read.
# If "XS:i:" field is present, assume the read was multimapped, otherwise assume it's unique

#ORIG CODE:
#if a read mapped to multiple genes, only keep the reads mapped to single gene
#and then split the multi-mapped reads into the several genes evenly
#or assign the multi-mapped reads to one of the genes randomly

import os
import sys
import random
import optparse
import re

#fname_sam='c:/Python26/sample.txt'
#fname_table='c:/Python26/isoform_gene_table_ref19.txt'
#output='MappingSingleGene.txt'
#output2='SplitbyGenelength.txt'
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-s", "--sam-mappingResult", action = "store", type = "string", dest = "mappingResult_sam_file",
		  help = "mapping result in sam format")
#JSH 2013-10-31 Leaving this in for now but it only applies to orig code
parser.add_option("-g", "--gtf-annotation", action = "store", type = "string", dest = "anotation_gtf_file",
		  help = "annotation in gtf format")		  
parser.add_option("-u", "--output-unique",  action = "store", type = "string", dest = "uniqueMap_sam_output",
		  help = "uniquely mapped reads result in sam format")		  
parser.add_option("-m", "--output-multiple", action = "store", type = "string", dest = "multiMap_txt_output",
		  help = "multiple mapped reads result in txt file")		  

(options, args) = parser.parse_args()

if (options.mappingResult_sam_file is None or
# JSH 2013-10-31 Commented out, only required for orig code
#    options.anotation_gtf_file is None or
    options.uniqueMap_sam_output is None or
    options.multiMap_txt_output is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_sam=options.mappingResult_sam_file # the original mapping file
#JSH 2013-10-31 Commented out, part of orig code
#fname_gtf=options.anotation_gtf_file # isoform_gene_gtf
fname_uniq=options.uniqueMap_sam_output # file for reads that map to single gene
fname_multi=options.multiMap_txt_output # file for reads that map to multiple genes

#JSH 2013-10-31 Commented out, part of orig code
#bk_gtf=open(fname_gtf)
#sh_gtf=bk_gtf.readlines()
#bk_gtf.close()

#JSH 2013-10-31 Commented out, part of orig code
##JSH 2013-08-26
#### create empty file if chr not present (required for PEGASUS workflow) ###
#if(sh_gtf==0):
#    print "SamSplitEvenly_and_Randomly_gencode_modify.py: No chr data for " + fname_gtf
#    try:
#        open(fname_uniq,'w').close()
#    except:
#        print prog_base + ": error: cannot create file " + fname_uniq
#        sys.exit(1)
#    try:
#        open(fname_multi,'w').close()
#    except:
#        print prog_base + ": error: cannot create file " + fname_multi
#        sys.exit(1)
#
#    sys.exit(0)
#
#chromosome=None

#isoform_gene={}
#gene_readsplit={}
#gene_readrandom={}
###########generate dictionary: key-isoform, value-gene########
#pre_transcript_id	='none'
#samline=0
#for v in range(0,len(sh_gtf)):
#  temp=sh_gtf[v][0:-1].split('\t')
#  temp1		=temp[8]
#  attributes	=temp1.split(';')
#  temp2		=attributes[0].split('"')
#  gene_id	=temp2[1]
#  temp3		=attributes[1].split('"')
#  transcript_id	=temp3[1]
#  #if gene_id==transcript_id:
#  #  print "Gene ID:"+gene_id+'\n'
#  #  print "Transcrit ID:"+transcript_id+'\n'
#  #  print "Gene ID and Transcript ID are the same, so please check the annot." 
#    #sys.exit(1)
#  if transcript_id!=pre_transcript_id:
#    isoform_gene[transcript_id]=gene_id
#    pre_gene_id=gene_id
#    pre_transcript_id=transcript_id

bk_sam=open(fname_sam)
#sh_sam=bk_sam.readlines()
#bk_sam.close()

#JSH 2013-10-31 Commented out, part of orig code
#read_alignments={}
#read_genes={}
#samline=0

#JSH 2013-10-31 Added the following 3 lines:
uniq_output=open(fname_uniq, 'w')
multi_output=open(fname_multi, 'w')
secondary=256 #secondary read tag from bowtie2 sam file

for v in bk_sam.xreadlines():
#  samline+=1
  if v[0]=='@':
    # sam header
    continue
  temp1=v[0:-1].split('\t')
#  readID=temp1[0]
  # Check that reference field is not "*"
  if re.match('^\*',temp1[2]):
      continue

  #JSH 2013-10-31 Added:
  # If field "XS:i:" is present, read is multimapped
  if re.search(r'XS:i:',v):
      multi_output.writelines(v)
  else:
      fields=v[0:-1].split('\t')
      flag=int(fields[1])
      if (flag&secondary):
          multi_output.writelines(v)
      else:
          uniq_output.writelines(v)

bk_sam.close()
uniq_output.close()
multi_output.close()

print "Done splitting %s for unique and multiple mapped reads" % fname_sam

# JSH 2013-10-31 THIS IS THE NEW END OF FILE 

#JSH 2013-10-31 Commented out the rest of the code, part of orig code
#  temp2=temp1[2].split('_',2)
#  #if len(temp2)!=3:
#  #  print "File: %s , line %d" % (fname_sam, samline)
#  #  print "Error:some information in transcripts title is missing."
#  #  print "please check the transcripts name in reference sequence(fa file or sam file)."
#  #  print "The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
#  #  sys.exit(1)
#  #if (temp2[0]=='') or (temp2[1]==''):
#  #  print "File: %s , line %d" % (fname_sam, samline)
#  #  print "Error:genomeName or AnnotationSource is missing."
#  #  print "please check the transcripts name in reference sequence(fa file or sam file). In line %d" % (samline)
#  #  print "The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
#  #  sys.exit(1)
#  temp3=temp2[2].split('=')
#  #if len(temp3)!=2:
#  #  print "File: %s , line %d" % (fname_sam, samline)
#  #  print "Error:transcriptID or chromosome information is missing."
#  #  print "please check the transcripts name in reference sequence(fa file or sam file). In line %d" % (samline)
#  #  print "The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
#  #  sys.exit(1)
#  #if temp3[0]=='':
#  #  print "File: %s , line %d" % (fname_sam, samline)
#  #  print "Error:transcriptID information is missing."
#  #  print "please check the transcripts name in reference sequence(fa file or sam file). In line %d" % (samline)
#  #  print "The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
#  #  sys.exit(1)
#  #if temp3[1]=='':
#  #  print "File: %s , line %d" % (fname_sam, samline)
#  #  print "Error:chromosome information is missing."
#  #  print "please check the transcripts name in reference sequence(fa file or sam file). In line %d" % (samline)
#  #  print "The correct format should be '$genomeName_$AnnotationSource_$TranscriptsID=$Chromosome:$Start-$End'."
#  #  sys.exit(1)
#  isoform=temp3[0]
#  temp4=temp3[1].split(':')
#  chromosome=temp4[0]
#  #JSH -- commented read_length since it is now handled in Get_ReadsMappingInformation.py
#  #temp6=temp1[5].split('M')
#  #read_length=int(temp6[0])
#  try:
#    gene=isoform_gene[isoform]
#  except:
#    #print "Warning:transcript ID in annotation and transcript ID in sam file are not consistent. File:%s, in line %d" % (fname_sam, samline)
#    #print "        If all transcripts in two files are not consistent, the result will be null."		
#    continue
#  try:
#    read_alignments[readID].append(v)
#    read_genes[readID].append(gene)
#  except:
#    read_alignments[readID]=[]
#    read_genes[readID]=[]
#    read_alignments[readID].append(v)
#    read_genes[readID].append(gene)
#bk_sam.close()
#############get read length and chromosome##########
##temp=sh_sam[-1].split('\t')
##temp1=temp[5].split('M')
##read_length=int(temp1[0])
##temp2=temp[2].split('=')
##temp3=temp2[1].split(':')
##chromosome=temp3[0]
#		
#file_output=open(fname_uniq, 'w')
#file_output2=open(fname_multi, 'w')
#
#read_count=0
#for readID in read_alignments:
#	read_count+=1
#	gene=[]
#	same_read_mapping=[]
#	gene=read_genes[readID]
#	same_read_mapping=read_alignments[readID]
#  	genelist=list(set(gene))
# 	if len(genelist)>1:
#    		gene_transcriptPosition={}
#    		for v in range(0, len(same_read_mapping)):
#      			temp=same_read_mapping[v].split('\t')
#      			temp1=temp[2].split('_',2)
#      			temp2=temp1[2].split('=')
#      			transcriptID=temp2[0]
#      			mapPosition=temp[3]
#      			transcript_position=str(transcriptID)+':'+str(mapPosition)
#      			geneID=isoform_gene[transcriptID]
#      			try:
#        			gene_transcriptPositions[geneID].append(transcript_position)
#      			except:
#        			gene_transcriptPosition[geneID]=[]
#        			gene_transcriptPosition[geneID].append(transcript_position)
#      		writeline=str(readID)+'\t'+chromosome+'\t'
#      		for geneID in gene_transcriptPosition:
#        		transcriptPosition_list=gene_transcriptPosition[geneID]
#        		length=len(transcriptPosition_list)
#       			randomly=random.randint(0,length-1) #create the random number to assign the read randomly
#       			writeline += str(transcriptPosition_list[randomly])
#			writeline += ','
#		file_output2.writelines(writeline+'\n')
#                
#  	else:
#  		length=len(same_read_mapping)
#  		randomly=random.randint(0,length-1)
#      		file_output.writelines(str(same_read_mapping[randomly]))
###########################################        
#        
#file_output.close()
#
#print "Uniquely mapped result and multiple mapped result Done for Current Chromosome"
#if chromosome is not None:
#    print "Uniquely mapped result and multiple mapped result Done for Chromosome: %s" % (chromosome)
