#!/usr/bin/python

'''
Created on 2012-08-19

@author: linliu
Modified by J.Herstein 2013-04-01
'''
import os
import sys
import optparse
import string
import re

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-i", "--input-file", action = "store", type = "string", dest = "input_sam_file",
		  help = "input file to be split (sam format)")
parser.add_option("-p", "--prefix-output", action = "store", type = "string", dest = "output_prefix",
                  help = "prefix of output files")
parser.add_option("-o", "--output-chrList", action = "store", type = "string", dest = "output_chromosome_list",
		  help = "chromosome list output file")	    	  
		  

(options, args) = parser.parse_args()

if (options.input_sam_file is None or
    options.output_prefix is None or
    options.output_chromosome_list is None):
    print prog_base + ": error: missing required command-line arguments"
    parser.print_help()
    sys.exit(1)

fname_input=options.input_sam_file
fname_output=options.output_chromosome_list
prefix_output=options.output_prefix

#JSH 2013-08-26
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
      open(prefix_output+"_"+complete_chr_list[v]+"_alignment.sam", 'w').close()


bk_input=open(fname_input)
#sh_input=bk_input.readlines()
#row_number_input=len(sh_input)
#bk_input.close()

############split input file##############
#chr_fileLine={}
chr_fileLine=[]
chr_out={}
samline=0

chrlist_output=open(fname_output, 'w')

for v in bk_input.xreadlines():
    samline+=1
    temp=v[0:-1].split('\t')
#    temp1=temp[2]
#    if (temp1 is not "*"):
#        temp2=temp1.split('_',2)
#        temp3=temp2[2].split('=',1)
#        temp4=temp3[1].split(':',1)
#        chromosome=temp4[0]

    #JSH 
    temp1=re.search(r'(chr.*?):',temp[2])
    if temp1:
       chromosome=temp1.group(1)

       if (not chromosome in chr_fileLine):
	    chr_fileLine.append(chromosome)
	    chr_out[chromosome]=open(prefix_output+'_'+chromosome+'_alignment.sam', 'w')
	
       chr_out[chromosome].writelines(v)

#    try:
#        chr_fileLine[chromosome].append(v)
#    except:
#        chr_fileLine[chromosome]=[]
#        chr_fileLine[chromosome].append(v)
bk_input.close()
for chr in chr_out:
   chr_out[chr].close()
for chr in chr_fileLine:
   chrlist_output.writelines(chr+'\n')
chrlist_output.close()
##########write to files#######################
#chrlist_output=open(fname_output, 'w')
#for chr in chr_fileLine:
#    chr_out=open(prefix_output+'_'+chr+'_alignment.sam', 'w')
#    for v in range(0, len(chr_fileLine[chr])):
#        insertline=chr_fileLine[chr][v]
#        chr_out.writelines(insertline)
#    chr_out.close()
#    chrlist_output.writelines(chr+'\n')
#chrlist_output.close()    

#print 'Split Done for Input Sam Files'
print "\nSplit done for", fname_input
