#!/usr/bin/env python

# this program will count the bases mapped to exons, and then count RPKM 
# to get the approximate exon's expression level 

import os
import sys
import optparse

# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-s", "--uniq-map", action = "store", type = "string", dest = "uniq_map",
		  help = "uniq map sam file")
parser.add_option("-m", "--mappable-exon-index", action = "store", type = "string", dest = "mappable_exon_index",
		  help = "mappable index no repeat regions")
parser.add_option("-e", "--exon-index", action = "store", type = "string", dest = "exon_index",
		  help = "original exon index")
parser.add_option("-l", "--exon-length", action = "store", type = "string", dest = "exon_length",
		  help = "exon length no repeat region")
parser.add_option("-c", "--mapping-count", action = "store", type = "string", dest = "mapping_count",
		  help = "mapping count")
parser.add_option("-o", "--output1", action = "store", type = "string", dest = "output1",
		  help = "output: exon expression level file")
#parser.add_option("-n", "--output2", action = "store", type = "string", dest = "output2",
                  #help = "output: read number in each exon")

(options, args) = parser.parse_args()

if (options.uniq_map is None or
    options.exon_index is None or
    options.mappable_exon_index is None or
    options.exon_length is None or
    options.mapping_count is None or
    options.output1 is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

fname_mapping_result = options.uniq_map
fname_exon = options.mappable_exon_index
fname_old_exon = options.exon_index
fname_exon_length = options.exon_length
fname_mapping_count = options.mapping_count
output_filename = options.output1


try:
    bk_exon=open(fname_exon)
except:
    print prog_base + ": error: cannot open file " + fname_exon
    sys.exit(1)
sh_exon=bk_exon.readlines()
row_number_exon=len(sh_exon)
bk_exon.close()

exon_express={}
exon_information={}
isoform_gene_pair={}
gene_strand={}
exon_read={}

for w in range(1,len(sh_exon)):###build dictory for exons for isoform
    temp1=sh_exon[w][0:-1].split('\t')
    if int(temp1[3])!=0:
        gene_name=temp1[1]
        strand=temp1[2]
        gene_strand[gene_name]=strand
        temp4=str(temp1[4])[2:-2].split('], [')
        temp5=[]
        exon_info=[]
        for vv in range(0,len(temp4)):
            temp6=temp4[vv].split(',')
            exon_express[str(gene_name)+':'+str(int(temp6[1]))]=0
            exon_read[str(gene_name)+':'+str(int(temp6[1]))]=0
            temp5=temp5+temp6
        exon_info=exon_info+temp5
        exon_information[temp1[0]]=exon_info
        isoform_gene_pair[temp1[0]]=gene_name

try:
    bk_old_exon=open(fname_old_exon)
except:
    print prog_base + ": error: cannot open file " + fname_old_exon
    sys.exit(1)
sh_old_exon=bk_old_exon.readlines()
bk_old_exon.close()

serial_location_pair={} #the table of gene:serial gene:locationstart:locationend

for w in range(1,len(sh_old_exon)):###build dictory for exons for isoform
    temp1=sh_old_exon[w][0:-1].split('\t')
    gene_name=temp1[1]
    temp4=str(temp1[4])[2:-3].split('), (')
    temp5=[]
    for vv in range(0,len(temp4)):
        temp6=temp4[vv].split(',')
        serial_location_pair[str(gene_name)+':'+str(int(temp6[1]))]=gene_name+':'+str(int(temp6[2]))+':'+str(int(temp6[3]))
    
try:
    bk_mapping_result=open(fname_mapping_result)
except:
    print prog_base + ": error: cannot open file " + fname_mapping_result
    sys.exit(1)
sh_mapping_result=bk_mapping_result.readlines()
row_number_mapping_result=len(sh_mapping_result)
bk_mapping_result.close()
if row_number_mapping_result==0:
    print "No reads mapped to only one genen for current chromosome, so no expression level can be caculated base on repeat region library."
    sys.exit(0)

try:
    bk_exon_length=open(fname_exon_length)
except:
    print prog_base + ": error: cannot open file " + fname_exon_length
    sys.exit(1)
sh_exon_length=bk_exon_length.readlines()
bk_exon_length.close()

exon_length={}

for i in range(1,len(sh_exon_length)):
    temp1=sh_exon_length[i][0:-1].split('\t')
    temp2=temp1[0].split(':')
    exon=temp2[0]+':'+temp2[1]
    exon_length[exon]=int(temp1[1])


try:
    bk_mapping_count=open(fname_mapping_count)
except:
    print prog_base + ": error: cannot open file " + fname_mapping_count
    sys.exit(1)
sh_mapping_count=bk_mapping_count.readlines()
row_number_mapping_count=len(sh_mapping_count)
bk_mapping_count.close()


###look for read_length and #mappedreads 
 
temp_position=sh_mapping_count[1].find(':')
temp3=sh_mapping_count[0].find(':')
if temp3!=-1 and temp_position!=-1:
   mapped_reads=int(sh_mapping_count[0][temp3+1:-1])
   read_length=int(sh_mapping_count[1][temp_position+1:-1])


#############################################
isoform_name=[]
chromosome = ''

##############################################################################
####### the mapping result
for v in range(0,row_number_mapping_result): #
     
    mapping_info=sh_mapping_result[v][0:-1].split('\t')
    temp2=mapping_info[2].split('_',2)
    temp3=temp2[2].split(':',1)
    temp4=temp3[0].split('=',1)
    chromosome=temp4[1]
    mapping_info[2]=temp4[0]                #####take out the isoform name
    gene_name=isoform_gene_pair.get(mapping_info[2],'none')
    isoform_name.append(mapping_info[2])
    if gene_name=='none':
        continue

    exon_info=[]
    exon_info=exon_information[isoform_name[v]]
    for p in range(0,(len(exon_info)/4)):
        if p>=0: 
            if int(mapping_info[3])>=int(exon_info[4*p]) and int(mapping_info[3])+read_length-1<=(int(exon_info[4*p])+int(exon_info[4*p+3])-int(exon_info[4*p+2])):
                temp22=gene_name+':'+str(int(exon_info[4*p+1]))
                temp11=exon_express[temp22]
                exon_express[temp22]=temp11+read_length

            elif int(mapping_info[3])>=int(exon_info[4*p]) and int(mapping_info[3])<=(int(exon_info[4*p])+int(exon_info[4*p+3])-int(exon_info[4*p+2])) and int(mapping_info[3])+read_length-1>=(int(exon_info[4*p])+int(exon_info[4*p+3])-int(exon_info[4*p+2])):
                temp22=gene_name+':'+str(int(exon_info[4*p+1]))
                temp11=exon_express[temp22]
                exon_express[temp22]=temp11+(int(exon_info[4*p])+int(exon_info[4*p+3])-int(exon_info[4*p+2]))-int(mapping_info[3])+1
            elif int(mapping_info[3])<=int(exon_info[4*p]) and int(mapping_info[3])+read_length-1>=int(exon_info[4*p]) and int(mapping_info[3])+read_length-1<=(int(exon_info[4*p])+int(exon_info[4*p+3])-int(exon_info[4*p+2])):
                temp22=gene_name+':'+str(int(exon_info[4*p+1]))
                temp11=exon_express[temp22]
                exon_express[temp22]=temp11+int(mapping_info[3])+read_length-int(exon_info[4*p])

            elif int(mapping_info[3])<=int(exon_info[4*p]) and int(mapping_info[3])+read_length-1>=(int(exon_info[4*p])+int(exon_info[4*p+3])-int(exon_info[4*p+2])):
                temp22=gene_name+':'+str(int(exon_info[4*p+1]))
                temp11=exon_express[temp22]
                exon_express[temp22]=temp11+int(exon_info[4*p+3])-int(exon_info[4*p+2])+1

for v in exon_express.keys():
    length=exon_length.get(v,0)

    if exon_express[v]!=0 and length!=0:
        temp=v.split(':')
	geneName=''
        for s in (0, len(temp)-2):
            geneName+=temp[s]
        rpkm=(float(exon_express[v])/read_length/length/mapped_reads)*1000000000
        read=float(exon_express[v])/read_length
        exon_express[v]=rpkm
        exon_read[v]=read



try:
    file_output=open(output_filename, 'w')
except:
    print prog_base + ": error: cannot open file " + output_filename
    sys.exit(1)
method_manual='# There are 4 methods to deal with multi-mapped reads:'+'\n'
method_manual+='# 1.Unique: remove the multi-mapped reads directly.'+'\n'
method_manual+='# 2.Random: assign the multi-mapped reads randomly.'+'\n'
method_manual+='# 3.Proportion: assign the multi-mapped reads according to the proportion of gene expression level.'+'\n'
method_manual+='# 4.Mappable Library: only count the reads falling into mappable regions, and it is based on the unmappable repetitive region library, which can be built automatically.'+'\n'
method_manual+='# You chose the 4th one. This method only counts the reads falling into mappable regions, and it is based on the unmappable repetitive region library, which can be built automatically.'+'\n'
#file_output.writelines(method_manual)
titles='#GeneID:ExonStartLocation:ExonEndLocation'+'\t'+'Chromosome'+'\t'+'Strand'+'\t'+'RPKM'+'\n'
file_output.writelines(titles)

for v in exon_express.keys():
    temp=serial_location_pair[v]
    temp1=temp.split(':')
    geneName=":".join(temp1[0:-2])
    insert_line=str(temp)+'\t'+chromosome+'\t'+gene_strand[geneName]+'\t'+format(exon_read[v],'.2f')+'\t'+str(exon_express[v])+'\n'
    file_output.writelines(insert_line)
file_output.close()

print 'ExonExpressionLevel Done for Chromosome %s' % (chromosome)
