#!/usr/bin/env python
'''
Created on 2009-10-30

@author: Administrator
'''
import sys
import copy
import optparse
import os

# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("--c1", action = "store", type = "string", dest = "condition1",
                  help = "specify the name of condition1")
parser.add_option("--c2", action = "store", type = "string", dest = "condition2",
                  help = "specify the name of condition2")
parser.add_option("-1", action = "store", type = "string", dest = "con1_exoncount",
                  help = "Exon Read Count from condition1")
parser.add_option("-2", action = "store", type = "string", dest = "con2_exoncount",
                  help = "Exon Read Count from condition2")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output_c1c2",
                  help = "Combine output of condition1 and condition2")


(options, args) = parser.parse_args()

if (options.condition1 is None or
    options.condition2 is None or
    options.con1_exoncount is None or
    options.con2_exoncount is None or
    options.output_c1c2 is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

con1 = options.condition1
con2 = options.condition2
fname1 = options.con1_exoncount
fname2 = options.con2_exoncount
outputFilename = options.output_c1c2



bk1=open(fname1)
sh1=bk1.readlines()
bk1.close()

bk2=open(fname2)
sh2=bk2.readlines()
bk2.close()


gene1={}
gene2={}

for v in range(1, len(sh1)):
    temp=sh1[v][0:-1].split('\t')
    gene1[temp[0]]=round(float(temp[3]))
  
for v in range(1, len(sh2)):
    temp=sh2[v][0:-1].split('\t')
    temp1=gene1.get(temp[0],'None')
    if temp1!='None':
        gene2[temp[0]]=round(float(temp[3]))

ws=open(outputFilename,'w')
title='gene'+'\t'+'exon'+'\t'+con1+'\t'+con2+'\n'
ws.writelines(title)

for v in gene2.keys():
    temp=v.split(':')
    ws.writelines(temp[0]+'\t'+temp[1]+':'+temp[2]+'\t'+str(int(gene1[v]))+'\t'+str(int(gene2[v]))+'\n')
ws.close()

    
