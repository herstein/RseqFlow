#!/usr/bin/env python
'''
Created on 2009-10-30

@author: Administrator
Modified by J.Herstein 2013-04-18
'''
import sys
import copy
import optparse
import os

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("--c1", action = "store", type = "string", dest = "condition1",
                  help = "specify the name of condition1")
parser.add_option("--c2", action = "store", type = "string", dest = "condition2",
                  help = "specify the name of condition2")
parser.add_option("--n1", action = "store", type = "string", dest = "n_sample1",
                  help = "the number of samples/datasets in condition 1")
parser.add_option("--n2", action = "store", type = "string", dest = "n_sample2",
                  help = "the number of samples/datasets in condition 2")
parser.add_option("-1", action = "store", type = "string", dest = "con1_files",
                  help = "the readcount files of each samples/datasets in condition 1")
parser.add_option("-2", action = "store", type = "string", dest = "con2_files",
                  help = "the readcount files of each samples/datasets in condition 2")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output_c1c2",
                  help = "Combine output of condition1 and condition2")


(options, args) = parser.parse_args()

if (options.condition1 is None or
    options.condition2 is None or
    options.n_sample1 is None or
    options.n_sample2 is None or
    options.con1_files is None or
    options.con2_files is None or
    options.output_c1c2 is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)

con1=options.condition1
con2=options.condition2

n1=int(options.n_sample1)
n2=int(options.n_sample2)

fname1_list=options.con1_files
fname2_list=options.con2_files
outputFilename1=options.output_c1c2

sample1=[]
fname1=fname1_list.split(',')
for i in range(0,len(fname1)):
    temp={}
    sample1.append(temp)

sample2=[]
fname2=fname2_list.split(',')
for i in range(0,len(fname2)):
    temp={}
    sample2.append(temp)


for i in range(0,len(fname1)):
    bk=open(fname1[i])
    sh=bk.readlines()
    for v in range(1,len(sh)):
        temp=sh[v][0:-1].split('\t')
        if temp[3]!='None(RR)':
            sample1[i][temp[0]]=round(float(temp[3]))
    bk.close()

for i in range(0,len(fname2)):
    bk=open(fname2[i])
    sh=bk.readlines()
    for v in range(1,len(sh)):
        temp=sh[v][0:-1].split('\t')
        if temp[3]!='None(RR)':
            sample2[i][temp[0]]=round(float(temp[3]))
    bk.close()


temp1=set(sample1[0].keys())
#for i in range(0,len(fname1)):
for i in range(1,len(fname1)):
    temp2=set(sample1[i].keys())
    temp3=temp1.intersection(temp2)
    temp1=temp3

temp4=set(sample2[0].keys())
#for i in range(0,len(fname2)):
for i in range(1,len(fname2)):
    temp5=set(sample2[i].keys())
    temp6=temp4.intersection(temp5)
    temp4=temp6

gene_list=list(temp1.intersection(temp4))

ws=open(outputFilename1, 'w')

title='gene'
#for i in range(0,n1):
#    title=title+'\t'+con1+'_sample'+str(i+1)
for i in fname1:
    base=os.path.basename(i)
    title=title+'\t'+con1+'_'+ os.path.splitext(base)[0]

#for i in range(0,n2):
#    title=title+'\t'+con2+'_sample'+str(i+1)
for i in fname2:
    base=os.path.basename(i)
    title=title+'\t'+con2+'_'+ os.path.splitext(base)[0]


ws.writelines(title+'\n')


for v in gene_list:
    line=v
    for i in range(0,n1):
        line=line+'\t'+str(int(sample1[i][v]))
    for i in range(0,n2):
        line=line+'\t'+str(int(sample2[i][v]))
    ws.writelines(line+'\n')
ws.close()

    
