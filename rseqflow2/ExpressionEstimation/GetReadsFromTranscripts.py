#!/usr/bin/env python

'''
Created on 2010-2-11

@author: Administrator
'''
#extract the read from transcripts
import os
import sys
import optparse

# Get our basename
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-s", "--transcriptome-sequence", action = "store", type = "string", dest = "raw_fa",
                  help = "raw transcritome reference sequence")
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output",
                  help = "output: cut setting lengths sequences from transcriptome sequences")
parser.add_option("-l", "--readlength", action = "store", type = "string", dest = "readlength",
                  help = "input the read length of the repeat finding")


(options, args) = parser.parse_args()

if (options.raw_fa is None or
    options.output is None or
    options.readlength is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(1)


fname_raw_fa=options.raw_fa
ReadLength=options.readlength
output=options.output

bk_raw=open(fname_raw_fa)
sh_raw=bk_raw.readlines()
row_raw=len(sh_raw)
bk_raw.close()

sequence=''
tran_seq={}
for j in range(0,row_raw):
    if sh_raw[j][0]=='>':
        temp1=sh_raw[j][0:-1].split(' ')
        for k in range(j+1,row_raw):
            sequence +=sh_raw[k][0:-1]
            if k!=row_raw-1:
                if sh_raw[k+1][0]=='>':
                    tran_seq[temp1[0]]=sequence
                    sequence=''
                    break
            tran_seq[temp1[0]]=sequence
Readfile=open(output,'w')

for j in tran_seq.keys():
    seq=tran_seq[j]
    seq_length=len(seq)
    for k in range(1,seq_length-int(ReadLength)+2):
        Read_ID=j+":"+str(k)+"\n"
        Read=seq[k-1:k+int(ReadLength)-1]
        Readfile.writelines(Read_ID)
        Readfile.writelines(Read+'\n')
Readfile.close()
