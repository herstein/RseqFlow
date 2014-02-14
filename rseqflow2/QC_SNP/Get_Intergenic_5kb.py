#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
Check reads distribution over intergenic
-------------------------------------------------------------------------------------------------'''

#import built-in modules
import os,sys
import re
import string
import optparse
import warnings
import string
import collections
import math
import sets

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
#from bx.binned_array import BinnedArray
from bx_extras.fpconst import isNaN
from bx.bitset_utils import *

prog_base = os.path.split(sys.argv[0])[1]
parser = optparse.OptionParser()

def getIntergenic(sh_annotation,direction='up', size=1000):
        '''get intergenic regions. direction=up or down or both.'''

        ret_lst=[]
        #bk_annotation=open(fname_annotation, 'r')
        for line in sh_annotation:
                if line.startswith(('#','track','browser')):
                        continue
                fields = line.split()
                chrom     = fields[0]
                tx_start  = int( fields[1] )
                tx_end    = int( fields[2] )
                strand    = fields[5]
                if(direction=="up" or direction=="both"):
                        if strand=='-':
                                region_st=tx_end
                                region_end=tx_end +size
                        else:
                                region_st = max(tx_start-size,0)
                                region_end=tx_start
                        ret_lst.append([chrom,region_st,region_end])
                if (direction=="down" or direction=="both"):
                        if strand == '-':
                                region_st = max(0,tx_start-size)
                                region_end = tx_start
                        else:
                                region_st = tx_end
                                region_end = tx_end+size
                        ret_lst.append([chrom,region_st,region_end])
        #self.f.seek(0)
        #bk_annotation.close()
        return ret_lst

def unionBed3(lst):
        '''Take the union of 3 column bed files. return a new list'''
        bitsets = binned_bitsets_from_list(lst)
        ret_lst=[]
        for chrom in bitsets:
                bits = bitsets[chrom]
                end = 0
                while 1:
                        start = bits.next_set( end )
                        if start == bits.size: break
                        end = bits.next_clear( start )
                        ret_lst.append([chrom, start, end])
        #bitsets=dict()
        bitsets=0
        #del lst[:]
        return ret_lst

def subtractBed3(lst1,lst2):
        '''subtrack lst2 from lst1'''
        bitsets1 = binned_bitsets_from_list(lst1)
        bitsets2 = binned_bitsets_from_list(lst2)
        ret_lst=[]
        for chrom in bitsets1:
                if chrom not in bitsets1:
                        continue
                bits1 = bitsets1[chrom]
                if chrom in bitsets2:
                        bits2 = bitsets2[chrom]
                        bits2.invert()
                        bits1.iand( bits2 )
                end=0
                while 1:
                        start = bits1.next_set( end )
                        if start == bits1.size: break
                        end = bits1.next_clear( start )
                        ret_lst.append([chrom,start,end])
        bitsets1=dict()
        bitsets2=dict()
        return ret_lst

def list_from_file(fname):
	'''get list from input file'''
	ret_lst=[]
	bk=open(fname,'r')
	for line in bk.xreadlines():
		fields=line[0:-1].split('\t')
		chrom=fields[0]
		start=int(fields[1])
		end=int(fields[2])
		ret_lst.append([chrom,start,end])
	bk.close()
	return ret_lst

def process_gene_model(sh_annotation,cds_exon,intron,utr_3,utr_5):
	intergenic_up_5kb = getIntergenic(sh_annotation,direction="up",size=5000)
        intergenic_down_5kb = getIntergenic(sh_annotation,direction="down",size=5000)

        print "merge integenic region"
        #merge integenic region
        intergenic_up_5kb=unionBed3(intergenic_up_5kb)
        intergenic_down_5kb=unionBed3(intergenic_down_5kb)

        print "purify intergenic region 5kb"
        #purify intergenic region
        intergenic_up_5kb=subtractBed3(intergenic_up_5kb,cds_exon)
        intergenic_up_5kb=subtractBed3(intergenic_up_5kb,utr_5)
        intergenic_up_5kb=subtractBed3(intergenic_up_5kb,utr_3)
        intergenic_up_5kb=subtractBed3(intergenic_up_5kb,intron)
        intergenic_down_5kb=subtractBed3(intergenic_down_5kb,cds_exon)
        intergenic_down_5kb=subtractBed3(intergenic_down_5kb,utr_5)
        intergenic_down_5kb=subtractBed3(intergenic_down_5kb,utr_3)
        intergenic_down_5kb=subtractBed3(intergenic_down_5kb,intron)

	return(intergenic_up_5kb,intergenic_down_5kb)

def main():

        usage="%prog [options]" + '\n' + __doc__ + "\n"
        #parser = OptionParser(usage,version="%prog " + __version__)
        parser.add_option("-a", "--annotation-file", action = "store", type = "string", dest = "annotation_bed_file",
                  help = "Reference gene model in bed fomat.")
	parser.add_option("-p", "--prefix", action = "store", type = "string", dest = "prefix_of_RegionFile",
                  help = "prefix of input Region files: Exon Intron UTR3 UTR5")
        parser.add_option("-o", "--output-prefix", action = "store", type = "string", dest = "output_prefix",
                  help = "prefix of output file(s). [required]")
        (options,args)=parser.parse_args()

        if (options.annotation_bed_file is None or
	    options.prefix_of_RegionFile is None or
            options.output_prefix is None ):
            print prog_base + ": error: missing required command-line argument."
            parser.print_help()
            sys.exit(0)

        fname_annotation=options.annotation_bed_file
	fname_utr3=options.prefix_of_RegionFile+".utr_3.txt"
	fname_utr5=options.prefix_of_RegionFile+".utr_5.txt"
	fname_exon=options.prefix_of_RegionFile+".cds_exon.txt"
	fname_intron=options.prefix_of_RegionFile+".intron.txt"
        bk_annotation=open(fname_annotation)
        sh_annotation=bk_annotation.readlines()
        row_number_annotation=len(sh_annotation)
        bk_annotation.close()
	
	cds_exon=list_from_file(fname_exon)
	intron=list_from_file(fname_intron)
	utr_3=list_from_file(fname_utr3)
	utr_5=list_from_file(fname_utr5)

	print "Processing BED file Done"
        (intergenic_up_5kb,intergenic_down_5kb) = process_gene_model(sh_annotation,cds_exon,intron,utr_3,utr_5)

	inter_up_5kb_OUT=open(options.output_prefix+".intergenic_up_5kb.txt",'w')
        for chrom,star,end in intergenic_up_5kb:
                print >>inter_up_5kb_OUT, "%s\t%s\t%s" % (chrom,star,end)
        inter_up_5kb_OUT.close()

	inter_down_5kb_OUT=open(options.output_prefix+".intergenic_down_5kb.txt",'w')
        for chrom,star,end in intergenic_down_5kb:
                print >>inter_down_5kb_OUT, "%s\t%s\t%s" % (chrom,star,end)
        inter_down_5kb_OUT.close()

if __name__ == '__main__':
        main()
