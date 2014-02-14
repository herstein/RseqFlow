#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
Check reads distribution over exon, intron, UTR, intergenic ... etc
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

def cal_size(list):
        '''calcualte bed list total size'''
        size=0
        for l in list:
                size += l[2] - l[1]
        return size

def foundone(chrom,ranges, st, end):
        found = 0
        if chrom in ranges:
                found = len(ranges[chrom].find(st,end))
        return found

def build_bitsets(list):
        '''build intevalTree from list'''
        ranges={}
        for l in list:
                chrom =l[0].upper()
                st = int(l[1])
                end = int(l[2])
                if chrom not in ranges:
                        ranges[chrom] = Intersecter()
                else:
                        ranges[chrom].add_interval( Interval( st, end ) )
        return ranges

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

def process_gene_model(cds_exon,intron,utr_3,utr_5,intergenic_up_1kb,intergenic_down_1kb,\
		       intergenic_up_5kb,intergenic_down_5kb,intergenic_up_10kb,intergenic_down_10kb):
	print "build intervalTree"
        #build intervalTree
        cds_exon_ranges = build_bitsets(cds_exon)
        utr_5_ranges = build_bitsets(utr_5)
        utr_3_ranges = build_bitsets(utr_3)
        intron_ranges = build_bitsets(intron)
        interg_ranges_up_1kb_ranges = build_bitsets(intergenic_up_1kb)
        interg_ranges_up_5kb_ranges = build_bitsets(intergenic_up_5kb)
        interg_ranges_up_10kb_ranges = build_bitsets(intergenic_up_10kb)
        interg_ranges_down_1kb_ranges = build_bitsets(intergenic_down_1kb)
        interg_ranges_down_5kb_ranges = build_bitsets(intergenic_down_5kb)
        interg_ranges_down_10kb_ranges = build_bitsets(intergenic_down_10kb)

        exon_size = cal_size(cds_exon)
        intron_size = cal_size(intron)
        utr3_size = cal_size(utr_3)
        utr5_size = cal_size(utr_5)
        int_up1k_size = cal_size(intergenic_up_1kb)
        int_up5k_size = cal_size(intergenic_up_5kb)
        int_up10k_size = cal_size(intergenic_up_10kb)
        int_down1k_size = cal_size(intergenic_down_1kb)
        int_down5k_size = cal_size(intergenic_down_5kb)
        int_down10k_size = cal_size(intergenic_down_10kb)


        print "Processing BED file Done"
        return (cds_exon_ranges,intron_ranges,utr_5_ranges,utr_3_ranges,\
                        interg_ranges_up_1kb_ranges,interg_ranges_up_5kb_ranges,interg_ranges_up_10kb_ranges,\
                        interg_ranges_down_1kb_ranges,interg_ranges_down_5kb_ranges,interg_ranges_down_10kb_ranges,\
                        exon_size,intron_size,utr5_size,utr3_size,\
                        int_up1k_size,int_up5k_size,int_up10k_size,\
                        int_down1k_size,int_down5k_size,int_down10k_size)
def main():

        usage="%prog [options]" + '\n' + __doc__ + "\n"
        #parser = OptionParser(usage,version="%prog " + __version__)
       	parser.add_option("-s", "--sam-file", action = "store", type = "string", dest = "sam_alignment_file",
                  help = "Input alignment file in SAM format")
        parser.add_option("-p", "--prefix", action = "store", type = "string", dest = "prefix_of_RegionFile",
                  help = "prefix of input Region files: Exon Intron UTR3 UTR5")
        parser.add_option("-o", "--output-prefix", action = "store", type = "string", dest = "output_prefix",
                  help = "refix of output file(s). [required]")
        (options,args)=parser.parse_args()

        if (options.sam_alignment_file is None or
            options.prefix_of_RegionFile is None or
            options.output_prefix is None ):
            print prog_base + ": error: missing required command-line argument."
            parser.print_help()
            sys.exit(0)

        fname_sam=options.sam_alignment_file
	bk_sam=open(fname_sam)
        fname_utr3=options.prefix_of_RegionFile+".utr_3.txt"
        fname_utr5=options.prefix_of_RegionFile+".utr_5.txt"
        fname_exon=options.prefix_of_RegionFile+".cds_exon.txt"
        fname_intron=options.prefix_of_RegionFile+".intron.txt"
	fname_int_up1k=options.prefix_of_RegionFile+".intergenic_up_1kb.txt"
	fname_int_down1k=options.prefix_of_RegionFile+".intergenic_down_1kb.txt"
	fname_int_up5k=options.prefix_of_RegionFile+".intergenic_up_5kb.txt"
	fname_int_down5k=options.prefix_of_RegionFile+".intergenic_down_5kb.txt"
	fname_int_up10k=options.prefix_of_RegionFile+".intergenic_up_10kb.txt"
	fname_int_down10k=options.prefix_of_RegionFile+".intergenic_down_10kb.txt"

        

        cds_exon_list=list_from_file(fname_exon)
        intron_list=list_from_file(fname_intron)
        utr_3_list=list_from_file(fname_utr3)
        utr_5_list=list_from_file(fname_utr5)
	int_up1k_list=list_from_file(fname_int_up1k)
	int_down1k_list=list_from_file(fname_int_down1k)
	int_up5k_list=list_from_file(fname_int_up5k)
	int_down5k_list=list_from_file(fname_int_down5k)
	int_up10k_list=list_from_file(fname_int_up10k)
	int_down10k_list=list_from_file(fname_int_down10k)
	
        print "Processing BED file Done"
	(cds_exon_r, intron_r, utr_5_r, utr_3_r,\
        intergenic_up_1kb_r,intergenic_up_5kb_r,intergenic_up_10kb_r,\
        intergenic_down_1kb_r,intergenic_down_5kb_r,intergenic_down_10kb_r,\
        cds_exon_base,intron_base,utr_5_base,utr_3_base,\
        intergenic_up1kb_base,intergenic_up5kb_base,intergenic_up10kb_base,\
        intergenic_down1kb_base,intergenic_down5kb_base,intergenic_down10kb_base)\
        = process_gene_model(cds_exon_list,intron_list,utr_3_list,utr_5_list,int_up1k_list,int_down1k_list,int_up5k_list,int_down5k_list,int_up10k_list,int_down10k_list)
	
	intron_read=0
        cds_exon_read=0
        utr_5_read=0
        utr_3_read=0

        intergenic_up1kb_read=0
        intergenic_down1kb_read=0
        intergenic_up5kb_read=0
        intergenic_down5kb_read=0
        intergenic_up10kb_read=0
        intergenic_down10kb_read=0

        totalReads=0
        totalFrags=0
        unAssignFrags=0

        R_qc_fail=0
        R_duplicate=0
        R_nonprimary=0
        R_unmap=0

        read_paired=1
        read_mapped=2
        read_unmapped=4
        mate_unmapped=8
        read_reverse=16
        mate_reverse=32
        is_read1=64
        is_read2=128
        secondary=256
        qc_fail=512
        duplicate=1024
        _cigar_split=re.compile(r'(\d+)[M|N]')
	
	print "processing " + fname_sam + " ..."
        for v in bk_sam.xreadlines():
                totalReads +=1
                fields=v.split('\t')
                flag=int(fields[1])
                if flag&qc_fail:                        #skip QC fail read
                        R_qc_fail +=1
                        continue
                if flag&duplicate:              #skip duplicate read
                        R_duplicate +=1
                        continue
                if flag&secondary:              #skip non primary hit
                        R_nonprimary +=1
                        continue
                if flag&read_unmapped:          #skip unmap read
                        R_unmap +=1
                        continue
                chrom = fields[2]
                chrom=chrom.upper()
                comb=[int(i) for i in _cigar_split.findall(fields[5])]  #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                #fragment_num += (len(comb) +1)/2
                #blockStart=[]
                #blockSize=[]
                chromStart=string.atoi(fields[3])-1
                exons=[]
                for i in range(0,len(comb),2):
                        ex_start=chromStart + sum(comb[:i])
                        ex_end=ex_start+comb[i]
                        exons.append([chrom, ex_start, ex_end])

                #cigar_str = cigar.list2str(aligned_read.cigar)
                #exons = cigar.fetch_exon(chrom, aligned_read.pos, cigar_str)
                totalFrags += len(exons)

		for exn in exons:
                        #print chrom + '\t' + str(exn[1]) + '\t' + str(exn[2])
                        mid = int(exn[1]) + int((int(exn[2]) - int(exn[1]))/2)
                        if foundone(chrom,cds_exon_r,mid,mid) > 0:
                                cds_exon_read += 1
                                continue
                        elif foundone(chrom,utr_5_r,mid,mid) >0 and foundone(chrom,utr_3_r,mid,mid) == 0:
                                utr_5_read += 1
                                continue
                        elif foundone(chrom,utr_3_r,mid,mid) >0 and foundone(chrom,utr_5_r,mid,mid) == 0:
                                utr_3_read += 1
                                continue
                        elif foundone(chrom,utr_3_r,mid,mid) >0 and foundone(chrom,utr_5_r,mid,mid) > 0:
                                unAssignFrags +=1
                                continue
                        elif foundone(chrom,intron_r,mid,mid) > 0:
                                intron_read += 1
                                continue
                        elif foundone(chrom,intergenic_up_10kb_r,mid,mid) >0 and foundone(chrom,intergenic_down_10kb_r,mid,mid) > 0:
                                unAssignFrags +=1
                                continue
                        elif foundone(chrom,intergenic_up_1kb_r,mid,mid) >0:
                                intergenic_up1kb_read += 1
                                intergenic_up5kb_read += 1
                                intergenic_up10kb_read += 1
                        elif foundone(chrom,intergenic_up_5kb_r,mid,mid) >0:
                                intergenic_up5kb_read += 1
                                intergenic_up10kb_read += 1
                        elif foundone(chrom,intergenic_up_10kb_r,mid,mid) >0:
                                intergenic_up10kb_read += 1

                        elif foundone(chrom,intergenic_down_1kb_r,mid,mid) >0:
                                intergenic_down1kb_read += 1
                                intergenic_down5kb_read += 1
                                intergenic_down10kb_read += 1
                        elif foundone(chrom,intergenic_down_5kb_r,mid,mid) >0:
                                intergenic_down5kb_read += 1
                                intergenic_down10kb_read += 1
                        elif foundone(chrom,intergenic_down_10kb_r,mid,mid) >0:
                                intergenic_down10kb_read += 1
                        else:
                                unAssignFrags +=1
        bk_sam.close()
	
	OUT=open(options.output_prefix+".read_distribution.txt",'w')
        total_AssignFrags=totalFrags-unAssignFrags
        print  >>OUT,"%-30s%.2f" % ("Percentage of Assigned Tags",total_AssignFrags/float(totalFrags))
        print  >>OUT,"====================================================================="
        print  >>OUT,"%-20s%-20s" % ('Region','Percentage of Tags')
        print  >>OUT,"%-10s%8.2f%%" % ('Exons',(cds_exon_read+utr_5_read+utr_3_read)*100/float(total_AssignFrags))
        print  >>OUT,"%-10s%8.2f%%" % ("Introns",intron_read*100/float(total_AssignFrags))

        print  >>OUT,"%-10s%8.2f%%" % ("TSS_up_1kb",intergenic_up1kb_read*100/float(total_AssignFrags))
        print  >>OUT,"%-10s%8.2f%%" % ("TSS_up_5kb",intergenic_up5kb_read*100/float(total_AssignFrags))
        print  >>OUT,"%-10s%8.2f%%" % ("TSS_up_10kb",intergenic_up10kb_read*100/float(total_AssignFrags))
        print  >>OUT,"%-10s%8.2f%%" % ("TES_down_1kb",intergenic_down1kb_read*100/float(total_AssignFrags))
        print  >>OUT,"%-10s%8.2f%%" % ("TES_down_5kb",intergenic_down5kb_read*100/float(total_AssignFrags))
        print  >>OUT,"%-10s%8.2f%%" % ("TES_down_10kb",intergenic_down10kb_read*100/float(total_AssignFrags))
        print  >>OUT,"====================================================================="

        print  >>OUT,"Note:"
        print  >>OUT,"If reads spliced once, it will be counted as 2 tags"
        print  >>OUT,"TSS= Transcription Start Site, TES=Transcription End Site, down=downstream, up=upstream"

        OUT.close()

if __name__ == '__main__':
        main()
