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
	
#############################BED subFunction##############################
def getUTR(sh_annotation,utr=35):
	'''Extract UTR regions from input bed file (must be 12-column). output is 6-column bed format.
	When utr=35 [default], extract both 5' and 3' UTR. When utr=3, only extract 3' UTR. When utr=5,
	only extract 5' UTR'''
		
	ret_lst=[]
	#bk_annotation=open(fname_annotation, 'r')
	for line in sh_annotation:
		if line.startswith('#'):continue
		if line.startswith('track'):continue
		if line.startswith('browser'):continue
		fields=line.rstrip('\r\n').split()
		chrom=fields[0]
		strand=fields[5]
		txStart=int(fields[1])
		txEnd=int(fields[2])
		cdsStart=int(fields[6])
		cdsEnd=int(fields[7])		
		exon_start=map(int,fields[11].rstrip(',').split(','))
		exon_start=map((lambda x: x + txStart),exon_start)
			
		exon_end=map(int,fields[10].rstrip(',').split(','))
		exon_end=map((lambda x,y:x+y),exon_start,exon_end)
		if (utr==35 or utr==5):
			for st,end in zip(exon_start,exon_end):
				if st < cdsStart:
					utr_st = st
					utr_end = min(end,cdsStart)
					ret_lst.append([chrom,utr_st,utr_end])					
		if (utr==35 or utr==3):
			for st,end in zip(exon_start,exon_end):
				if end > cdsEnd:
					utr_st = max(st, cdsEnd)
					utr_end = end
					ret_lst.append([chrom,utr_st,utr_end])
	#self.f.seek(0)
	#bk_annotation.close()
	return ret_lst
		
def getCDSExon(sh_annotation):	
	'''Extract CDS exon regions from input bed file (must be 12-column).'''		
	ret_lst=[]
	#bk_annotation=open(fname_annotation, 'r')
        for f in sh_annotation:
		f = f.strip().split()
		chrom = f[0]
		chrom_start = int(f[1])
		name = f[4]
		strand = f[5]
		cdsStart = int(f[6])
		cdsEnd = int(f[7])
		blockCount = int(f[9])
		blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
		blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]
                # grab cdsStart - cdsEnd
		cds_exons = []
		cds_seq = ''
		genome_seq_index = []
		for base,offset in zip( blockStarts, blockSizes ):
			if (base + offset) < cdsStart: continue
			if base > cdsEnd: continue
			exon_start = max( base, cdsStart )
			exon_end = min( base+offset, cdsEnd ) 
			#cds_exons.append( (exon_start, exon_end) )
			ret_lst.append([chrom,exon_start,exon_end])	
	#self.f.seek(0)
	#bk_annotation.close()
	return ret_lst
		
def getIntron(sh_annotation):
	ret_lst=[]
	#bk_annotation=open(fname_annotation, 'r')
        for line in sh_annotation:
		try:
			if line.startswith('#'):continue
			if line.startswith('track'):continue
			if line.startswith('browser'):continue   
            	# Parse fields from gene tabls
			fields = line.split()
			chrom     = fields[0]
			tx_start  = int( fields[1] )
			tx_end    = int( fields[2] )
			geneName      = fields[3]
			strand    = fields[5].replace(" ","_")
			cds_start = int( fields[6] )
			cds_end   = int( fields[7] )
			if int(fields[9] ==1):
				continue
    	
			exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
			intron_start = exon_ends[:-1]
			intron_end=exon_starts[1:]
			
			if(strand == '-'):
				intronNum=len(intron_start)
				for st,end in zip(intron_start,intron_end):
					#FO.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t0\t" + strand + '\n')
					#intronNum -= 1
					ret_lst.append([chrom,st,end])	
			else:
				intronNum=1
				for st,end in zip(intron_start,intron_end):
					#FO.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t0\t" + strand + '\n')
					#intronNum += 1
					ret_lst.append([chrom,st,end])
		except:
			print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
			continue
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
	bitsets=0
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
	
#############################BED main processing function##############################
def process_gene_model(sh_annotation):
  	print "start to process UTR CDS_exon Intron"
	utr_3 = getUTR(sh_annotation,utr=3)
	utr_5 = getUTR(sh_annotation,utr=5)
	cds_exon = getCDSExon(sh_annotation)
	intron = getIntron(sh_annotation)
	
	intron = unionBed3(intron)
	cds_exon=unionBed3(cds_exon)
	utr_5 = unionBed3(utr_5)
	utr_3 = unionBed3(utr_3)
	
	utr_5 = subtractBed3(utr_5,cds_exon)
	utr_3 = subtractBed3(utr_3,cds_exon)
	#print cal_size(intron)
	intron = subtractBed3(intron,cds_exon)
	#print cal_size(intron)
	intron = subtractBed3(intron,utr_5)
	intron =subtractBed3(intron,utr_3)

	#print "Processing BED file Done"
	return (cds_exon,intron,utr_5,utr_3)
	
def main():
  
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	#parser = OptionParser(usage,version="%prog " + __version__)
        parser.add_option("-a", "--annotation-file", action = "store", type = "string", dest = "annotation_bed_file",
		  help = "Reference gene model in bed fomat.")
        parser.add_option("-o", "--output-prefix", action = "store", type = "string", dest = "output_prefix",
		  help = "prefix of output file(s). [required]")	
	(options,args)=parser.parse_args()
	
        if (options.annotation_bed_file is None or
            options.output_prefix is None ):
            print prog_base + ": error: missing required command-line argument."
            parser.print_help()
            sys.exit(0)
    
        fname_annotation=options.annotation_bed_file
        bk_annotation=open(fname_annotation)
        sh_annotation=bk_annotation.readlines()
        row_number_annotation=len(sh_annotation)
        bk_annotation.close()
       
	#build bitset
	
        #print "processing %s ..." % (fname_annotation)
	print "Processing BED file Done"
	(cds_exon,intron,utr_5,utr_3) = process_gene_model(sh_annotation)
        utr3_OUT=open(options.output_prefix+".utr_3.txt",'w')
        for chrom,star,end in utr_3:
                print >>utr3_OUT, "%s\t%s\t%s" % (chrom,star,end)
        utr3_OUT.close()

        utr5_OUT=open(options.output_prefix+".utr_5.txt",'w')
        for chrom,star,end in utr_5:
                print >>utr5_OUT, "%s\t%s\t%s" % (chrom,star,end)
        utr5_OUT.close()

        exon_OUT=open(options.output_prefix+".cds_exon.txt",'w')
        for chrom,star,end in cds_exon:
                print >>exon_OUT, "%s\t%s\t%s" % (chrom,star,end)
        exon_OUT.close()

        intron_OUT=open(options.output_prefix+".intron.txt",'w')
        for chrom,star,end in intron:
                print >>intron_OUT, "%s\t%d\t%d" % (chrom,star,end)
        intron_OUT.close()
	
if __name__ == '__main__':
	main()
