#!/bin/bash
USAGE=$'Usage:

ExpressionEstimation.sh {-f <r1.fastq.gz> or -1 <r1.fastq.gz> -2 <r2.fastq.gz>} [options]

Required arugments:

-f/--fastq <reads.fastq.gz>            fastq or fastq.gz file for single end data
-- OR --
-1/--read1 <reads.fastq.gz>            the first read fastq file for paired end data
-2/--read2 <reads.fastq.gz>            the second read fastq file for paired end data
-a/--annotation <annotation.gtf>    reference annotation in GTF format
-o/--output-prefix <outputPrefix>   prefix of output files


Optional arguments:

-c/--transcriptome [ ref.fa ]       transcriptome reference sequences
--mp [ MAX,MIN ]                    both integers. bowtie2 option to set the mistch penalty, default is 6,2
--score-min [ Function,a,b ]        bowtie2 option to set a function of read length for the minimum alignment score
                                    necessary for an alignment to be considered valid.
                                    Default is L,0,-0.6 which is defined as Linear function:
                                    f(x) = 0 + -0.6 * x  where x is the read length.
                                    Available function types are constant (C), linear (L), square-root (S), and
                                    natural log (G). The parameters are specified as F,B,A - the function type,
                                    the constant term, and the coefficient separated by commas with no whitespace.
                                    The constant term and coefficient may be negative and/or floating-point numbers.
                                    For more info see the Bowtie2 manual.
--tSam [ alignment.sam ]            alignments to transcriptome in SAM format
--cleanup                           delete temporary files
-h/--help                           print this usage message'

## 2013-06-11 -- "unique" is the only option offered right now but "proportion" and "random" may be used at a later time
## Bowtie2 takes too long to run in order to obtain multiple mappings that make the random and proportion options worthwhile, hence they are being removed for the time being.
# Removed proportion and random from Usage for now, unique is the default, but bowtie2 maps reads using the "best" alignment so 
# it doesn't guarantee uniqueness, it means only that it takes the best alignment which may or not be unique. 
# The "XS:" field in the samfile will be present if a read was multimapped and absent if a read was uniquely mapped
#-u/--unique                         use data of reads mapped to only one gene.
#-p/--proportion                     assign reads according to proportion of genes expression level.
#-r/--random                         assign reads randomly.


echo ""
echo "You are running: $VERSION"
echo ""

if [ $# -eq 0 ]; then
        echo "No arguments or options!"
	echo "$USAGE"
        exit 1
fi

is_unique=true
is_proportion=false
is_random=false
is_cleanup=false
MP='6,2'
Score_Min='L,0,-0.6'
TOOL='bowtie2'
OUTPUT='Expression_output'

declare -a ARGS

ARGS=($@)

for ((i=0;i<$#;++i))
do
        if [[ ${ARGS[i]} = '-f' || ${ARGS[i]} = '--fastq' ]]; then # Reads Input Fastq file (single end data)
                FASTQ=${ARGS[(i+1)]}
                i=($i+1)
	elif [[ ${ARGS[i]} = '-1' || ${ARGS[i]} = '--read1' ]]; then # Reads Input Fastq file (paired end data)
                READ1=${ARGS[(i+1)]}
                i=($i+1)
	elif [[ ${ARGS[i]} = '-2' || ${ARGS[i]} = '--read2' ]]; then # Reads Input Fastq file (paired end data)
                READ2=${ARGS[(i+1)]}
                i=($i+1)
        elif [[ ${ARGS[i]} = '-c' || ${ARGS[i]} = '--transcriptome' ]]; then # Transcriptome Reference Sequences
                TRANSCRIPTOME=${ARGS[(i+1)]}
                i=($i+1)
        elif [[ ${ARGS[i]} = '-a' || ${ARGS[i]} = '--annotation' ]]; then # Reference Annotation
                ANNOTATION=${ARGS[(i+1)]}
                i=($i+1)
        elif [[ ${ARGS[i]} = '-o' || ${ARGS[i]} = '--output-prefix' ]]; then # Output file Prefix
                OUTPUT=${ARGS[(i+1)]}
                i=($i+1)
        elif [[ ${ARGS[i]} = '--mp' ]]; then # max mismatch penalty
                MP=${ARGS[(i+1)]}
                i=($i+1)
	elif [[ ${ARGS[i]} = '--score-min' ]]; then # minimum alignment score
                Score_Min=${ARGS[(i+1)]}
                i=($i+1)        
## 2013-06-11 -- "unique" is the only option offered right now but "proportion" and "random" may be used at a later time
	elif [[ ${ARGS[i]} = '-u' || ${ARGS[i]} = '--unique' ]]; then # Just base on reads uniquely mapped to one gene
                is_unique=true
#        elif [[ ${ARGS[i]} = '-p' || ${ARGS[i]} = '--proportion' ]]; then # Assign Reads base on Gene Expression Level Proportion
#                is_proportion=true
#	elif [[ ${ARGS[i]} = '-r' || ${ARGS[i]} = '--random' ]]; then # Assign Reads Randomly
#                is_random=true
	elif [  ${ARGS[i]} = '--tSam'  ]; then #input transcriptome sam file
                TrantoSam=${ARGS[(i+1)]}
                i=($i+1)
	elif [  ${ARGS[i]} = '--cleanup'  ]; then
		is_cleanup=true
	elif [[ ${ARGS[i]} = '-h' || ${ARGS[i]} = '--help' ]]; then # Help Information
		echo "$USAGE"
             	exit 0
        else
                UNKNOWN=${ARGS[i]}
                echo "No switch encountered: $UNKNOWN"
		echo "$USAGE"
                exit 1
        fi
done
########################################Check Required Arguments#########################################
if [[ -z $TrantoSam ]]; then
	if [[ -z $FASTQ && -z $READ1 && -z $READ2 || -z $TRANSCRIPTOME || -z $ANNOTATION ]]; then
        	echo "Error: required input files not specified!"
		echo "$USAGE"
        	exit 1
	fi
	if [[ -n $FASTQ && -n $READ1 ]] || [[ -n $FASTQ && -n $READ2 ]]; then
                echo "Error: single end data and paired end data can't be given together!"
                echo "If you want to run with single end data, try option: -f <reads.fastq.gz>"
                echo "If you want to run with paired end data, try options: -1 <read1.fastq.gz> -2 <read2.fastq.gz>"
                exit 1
        fi

        if [ -z $FASTQ ] && [[ -z $READ1 || -z $READ2 ]]; then
                echo "Error: for paired end data, read1 and read2 must be given together!"
                echo "If you want to run with paired end data, try options: -1 <read1.fastq.gz> -2 <read2.fastq.gz>"
                exit 1
        fi
else 
	if [[ -z $ANNOTATION ]]; then
		echo "Error: you input a sam file without the annotation file!"
		echo "Try option: -a <annotation.gtf>"
		exit 1
	fi
fi

if [ -z $OUTPUT ]; then
        echo "Error: missing output prefix! Please give the prefix of output files."
	echo "Try this option to specify the output prefix: - --output-prefix < out.prefix >" 
        exit 1
fi
##########################################Check Method##########################################
is_only_alignment=false

## 2013-06-11 -- "unique" is the only option offered right now but "proportion" and "random" may be used at a later time
#if  ! ($is_unique || $is_proportion || $is_random); then
#	echo "Warning: you didn't choose any analysis to do, so only alignment will be done."
#	echo "try one of the following options:"
#	echo "-u/--unique     #Description: use data of reads mapped to only one gene."
#        echo "-p/--proportion #Description: assign reads according to proportion of genes expression level."
#        echo "-r/--random     #Description: assign reads randomly."
#        is_only_alignment=true
#        #exit 1
#fi

########################################Check input files#######################################
echo "Checking input files..."

if [[ -z $TrantoSam ]]; then
	echo "Check_for_reference_annotation_withoutGenome.py -t $TRANSCRIPTOME -a $ANNOTATION"
	Check_for_reference_annotation_withoutGenome.py -t $TRANSCRIPTOME -a $ANNOTATION
	ERR=$?
	if [ $ERR -ne 0 ]; then
             	echo "Error: input files failed to pass the check. There are some errors in the input files, please check them!"
               	exit 1
        fi
#fi

        if [ -z $FASTQ ]; then
	     if [[ -n $READ1 ]]; then
        	     echo "Check_for_reads_file.py -r $READ1"
        	     Check_for_reads_file.py -r $READ1
		     ERR=$?
		     if [ $ERR -ne 0 ]; then
                	     echo "Error: input files failed to pass the check. There are some errors in the input files, please check them!"
                	     exit 1
        	     fi
        	     echo "Check_for_reads_file.py -r $READ2"
        	     Check_for_reads_file.py -r $READ2
		     ERR=$?
		     if [ $ERR -ne 0 ]; then
                	     echo "Error: input files failed to pass the check. There are some errors in the input files, please check them!"
                	     exit 1
        	     fi
	     fi
        else
             echo "Check_for_reads_file.py -r $FASTQ"
             Check_for_reads_file.py -r $FASTQ
	     ERR=$?
	     if [ $ERR -ne 0 ]; then
                     echo "Error: input files failed to pass the check. There are some errors in the input files, please check them!"
                     exit 1
             fi
        fi
fi
#####################################################################Alignment#####################################################################
if [[ -n $TrantoSam ]]; then
	echo "You input a sam file, so alignment will be skipped."
	echo "Check_for_transcriptomeSam.py -s $TrantoSam"
	Check_for_transcriptomeSam.py -s $TrantoSam
        ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Error: input files failed to pass the check. There are some errors in the input files, please check them!"
                exit 1
        fi
else 
        #--------------------------------------------------Bowtie2:Alignment to Transcriptome------------------------------------------------------#
	TranId=${TRANSCRIPTOME%.fa} #get full path of filename without extension
	TrantoIndexSamp=${TranId##*/} #get basename of file
        TrantoSam="$OUTPUT"_Bowtie2_transcriptome.sam

        # Check if bowtie2 indexes already exist in reference path, if not run bowtie2-build
        # Can keep reference index so it doesn't need to be created every time if mapping multiple samples.
	
	echo "Checking for existing bowtie2 indexes for $TRANSCRIPTOME"
        if [[ -f $TranId.1.bt2 && -f $TranId.2.bt2 && -f $TranId.3.bt2 && -f $TranId.4.bt2 && -f $TranId.rev.1.bt2 && -f $TranId.rev.2.bt2 ]]; then
	    echo "Existing bowtie2 indexes found for $TRANSCRIPTOME"
            TrantoIndex=$TranId
        else
	    echo "Start building bowtie2 indexes for $TRANSCRIPTOME" 
            TrantoIndex=$TrantoIndexSamp
	    echo "bowtie2-build $TRANSCRIPTOME $TrantoIndex"
            bowtie2-build $TRANSCRIPTOME $TrantoIndex
            ERR=$?
            if [ $ERR -ne 0 ]; then
                echo "Errors when building bowtie2 index for $TRANSCRIPTOME, stopping the pipeline!"
		if [ -f $TrantoIndexSamp ]; then
                    rm $TrantoIndexSamp* -f
                fi
                exit 1
            fi
        fi

        echo "Starting alignment with bowtie2"
	# Output: $TrantoSam ($OUTPUT"_Bowtie2_transcriptome.sam)
        if [[ -n $FASTQ ]]; then
		echo "bowtie2 -x $TrantoIndex -U $FASTQ --mp $MP --score-min $Score_Min --sam-no-hd -S $TrantoSam"
                bowtie2 -x $TrantoIndex -U $FASTQ --mp $MP --score-min $Score_Min --sam-no-hd -S $TrantoSam
        else
		echo "bowtie2 -x $TrantoIndex -1 $READ1 -2 $READ2 --mp $MP --score-min $Score_Min --sam-no-hd -S $TrantoSam"
                bowtie2 -x $TrantoIndex -1 $READ1 -2 $READ2 --mp $MP --score-min $Score_Min --sam-no-hd -S $TrantoSam
        fi

        if [ $ERR -ne 0 ]; then                
		echo "Errors when running bowtie2, stopping the pipeline!"
                #rm $TrantoIndexSamp* -f
                exit 1
        fi
	if [[ $is_cleanup && -f $TrantoIndexSamp ]]; then
		rm $TrantoIndexSamp* -f
        fi

        echo "bowtie2 has finished!"
fi

if $is_only_alignment; then
        exit 0
fi
#############################################################Split input files by chrom############################################################
echo "Split files by chromosome..."
#-------split sam file------#
if [[ -n $TrantoSam ]]; then
	# Output: $OUTPUT"_chrList_sam.txt
	echo "SplitByChromosome_for_transcriptomeSamFile.py -i $TrantoSam -p $OUTPUT -o "$OUTPUT"_chrList_sam.txt"
	SplitByChromosome_for_transcriptomeSamFile.py -i $TrantoSam -p $OUTPUT -o "$OUTPUT"_chrList_sam.txt
	ERR=$?
	if [ $ERR -ne 0 ]; then
                echo "Errors when running SplitByChromosome_for_transcriptomeSamFile.py, stopping the pipeline!"
                exit 1
        fi
else
	echo "$TrantoSam: file not found, stopping the pipeline!"
	exit 1
fi

#-------split fa file-------#
# Reference file is being split up, not the sample fasta
#if [[ -n $TRANSCRIPTOME ]]; then
#	# Output: $OUTPUT"_chrList_fa.txt
#	echo "SplitByChromosome_for_transcriptomeSequenceFaFile.py -i $TRANSCRIPTOME -p $OUTPUT -o "$OUTPUT"_chrList_fa.txt"
#	SplitByChromosome_for_transcriptomeSequenceFaFile.py -i $TRANSCRIPTOME -p $OUTPUT -o "$OUTPUT"_chrList_fa.txt
#	ERR=$?
#	if [ $ERR -ne 0 ]; then
#                echo "Errors when running SplitByChromosome_for_transcriptomeSequenceFaFile.py, stopping the pipeline!"
#                exit 1
#        fi
#else
#	echo "$TRANSCRIPTOME: file not found, stopping the pipeline!"
#        exit 1
#fi

#------split gtf file-------#
if [[ -n $ANNOTATION ]]; then
	# Output: $OUTPUT"_chrList_gtf.txt
	echo "SplitByChromosome_for_annotationGtfFile.py -i $ANNOTATION -p $OUTPUT -o "$OUTPUT"_chrList_gtf.txt"
 	SplitByChromosome_for_annotationGtfFile.py -i $ANNOTATION -p $OUTPUT -o "$OUTPUT"_chrList_gtf.txt
	ERR=$?
	if [ $ERR -ne 0 ]; then
                echo "Errors when running SplitByChromosome_for_annotationGtfFile.py, stopping the pipeline!"
                exit 1
        fi
else
        echo "$ANNOTATION: file not found, stopping the pipeline!"
        exit 1
fi
############################################################Build Index and combination############################################################

chrList=`cat "$OUTPUT"_chrList_sam.txt`
MAP_LOG="$OUTPUT"_mapping_info.log
i=0
for chr in $chrList
do
	AN[$i]="$OUTPUT"_"$chr"_annotation.gtf
	#TS[$i]="$OUTPUT"_"$chr"_sequence.fa
	EC[$i]="$OUTPUT"_ExonCombination_"$chr".txt
	EI[$i]="$OUTPUT"_ExonIndex_"$chr".txt
	JC[$i]="$OUTPUT"_JunctionCombination_"$chr".txt
	JI[$i]="$OUTPUT"_JunctionIndex_"$chr".txt
	CSAM[$i]="$OUTPUT"_"$chr"_alignment.sam
	USAM[$i]="$OUTPUT"_"$chr"_uniqGene.sam
	MUL[$i]="$OUTPUT"_"$chr"_reads_multiGene.txt

	GEU[$i]="$OUTPUT"_"$chr"_GeneExpressionLevel_unique.txt
	EEU[$i]="$OUTPUT"_"$chr"_ExonExpressionLevel_unique.txt
	JEU[$i]="$OUTPUT"_"$chr"_JunctionExpressionLevel_unique.txt

	#RFG[$i]="$OUTPUT"_ReadsFromGene_"$chr".fa
	#TRR[$i]="$OUTPUT"_RR_trans_"$chr".txt
	#GRR[$i]="$OUTPUT"_RR_chrs_"$chr".txt
	#EL[$i]="$OUTPUT"_ExonLength_"$chr".txt
	#RTSAM[$i]="$OUTPUT"_transcriptome_ReadsFromGene_"$chr".sam
	#RGSAM[$i]="$OUTPUT"_genome_ReadsFromGene_"$chr".sam

## 2013-06-11 -- "unique" is the only option offered right now but "proportion" and "random" may be used at a later time
#	GEP[$i]="$OUTPUT"_"$chr"_GeneExpressionLevel_proportion.txt
#        EEP[$i]="$OUTPUT"_"$chr"_ExonExpressionLevel_proportion.txt
#        JEP[$i]="$OUTPUT"_"$chr"_JunctionExpressionLevel_proportion.txt
#	GEP_Merge[$i]="$OUTPUT"_"$chr"_GeneExpressionLevel_proportion_merge.txt
#        EEP_Merge[$i]="$OUTPUT"_"$chr"_ExonExpressionLevel_proportion_merge.txt
#        JEP_Merge[$i]="$OUTPUT"_"$chr"_JunctionExpressionLevel_proportion_merge.txt

#	GER[$i]="$OUTPUT"_"$chr"_GeneExpressionLevel_random.txt
#        EER[$i]="$OUTPUT"_"$chr"_ExonExpressionLevel_random.txt
#        JER[$i]="$OUTPUT"_"$chr"_JunctionExpressionLevel_random.txt
#	GER_Merge[$i]="$OUTPUT"_"$chr"_GeneExpressionLevel_random_merge.txt
#        EER_Merge[$i]="$OUTPUT"_"$chr"_ExonExpressionLevel_random_merge.txt
#        JER_Merge[$i]="$OUTPUT"_"$chr"_JunctionExpressionLevel_random_merge.txt

	#CHR[$i]=$chr

	i=$i+1
done

if $is_cleanup; then
#	rm "$OUTPUT"_chrList_sam.txt "$OUTPUT"_chrList_fa.txt "$OUTPUT"_chrList_gtf.txt -f
	rm "$OUTPUT"_chrList_sam.txt "$OUTPUT"_chrList_gtf.txt -f
fi

l=${#AN[@]}

echo "Build Exon and Junction Index..."
#-------------build Exon Index-------------#
for ((i=0;i<l;++i))
do
	# Output: ${EC[i]}  ("$OUTPUT"_ExonCombination_"$chr".txt)
	echo "ExonCombination.py -g ${AN[i]} -o ${EC[i]}"
	ExonCombination.py -g ${AN[i]} -o ${EC[i]} 
	ERR=$?
	if [ $ERR -ne 0 ]; then
                echo "Errors when running ExonCombination.py, stopping the pipeline!"
                exit 1
        fi

	# Output: ${EI[i]}  ("$OUTPUT"_ExonIndex_"$chr".txt)
	echo "ExonIndex.py -g ${AN[i]} -e ${EC[i]} -o ${EI[i]}"
	ExonIndex.py -g ${AN[i]} -e ${EC[i]} -o ${EI[i]}
	ERR=$?
	if [ $ERR -ne 0 ]; then
                echo "Errors when running ExonIndex.py, stopping the pipeline!"
                exit 1
        fi
done
#-----------build Junction Index-----------#
for ((i=0;i<l;++i))
do
	# Output: ${JC[i]} ("$OUTPUT"_JunctionCombination_"$chr".txt)
	echo "JunctionCombination.py -g ${AN[i]} -o ${JC[i]}"
	JunctionCombination.py -g ${AN[i]} -o ${JC[i]}
	ERR=$?
	if [ $ERR -ne 0 ]; then
                echo "Errors when running JunctionCombination.py, stopping the pipeline!"
                exit 1
        fi

	# Output: ${JI[i]} ("$OUTPUT"_JunctionIndex_"$chr".txt)
	echo "JunctionIndex.py -g ${AN[i]} -j ${JC[i]} -o ${JI[i]}"
	JunctionIndex.py -g ${AN[i]} -j ${JC[i]} -o ${JI[i]}
	ERR=$?
	if [ $ERR -ne 0 ]; then
                echo "Errors when running JunctionIndex.py, stopping the pipeline!"
                exit 1
        fi
done

###################Split Result, one is based on uniquely mapped reads and another is based on multiple mapped reads###############################
echo "start to get reads mapping information"
# Output: $MAP_LOG ("$OUTPUT"_mapping_info.log)
echo "Get_ReadsMappingInformation.py -s $TrantoSam -l $MAP_LOG"
Get_ReadsMappingInformation.py -s $TrantoSam -l $MAP_LOG
echo "Start to split alignment results..."
for ((i=0;i<l;++i))
do
	# Output: ${USAM[i]} and ${MUL[i]}  ("$OUTPUT"_"$chr"_uniqGene.sam, "$OUTPUT"_"$chr"_reads_multiGene.txt)
	echo "SamSplitEvenly_and_Randomly_gencode_modify.py -s ${CSAM[i]} -g ${AN[i]} -u ${USAM[i]} -m ${MUL[i]}"
        SamSplitEvenly_and_Randomly_gencode_modify.py -s ${CSAM[i]} -g ${AN[i]} -u ${USAM[i]} -m ${MUL[i]}
        ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running SamSplitEvenly_and_Randomly_gencode_modify.py, stopping the pipeline!"
                exit 1
        fi
done

if $is_cleanup; then
	rm ${AN[@]} -f
	rm ${CSAM[@]} -f
fi
##############################################################UniqueMap ExpressionLevel############################################################
if ( $is_unique || $is_proportion || $is_random ); then
	echo "Start to estimate expression level by using reads uniquely mapped to only one gene..."
	for ((i=0;i<l;++i))
	do
		# Output: ${GEU[i]} ("$OUTPUT"_"$chr"_GeneExpressionLevel_unique.txt)
		echo "GeneExpressionLevel.py -u ${USAM[i]} -i ${EI[i]} -c ${EC[i]} -l ${MAP_LOG} -o ${GEU[i]}"
		GeneExpressionLevel.py -u ${USAM[i]} -i ${EI[i]} -c ${EC[i]} -l ${MAP_LOG} -o ${GEU[i]}
		ERR=$?
		if [ $ERR -ne 0 ]; then
        	        echo "Errors when running GeneExpressionLevel.py, stopping the pipeline!"
               		exit 1
        	fi

                # Output: ${EEU[i]} ("$OUTPUT"_"$chr"_ExonExpressionLevel_unique.txt)
		echo "ExonExpressionLevel.py -u ${USAM[i]} -i ${EI[i]} -l ${MAP_LOG} -o ${EEU[i]}"
		ExonExpressionLevel.py -u ${USAM[i]} -i ${EI[i]} -l ${MAP_LOG} -o ${EEU[i]}
		ERR=$?
		if [ $ERR -ne 0 ]; then
                	echo "Errors when running ExonExpressionLevel.py, stopping the pipeline!"
                	exit 1
        	fi

                # Output: ${JEU[i]} ("$OUTPUT"_"$chr"_JunctionExpressionLevel_unique.txt)
		echo "JunctionExpressionLevel.py -u ${USAM[i]} -i ${JI[i]} -l ${MAP_LOG} -o ${JEU[i]}"
		JunctionExpressionLevel.py -u ${USAM[i]} -i ${JI[i]} -l ${MAP_LOG} -o ${JEU[i]}
		ERR=$?
		if [ $ERR -ne 0 ]; then
                	echo "Errors when running JunctionExpressionLevel.py, stopping the pipeline!"
                	exit 1
        	fi
	done
	if $is_unique; then

		GEU_Merge_Whole="$OUTPUT"_whole_GeneExpressionLevel_unique.txt
        	EEU_Merge_Whole="$OUTPUT"_whole_ExonExpressionLevel_unique.txt
        	JEU_Merge_Whole="$OUTPUT"_whole_JunctionExpressionLevel_unique.txt
	
		#echo "#Unique method: remove the multi-mapped reads directly" >$GEU_Merge_Whole
		grep 'Strand' ${GEU[0]} >>$GEU_Merge_Whole
		grep -v 'Strand' -h ${GEU[@]} >>$GEU_Merge_Whole
		#echo "#Unique method: remove the multi-mapped reads directly" >$EEU_Merge_Whole
		grep 'Strand' ${EEU[0]} >>$EEU_Merge_Whole
		grep -v 'Strand' -h ${EEU[@]} >>$EEU_Merge_Whole
		#echo "#Unique method: remove the multi-mapped reads directly" >$JEU_Merge_Whole
		grep 'Strand' ${JEU[0]} >>$JEU_Merge_Whole
		grep -v 'Strand' -h ${JEU[@]} >>$JEU_Merge_Whole
		echo "Expression Level Estimation (Unique Method) Done!"
	fi
fi
#####################################################MultipleMap ExpressionLevel and Merge two parts###############################################

## 2013-06-11 -- "unique" is the only option offered right now but "proportion" and "random" may be used at a later time

#if $is_proportion; then
#	echo "Start to estimation expression level by assigning reads according to proportion of genes expression level..."
#	GEP_Merge_Whole="$OUTPUT"_whole_GeneExpressionLevel_proportion.txt
#	EEP_Merge_Whole="$OUTPUT"_whole_ExonExpressionLevel_proportion.txt
#	JEP_Merge_Whole="$OUTPUT"_whole_JunctionExpressionLevel_proportion.txt
#	for ((i=0;i<l;++i))
#        do
#		# Output: ${GEP[i]} ("$OUTPUT"_"$chr"_GeneExpressionLevel_proportion.txt)
#		echo "GeneExpressionLevel_proportionAssign.py -m ${MUL[i]} -i ${EI[i]} -c ${EC[i]} -u ${GEU[i]} -l $MAP_LOG -o ${GEP[i]}"
#		GeneExpressionLevel_proportionAssign.py -m ${MUL[i]} -i ${EI[i]} -c ${EC[i]} -u ${GEU[i]} -l $MAP_LOG -o ${GEP[i]}
#		ERR=$?
#                if [ $ERR -ne 0 ]; then
#                        echo "Errors when running GeneExpressionLevel_proportionAssign.py, stopping the pipeline!"
#                        exit 1
#                fi
#
#		# Output: ${EEP[i]} ("$OUTPUT"_"$chr"_ExonExpressionLevel_proportion.txt)
#		echo "ExonExpressionLevel_proportionAssign.py -m ${MUL[i]} -i ${EI[i]} -u ${GEU[i]} -l $MAP_LOG -o ${EEP[i]}"
#		ExonExpressionLevel_proportionAssign.py -m ${MUL[i]} -i ${EI[i]} -u ${GEU[i]} -l $MAP_LOG -o ${EEP[i]}
#		ERR=$?
#		if [ $ERR -ne 0 ]; then
#                        echo "Errors when running ExonExpressionLevel_proportionAssign.py, stopping the pipeline!"
#                        exit 1
#                fi
#
#		# Output: ${JEP[i]} ("$OUTPUT"_"$chr"_JunctionExpressionLevel_proportion.txt)
#		echo "JunctionExpressionLevel_proportionAssign.py -m ${MUL[i]} -i ${JI[i]} -u ${GEU[i]} -l $MAP_LOG -o ${JEP[i]}"
#		JunctionExpressionLevel_proportionAssign.py -m ${MUL[i]} -i ${JI[i]} -u ${GEU[i]} -l $MAP_LOG -o ${JEP[i]}
#		ERR=$?
#                if [ $ERR -ne 0 ]; then
#                        echo "Errors when running JunctionExpressionLevel_proportionAssign.py, stopping the pipeline!"
#                        exit 1
#                fi
#
#		# Output: ${GEP_Merge[i]} ("$OUTPUT"_"$chr"_GeneExpressionLevel_proportion_merge.txt)
#		echo "Merge_unique_mulitple.py -u ${GEU[i]} -m ${GEP[i]} -o ${GEP_Merge[i]} -t gene"
#       		Merge_unique_mulitple.py -u ${GEU[i]} -m ${GEP[i]} -o ${GEP_Merge[i]} -t gene
#
#		# Output: ${EEP_Merge[i]} ("$OUTPUT"_"$chr"_ExonExpressionLevel_proportion_merge.txt)
#		echo "Merge_unique_mulitple.py -u ${EEU[i]} -m ${EEP[i]} -o ${EEP_Merge[i]} -t exon"
#		Merge_unique_mulitple.py -u ${EEU[i]} -m ${EEP[i]} -o ${EEP_Merge[i]} -t exon
#
#		# Output: ${JEP_Merge[i]} ("$OUTPUT"_"$chr"_JunctionExpressionLevel_proportion_merge.txt)
#		echo "Merge_unique_mulitple.py -u ${JEU[i]} -m ${JEP[i]} -o ${JEP_Merge[i]} -t junction"
#		Merge_unique_mulitple.py -u ${JEU[i]} -m ${JEP[i]} -o ${JEP_Merge[i]} -t junction
#	done
#
#	#echo "#Proportion method: assigns the multi-mapped reads according to the proportion of gene expression level" >$GEP_Merge_Whole
#	grep 'Strand' ${GEP_Merge[0]} >>$GEP_Merge_Whole
#	grep -v 'Strand' -h ${GEP_Merge[@]} >>$GEP_Merge_Whole
#	#echo "#Proportion method: assigns the multi-mapped reads according to the proportion of gene expression level" >$EEP_Merge_Whole
#	grep 'Strand' ${EEP_Merge[0]} >>$EEP_Merge_Whole
#	grep -v 'Strand' -h ${EEP_Merge[@]} >>$EEP_Merge_Whole
#	#echo "#Proportion method: assigns the multi-mapped reads according to the proportion of gene expression level" >$JEP_Merge_Whole
#	grep 'Strand' ${JEP_Merge[0]} >>$JEP_Merge_Whole
#	grep -v 'Strand' -h ${JEP_Merge[@]} >>$JEP_Merge_Whole
#	
#	if $is_cleanup; then
#		rm ${GEP[@]} ${EEP[@]} ${JEP[@]} -f
#		rm ${GEP_Merge[@]} ${EEP_Merge[@]} ${JEP_Merge[@]} -f
#	fi
#	echo "Expression Level Estimation (Proportion Method) Done!"
#fi
#
#if $is_random; then
#	echo "Start to estimation expression level by assigning reads randomly..."
#	GER_Merge_Whole="$OUTPUT"_whole_GeneExpressionLevel_random.txt
#        EER_Merge_Whole="$OUTPUT"_whole_ExonExpressionLevel_random.txt
#        JER_Merge_Whole="$OUTPUT"_whole_JunctionExpressionLevel_random.txt
#	for ((i=0;i<l;++i))
#	do
#		# Output: ${GER[i]} and ${EER[i]} ("$OUTPUT"_"$chr"_GeneExpressionLevel_random.txt, "$OUTPUT"_"$chr"_ExonExpressionLevel_random.txt)
#		echo "Gene_Exon_ExpressionLevel_randomAssign.py -m ${MUL[i]} -i ${EI[i]} -c ${EC[i]}  -l ${MAP_LOG} -g ${GER[i]} -e ${EER[i]}"
#		Gene_Exon_ExpressionLevel_randomAssign.py -m ${MUL[i]} -i ${EI[i]} -c ${EC[i]}  -l ${MAP_LOG} -g ${GER[i]} -e ${EER[i]}
#		ERR=$?
#                if [ $ERR -ne 0 ]; then
#                        echo "Errors when running Gene_Exon_ExpressionLevel_randomAssign.py, stopping the pipeline!"
#                        exit 1
#                fi
#
#                # Output: ${JER[i]} ("$OUTPUT"_"$chr"_JunctionExpressionLevel_random.txt)
#		echo "JunctionExpressionLevel_randomAssign.py -m ${MUL[i]} -i ${JI[i]}  -l ${MAP_LOG} -j ${JER[i]}"
#		JunctionExpressionLevel_randomAssign.py -m ${MUL[i]} -i ${JI[i]}  -l ${MAP_LOG} -j ${JER[i]}
#		ERR=$?
#                if [ $ERR -ne 0 ]; then
#                        echo "Errors when running JunctionExpressionLevel_randomAssign.py, stopping the pipeline!"
#                        exit 1
#                fi
#
#		# Output: ${GER_Merge[i]} ("$OUTPUT"_"$chr"_GeneExpressionLevel_random_merge.txt)
#		echo "Merge_unique_mulitple.py -u ${GEU[i]} -m ${GER[i]} -o ${GER_Merge[i]} -t gene"
#		Merge_unique_mulitple.py -u ${GEU[i]} -m ${GER[i]} -o ${GER_Merge[i]} -t gene
#
#		# Output: ${EER_Merge[i]} ("$OUTPUT"_"$chr"_ExonExpressionLevel_random_merge.txt)
#		echo "Merge_unique_mulitple.py -u ${EEU[i]} -m ${EER[i]} -o ${EER_Merge[i]} -t exon"
#        	Merge_unique_mulitple.py -u ${EEU[i]} -m ${EER[i]} -o ${EER_Merge[i]} -t exon
#
#		# Output: ${JER_Merge[i]} ("$OUTPUT"_"$chr"_JunctionExpressionLevel_random_merge.txt)
#		echo "Merge_unique_mulitple.py -u ${JEU[i]} -m ${JER[i]} -o ${JER_Merge[i]} -t junction"
#	        Merge_unique_mulitple.py -u ${JEU[i]} -m ${JER[i]} -o ${JER_Merge[i]} -t junction
#	done
#        
#	#echo "#Random method: assigns the multi-mapped reads randomly" >$GER_Merge_Whole
#	grep 'Strand' ${GER_Merge[0]} >>$GER_Merge_Whole
#	grep -v 'Strand' -h ${GER_Merge[@]} >>$GER_Merge_Whole
#	#echo "#Random method: assigns the multi-mapped reads randomly" >$EER_Merge_Whole
#	grep 'Strand' ${EER_Merge[0]} >>$EER_Merge_Whole
#	grep -v 'Strand' -h ${EER_Merge[@]} >>$EER_Merge_Whole
#	#echo "#Random method: assigns the multi-mapped reads randomly" >$JER_Merge_Whole
#	grep 'Strand' ${JER_Merge[0]} >>$JER_Merge_Whole
#	grep -v 'Strand' -h ${JER_Merge[@]} >>$JER_Merge_Whole
#	
#	if $is_cleanup; then
#		rm ${GER[@]} ${EER[@]} ${JER[@]} -f
#		rm ${GER_Merge[@]} ${EER_Merge[@]} ${JER_Merge[@]} -f
#	fi
#	echo "Expression Level Estimation (Random Method) Done!"
#fi
#############################################################Done####################################################################################################
if $is_cleanup; then
        rm ${GEU[@]} -f
        rm ${EEU[@]} -f
        rm ${JEU[@]} -f
fi
if $is_cleanup; then
	rm ${TS[@]} -f
        rm ${USAM[@]} ${MUL[@]} ${MAP_LOG} -f
	rm ${EI[@]} ${EC[@]} ${JI[@]} ${JC[@]} -f
fi
#rm ${GEU[@]} ${EEU[@]} ${JEU[@]} -f
echo "Expression Level Estimation Done!"
#####################################################################################################################################################################
