#!/bin/bash
USAGE=$'Usage:

QC_SNP.sh {-f <r1.fastq.gz> or -1 <r1.fastq.gz> -2 <r2.fastq.gz>} [options]

Required arguments:

-f/--fastq <reads.fastq.gz>	reads fastq or fastq.gz file for single end data
-- OR --
-1/--read1 <reads.fastq.gz>	the first read fastq file for paired end data
-2/--read2 <reads.fastq.gz>	the second read fastq file for paired end data
-g/--genome <genome_ref.fa>     genome reference sequences
-c/--transcriptome <transcriptome_ref.fa>	transcriptome reference sequences
-a/--annotation <anno.gtf>	reference annotation in GTF format
-o/--output-prefix <outputPrefix>	prefix of output files

Optional arguments:

--mp [ MAX,MIN ]                both integers. bowtie2 option to set the mistch penalty, default is 6,2
--score-min [ Function,a,b ]    bowtie2 option to set a function of read length for the minimum alignment score
                                necessary for an alignment to be considered valid.
                                Default is L,0,-0.6 which is defined as Linear function:
                                f(x) = 0 + -0.6 * x  where x is the read length.
                                Available function types are constant (C), linear (L), square-root (S), and
                                natural log (G). The parameters are specified as F,B,A - the function type,
                                the constant term, and the coefficient separated by commas with no whitespace.
                                The constant term and coefficient may be negative and/or floating-point numbers.
                                For more info see the Bowtie2 manual.

-p/--pre-QC			generate QC reports before alignment. *Note: this option is ignored if -q/--QC is set
-q/--QC				generate QC reports after alignment
-s/--SNP			SNP calling
--ribo [ ribo_ref.fa ]		ribosomal RNA reference sequences
--mito [ mito_ref.fa ]		mitochondrial reference sequences
--gSam [ aln.sam ]		alignments to genome in SAM format
--tSam [ aln.sam ]		alignments to transcriptome in SAM format
--cleanup			delete temporary files
-h/--help			print this usage message'


echo ""
echo "You are running: $VERSION"
echo ""
if [ $# -eq 0 ]; then
        echo "No arguments or options!"
	echo "$USAGE"
        exit 1
fi

is_preQC=false
is_QC=false
is_SNP=false
is_cleanup=false
MP='6,2'
Score_Min='L,0,-0.6'
OUTPUT="QC_SNP_output"

declare -a ARGS

ARGS=($@)

for ((i=0;i<$#;++i))
do
        if [[ ${ARGS[i]} = '-f' || ${ARGS[i]} = '--fastq' ]]; then # Reads Input Fastq file
                FASTQ=${ARGS[(i+1)]}
                i=($i+1)
	elif [[ ${ARGS[i]} = '-1' || ${ARGS[i]} = '--read1' ]]; then # Reads Input Fastq file
                READ1=${ARGS[(i+1)]}
                i=($i+1)
        elif [[ ${ARGS[i]} = '-2' || ${ARGS[i]} = '--read2' ]]; then # Reads Input Fastq file
                READ2=${ARGS[(i+1)]}
                i=($i+1)
        elif [[ ${ARGS[i]} = '-g' || ${ARGS[i]} = '--genome' ]]; then # Genome Reference Sequences
                GENOME=${ARGS[(i+1)]}
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
	elif [[ ${ARGS[i]} = '-p' || ${ARGS[i]} = '--pre-QC' ]]; then # Choose pre-QC
         	is_preQC=true
	elif [[ ${ARGS[i]} = '-q' || ${ARGS[i]} = '--QC' ]]; then # Choose QC
                is_QC=true
	elif [[ ${ARGS[i]} = '-s' || ${ARGS[i]} = '--SNP' ]]; then # Choose SNP Calling
                is_SNP=true
	elif [  ${ARGS[i]} = '--ribo'  ]; then # ribosomal RNA
                RIBOSOME=${ARGS[(i+1)]}
                i=($i+1)
	elif [  ${ARGS[i]} = '--mito'  ]; then # mitochondrial
                MITOCHONDRIAL=${ARGS[(i+1)]}
                i=($i+1)
	elif [  ${ARGS[i]} = '--gSam'  ]; then # genome sam
                GenoSam=${ARGS[(i+1)]}
                i=($i+1)
	elif [  ${ARGS[i]} = '--tSam'  ]; then # transcriptome sam
                TranSam=${ARGS[(i+1)]}
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
if $is_preQC; then
	if [[ -z $FASTQ && -z $READ1 && -z $READ2 ]]; then
		echo "Error: You chose to do pre-QC, but there is no reads file!"
                echo "If you want to run with single end data, try option: -f < reads.fastq.gz >"
                echo "If you want to run with paired end data, try options: -1 < read1.fastq.gz > -2 < read2.fastq.gz >"
		exit 1
	fi

	if [[ -n $FASTQ && -n $READ1 ]] || [[ -n $FASTQ && -n $READ2 ]]; then
                echo "Error: single end data and paired end data can't be given together!"
                echo "If you want to run with single end data, try option: -f < reads.fastq.gz >"
                echo "If you want to run with paired end data, try options: -1 < read1.fastq.gz > -2 < read2.fastq.gz >"
                exit 1
        fi

        if [ -z $FASTQ ] && [[ -z $READ1 || -z $READ2 ]]; then
                echo "Error: for paired end data, read1 and read2 must be given together!"
                echo "If you want to run with paired end data, try options: -1 < read1.fastq.gz > -2 < read2.fastq.gz >"
                exit 1
        fi
fi
if $is_QC; then
	if [[ -z $GenoSam ]]; then
		if [[ -z $FASTQ && -z $READ1 && -z $READ2 || -z $GENOME || -z $TRANSCRIPTOME || -z $ANNOTATION ]]; then
			echo "Error: required input files not specified!"
			echo "$USAGE"
                	exit 1
		fi
	else
		if  [[ -z $FASTQ && -z $READ1 && -z $READ2 || -z $ANNOTATION ]]; then
			echo "Error: required input files not specified!"
			echo "$USAGE"
                	exit 1
		fi
	fi
	if [[ -n $FASTQ && -n $READ1 ]] || [[ -n $FASTQ && -n $READ2 ]]; then
       		echo "Error: single end data and paired end data can't be given together!"
        	echo "If you want to run with single end data, try option: -f < reads file >"
        	echo "If you want to run with paired end data, try options: -1 < read1 file > -2 < read2 file >"
        	exit 1
	fi

	if [ -z $FASTQ ] && [[ -z $READ1 || -z $READ2 ]]; then
        	echo "Error: for paired end data, read1 and read2 must be given together!"
                echo "If you want to run with paired end data, try options: -1 < read1.fastq.gz > -2 < read2.fastq.gz >"
        	exit 1
	fi

fi
if $is_SNP; then
	if [[ -z $GenoSam ]]; then
		if [[ -z $FASTQ && -z $READ1 && -z $READ2 || -z $GENOME || -z $TRANSCRIPTOME ]]; then
			echo "Error: required input files not specified!"
			echo "$USAGE"
			exit 1
		fi
		if [[ -n $FASTQ && -n $READ1 ]] || [[ -n $FASTQ && -n $READ2 ]]; then
  			echo "Error: single end data and paired end data can't be given together!"
			echo "If you want to run with single end data, try option: -f < reads.fastq.gz >"
			echo "If you want to run with paired end data, try options: -1 < read1.fastq.gz > -2 < read2.fastq.gz >"

        		exit 1
		fi

		if [ -z $FASTQ ] && [[ -z $READ1 || -z $READ2 ]]; then
        		echo "Error: for paired end data, read1 and read2 must be given together!"
                	echo "If you want to run with paired end data, try options: -1 < read1.fastq.gz > -2 < read2.fastq.gz >"
        		exit 1
		fi
	else
		if [[ -z $GENOME ]]; then
                        echo "Error: required Genome reference input file not specified!"
                        echo "$USAGE"
                        exit 1
		fi
	fi
fi

##########################################Check Analysis Options##########################################
is_only_alignment=false

if ! ($is_preQC || $is_QC || $is_SNP ); then
	echo "Warning: you didn't choose any analysis to do, so only alignment will be done."
	echo "If you wish, you may terminate the job and try re-running with one of the following options:"
	echo "-p/--pre-QC	 Description: generate QC reports before alignment. *Note: if --QC is set, --pre-QC will be ignored"
        echo "-q/--QC		Description: generate QC reports before and after alignment."
        echo "-s/--SNP		Description: SNP Calling "
        is_only_alignment=true
fi

#################################################################Check input files#################################################################
echo "Checking input files..."

if [[ -n $GENOME && -n $TRANSCRIPTOME && -n $ANNOTATION ]]; then
	echo "Check_for_reference_annotation_withinGenome.py -g $GENOME -t $TRANSCRIPTOME -a $ANNOTATION"
	Check_for_reference_annotation_withinGenome.py -g $GENOME -t $TRANSCRIPTOME -a $ANNOTATION
	ERR=$?
	if [ $ERR -ne 0 ]; then
        	echo "Error: input files failed to pass the check. There are some errors in the input files, please check them!"
        	exit 1
	fi
fi
if [[ -n $FASTQ ]]; then
	echo "Check_for_reads_file.py -r $FASTQ"
	Check_for_reads_file.py -r $FASTQ
	ERR=$?
	if [ $ERR -ne 0 ]; then
        	echo "Error: input files failed to pass the check. There are some errors in the input files, please check them!"
        	exit 1
	fi
else
	echo "Check_for_reads_file.py -r $READ1"
        Check_for_reads_file.py -r $READ1
	ERR=$?
	if [ $ERR -ne 0 ]; then
        	echo "Error: input files failed to pass the check. There are some errors in the input files, please check them!"
        	exit 1
	fi
        echo " Check_for_reads_file.py -r $READ2"
        Check_for_reads_file.py -r $READ2
	ERR=$?
	if [ $ERR -ne 0 ]; then
        	echo "Error: input files failed to pass the check. There are some errors in the input files, please check them!"
        	exit 1
	fi

fi
if [[ -n $GenoSam ]]; then
        echo "Check_for_genomeSam.py -s $GenoSam"
        Check_for_genomeSam.py -s $GenoSam
        ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Error: input files failed to pass the check. There are some errors in the input files, please check them!"
                exit 1
        fi
fi
if [[ -n $TranSam ]]; then
        echo "Check_for_transcriptomeSam.py -s $TranSam"
        Check_for_transcriptomeSam.py -s $TranSam
        ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Error: input files failed to pass the check. There are some errors in the input files, please check them!"
                exit 1
        fi
fi

#####################################################################Alignment#####################################################################
if ($is_QC) || ($is_SNP); then
    if [[ -z $GenoSam ]]; then
	#---------------------------------------------------------Alignment to Genome------------------------------------------------------------#
	GenoId=${GENOME%.fa} #get full path of filename without extension
	GenoIndexSamp=${GenoId##*/} #get basename of file
	GenoSam="$OUTPUT"_Bowtie2_genome.sam

        # Check if bowtie2 indexes already exist in reference path, if not run bowtie2-build
        # Can keep reference index so it doesn't need to be created every time if mapping multiple samples.

        echo "Checking for existing bowtie2 indexes for $GENOME"
	if [[ -f $GenoId.1.bt2 && -f $GenoId.2.bt2 && -f $GenoId.3.bt2 && -f $GenoId.4.bt2 && -f $GenoId.rev.1.bt2 && -f $GenoId.rev.2.bt2 ]]
        then
            echo "Existing bowtie2 indexes found for $GENOME"
            GenoIndex=$GenoId
        else
            echo "Start building bowtie2 indexes for $GENOME"
            GenoIndex=$GenoIndexSamp
	    echo "bowtie2-build $GENOME $GenoIndex"
            bowtie2-build $GENOME $GenoIndex
            ERR=$?
            if [ $ERR -ne 0 ]; then
                echo "Errors when building bowtie2 index for $GENOME, stopping the pipeline!"
                if [ -f $GenoIndexSamp ]; then
                    rm $GenoIndexSamp* -f
                fi
                exit 1
            fi
        fi

        echo "Starting genome alignment with bowtie2"
	# Output: $GenoSam  ("$OUTPUT"_Bowtie2_genome.sam)
	if [[ -n $FASTQ ]]; then
		echo "bowtie2 -x $GenoIndex -U $FASTQ --mp $MP --score-min $Score_Min --sam-no-hd -S $GenoSam"
	        bowtie2 -x $GenoIndex -U $FASTQ --mp $MP --score-min $Score_Min --sam-no-hd -S $GenoSam
	else
		echo "bowtie2 -x $GenoIndex -1 $READ1 -2 $READ2 --mp $MP --score-min $Score_Min --sam-no-hd -S $GenoSam"
		bowtie2 -x $GenoIndex -1 $READ1 -2 $READ2 --mp $MP --score-min $Score_Min --sam-no-hd -S $GenoSam
	fi
	ERR=$?
	if [ $ERR -ne 0 ]; then
              	echo "Errors when running Bowtie2 (mapping to genome), stopping the pipeline!"
              	exit 1
   	fi
	if [[ $is_cleanup && -f $GenoIndexSamp ]]; then
                rm $GenoIndexSamp* -f
        fi

    fi

    if [[ -z $TranSam ]]; then
        #-----------------------------------------------------Alignment to Transcriptome---------------------------------------------------------#
        TranId=${TRANSCRIPTOME%.fa} #get full path of filename without extension
        TrantoIndexSamp=${TranId##*/} #get basename of file
        TranSam="$OUTPUT"_Bowtie2_transcriptome.sam

	# Check if bowtie2 indexes already exist in reference path, if not run bowtie2-build
	# Can keep reference index so it doesn't need to be created every time if mapping multiple samples.

        echo "Checking for existing bowtie2 indexes for $TRANSCRIPTOME"
        if [[ -f $TranId.1.bt2 && -f $TranId.2.bt2 && -f $TranId.3.bt2 && -f $TranId.4.bt2 && -f $TranId.rev.1.bt2 && -f $TranId.rev.2.bt2 ]]
        then
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

	echo "Starting transcriptome alignment with bowtie2"
	# Output: $TranSam ("$OUTPUT"_Bowtie2_transcriptome.sam)
	if [[ -n $FASTQ ]]; then
		echo "bowtie2 -x $TrantoIndex -U $FASTQ --mp $MP --score-min $Score_Min --sam-no-hd -S $TranSam"
	        bowtie2 -x $TrantoIndex -U $FASTQ --mp $MP --score-min $Score_Min --sam-no-hd -S $TranSam
	else
		echo "bowtie2 -x $TrantoIndex -1 $READ1 -2 $READ2 --mp $MP --score-min $Score_Min --sam-no-hd -S $TranSam"
		bowtie2 -x $TrantoIndex -1 $READ1 -2 $READ2 --mp $MP --score-min $Score_Min --sam-no-hd -S $TranSam
	fi
	ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running Bowtie2 (mapping to transcriptome), stopping the pipeline!"
                exit 1
        fi
	if [[ $is_cleanup && -f $TrantoIndexSamp ]]; then
                rm $TrantoIndexSamp* -f
        fi
    fi
	#-------------------------------------------------------Alignment to rRNA-----------------------------------------------------------#
    if ($is_QC); then
       	if [[ -n $RIBOSOME ]]; then
		RiboId=${RIBOSOME%.fa} #get full path of filename without extension
		RiboIndexSamp=${RiboId##*/} #get basename of file
           	RiboSam="$OUTPUT"_Bowtie2_ribosome.sam

	        # Check if bowtie2 indexes already exist in reference path, if not run bowtie2-build
		# Can keep reference index so it doesn't need to be created every time if mapping multiple samples.

	        echo "Checking for existing bowtie2 indexes for $RIBOSOME"
		if [[ -f $RiboId.1.bt2 && -f $RiboId.2.bt2 && -f $RiboId.3.bt2 && -f $RiboId.4.bt2 && -f $RiboId.rev.1.bt2 && -f $RiboId.rev.2.bt2 ]]
		then
		    echo "Existing bowtie2 indexes found for $RIBOSOME"
		    RiboIndex=$RiboId
		else
		    echo "Start building bowtie2 indexes for $RIBOSOME"
		    RiboIndex=$RiboIndexSamp
	            echo "bowtie2-build $RIBOSOME $RiboIndex"
		    bowtie2-build $RIBOSOME $RiboIndex
		    ERR=$?
		    if [ $ERR -ne 0 ]; then
			echo "Errors when building bowtie2 index for $RIBOSOME, stopping the pipeline!"
			if [ -f $RiboIndexSamp ]; then
			    rm $RiboIndexSamp* -f
			fi
			exit 1
		    fi
		fi

		echo "Starting ribosome alignment with bowtie2"
		# Output: $RiboSam ("$OUTPUT"_Bowtie2_ribosome.sam)
       	        if [[ -n $FASTQ ]]; then
			echo "bowtie2 -x $RiboIndex -U $FASTQ --mp $MP --score-min $Score_Min --sam-no-hd -S $RiboSam"
                       	bowtie2 -x $RiboIndex -U $FASTQ --mp $MP --score-min $Score_Min --sam-no-hd -S $RiboSam
               	else
			echo "bowtie2 -x $RiboIndex -1 $READ1 -2 $READ2 --mp $MP --score-min $Score_Min --sam-no-hd -S $RiboSam"
                  	bowtie2 -x $RiboIndex -1 $READ1 -2 $READ2 --mp $MP --score-min $Score_Min --sam-no-hd -S $RiboSam
               	fi
		ERR=$?
		if [ $ERR -ne 0 ]; then
               		echo "Errors when running Bowtie2 (mapping to ribosome), stopping the pipeline!"
	             	exit 1
        	fi
		if [[ $is_cleanup && -f $RiboIndexSamp ]]; then
       		        rm $RiboIndexSamp* -f
       		fi
        fi

        #-------------------------------------------------------Alignment to chrM-----------------------------------------------------------#
	if [[ -n $MITOCHONDRIAL ]]; then
		MitoId=${MITOCHONDRIAL%.fa} #get full path of filename without extension
	        MitoIndexSamp=${MitoId##*/} #get basename of file
		MitoSam="$OUTPUT"_Bowtie2_mitochondrial.sam

        	# Check if bowtie2 indexes already exist in reference path, if not run bowtie2-build
        	# Can keep reference index so it doesn't need to be created every time if mapping multiple samples.

        	echo "Checking for existing bowtie2 indexes for $MITOCHONDRIAL"
	        MitoId=${MITOCHONDRIAL%.fa} #get filename without extension
        	if [[ -f $MitoId.1.bt2 && -f $MitoId.2.bt2 && -f $MitoId.3.bt2 && -f $MitoId.4.bt2 && -f $MitoId.rev.1.bt2 && -f $MitoId.rev.2.bt2 ]]
        	then
            	    echo "Existing bowtie2 indexes found for $MITOCHONDRIAL"
            	    MitoIndex=$MitoId
        	else
            	    echo "Start building bowtie2 indexes for $MITOCHONDRIAL"
                    MitoIndex=$MitoIndexSamp
		    echo "bowtie2-build $MITOCHONDRIAL $MitoIndex"
                    bowtie2-build $MITOCHONDRIAL $MitoIndex
                    ERR=$?
                    if [ $ERR -ne 0 ]; then
                        echo "Errors when building bowtie2 index, stopping the pipeline!"
                        if [ -f $MitoIndexSamp ]; then
                            rm $MitoIndexSamp* -f
                        fi
                        exit 1
                   fi
                fi

                echo "Starting mitochondrial alignment with bowtie2"
		# Output: $MitoSam ("$OUTPUT"_Bowtie2_mitochondrial.sam)
		if [[ -n $FASTQ ]]; then
			echo "bowtie2 -x $MitoIndex -U $FASTQ --mp $MP --score-min $Score_Min --sam-no-hd -S $MitoSam"
			bowtie2 -x $MitoIndex -U $FASTQ --mp $MP --score-min $Score_Min --sam-no-hd -S $MitoSam
		else
			echo "bowtie2 -x $MitoIndex -1 $READ1 -2 $READ2 --mp $MP --score-min $Score_Min --sam-no-hd -S $MitoSam"
			bowtie2 -x $MitoIndex -1 $READ1 -2 $READ2 --mp $MP --score-min $Score_Min --sam-no-hd -S $MitoSam
		fi
		ERR=$?
                if [ $ERR -ne 0 ]; then
                        echo "Errors when running Bowtie2 (mapping to mitochondrial), stopping the pipeline!"
                	exit 1
		fi
		if [[ $is_cleanup && -f $MitoIndexSamp ]]; then
                	rm $MitoIndexSamp* -f
        	fi
	fi
     fi
fi

if $is_only_alignment; then
        exit 0
fi

###############################################################Merge Mapping Result################################################################
#MergeSam="$OUTPUT"_merge.sam
UniqueSam="$OUTPUT"_Unique.sam
MergeSam=$UniqueSam
if [[ -z $TranSam && -n $GenoSam ]]; then
# JSH 2013-08-27
#        cat $GenoSam >$MergeSam
#	cat $GenoSam >"$OUTPUT"_Unique.sam
	ln -s $GenoSam $UniqueSam
elif [[ -n $TranSam && -n $GenoSam ]]; then
	# Output: "$OUTPUT"_chrList_gtf.txt
	echo "SplitByChromosome_for_annotationGtfFile.py -i $ANNOTATION -p "$OUTPUT" -o "$OUTPUT"_chrList_gtf.txt"
	SplitByChromosome_for_annotationGtfFile.py -i $ANNOTATION -p "$OUTPUT" -o "$OUTPUT"_chrList_gtf.txt
	ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running SplitByChromosome_for_annotationGtfFile.py, stopping the pipeline!"
                exit 1
        fi
	chrList=`cat "$OUTPUT"_chrList_gtf.txt`
	for i in $chrList
	do
		# Output: "$OUTPUT"_JunctionCombination_"$i".txt
		echo "JunctionCombination.py -g "$OUTPUT"_"$i"_annotation.gtf -o "$OUTPUT"_JunctionCombination_"$i".txt"
		JunctionCombination.py -g "$OUTPUT"_"$i"_annotation.gtf -o "$OUTPUT"_JunctionCombination_"$i".txt
		ERR=$?
        	if [ $ERR -ne 0 ]; then
                	echo "Errors when running JunctionCombination.py, stopping the pipeline!"
                	exit 1
        	fi

		# Output: "$OUTPUT"_JunctionIndex_"$i".txt
                echo "JunctionIndex.py -g "$OUTPUT"_"$i"_annotation.gtf -j "$OUTPUT"_JunctionCombination_"$i".txt -o "$OUTPUT"_JunctionIndex_"$i".txt"
		JunctionIndex.py -g "$OUTPUT"_"$i"_annotation.gtf -j "$OUTPUT"_JunctionCombination_"$i".txt -o "$OUTPUT"_JunctionIndex_"$i".txt
		ERR=$?
		if [ $ERR -ne 0 ]; then
                        echo "Errors when running JunctionIndex.py, stopping the pipeline!"
                        exit 1
                fi
	done
	cat "$OUTPUT"_JunctionIndex_chr*.txt >"$OUTPUT"_JunctionIndex.txt
	if [[ -z $FASTQ ]]; then
		# Output: "$OUTPUT"_Unique.sam and "$OUTPUT"_Multiple.sam
		echo "Merge_genome_transcriptome_pairEnd "$OUTPUT"_JunctionIndex.txt $TranSam $GenoSam "$OUTPUT"_Unique.sam "$OUTPUT"_Multiple.sam"
		Merge_genome_transcriptome_pairEnd "$OUTPUT"_JunctionIndex.txt $TranSam $GenoSam "$OUTPUT"_Unique.sam "$OUTPUT"_Multiple.sam
		ERR=$?
                if [ $ERR -ne 0 ]; then
                 	echo "Errors when running Merge_genome_transcriptome_pairEnd, stopping the pipeline!"
                      	exit 1
              	fi
	else
		# Output: "$OUTPUT"_Unique.sam and "$OUTPUT"_Multiple.sam
		echo "Merge_genome_transcriptome_singleEnd "$OUTPUT"_JunctionIndex.txt $TranSam $GenoSam "$OUTPUT"_Unique.sam "$OUTPUT"_Multiple.sam"
		Merge_genome_transcriptome_singleEnd "$OUTPUT"_JunctionIndex.txt $TranSam $GenoSam "$OUTPUT"_Unique.sam "$OUTPUT"_Multiple.sam
		ERR=$?
                if [ $ERR -ne 0 ]; then
                        echo "Errors when running Merge_genome_transcriptome_singleEnd, stopping the pipeline!"
                        exit 1
                fi
	fi
#JSH 2013-08-27
#	cat "$OUTPUT"_Unique.sam "$OUTPUT"_Multiple.sam >$MergeSam
	if $is_cleanup; then
		rm "$OUTPUT"_chr*_annotation.gtf "$OUTPUT"_chrList_gtf.txt
		rm "$OUTPUT"_JunctionCombination_chr*.txt -f
		rm "$OUTPUT"_JunctionIndex_chr*.txt -f
		rm "$OUTPUT"_JunctionIndex.txt  "$OUTPUT"_Multiple.sam -f
	fi
fi
######################################################################Pre-QC#######################################################################
if [[ -z $FASTQ ]]; then
	if [[ $READ1 =~ fastq.gz$ ]]; then
	    FASTQ="$OUTPUT"_read_1_2.fastq.gz
	else
            FASTQ="$OUTPUT"_read_1_2.fastq
	fi
        cat $READ1 $READ2 >$FASTQ
fi
if $is_preQC; then
	#--------------------Read Quality-------------------#
	echo "Calculate reads quality"
        # Output: "$OUTPUT".read_qual.r
	read_quality.py -i $FASTQ -o $OUTPUT
	ERR=$?
        if [ $ERR -ne 0 ]; then
               	echo "Errors when running read_quality.py, stopping the pipeline!"
                exit 1
        fi
	#------------------------ATCG-----------------------#
	#echo "Calculate squence(ATCG & GC) content across all bases"
	#fastx_quality_stats -i $FASTQ -o "$OUTPUT".reads_stat.txt
	#read_stat_plot.py -i "$OUTPUT".reads_stat.txt
	#---------------------GC Content--------------------#
	echo "Calculate squence(ATCG & GC) content across read length"
	echo "Calculate GC distribution over all sequences"
	# Output: "$OUTPUT".ATCG_GC_Percentage_plot.r, "$OUTPUT".GC_hist_plot.r, "$OUTPUT".GC_hist.xls, "$OUTPUT".ACGT_Percentage.txt, "$OUTPUT".GC_Percentage.txt
	read_ACGT_GC.py -i $FASTQ -o $OUTPUT
	ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running read_ACGT_GC.py, stopping the pipeline!"
                exit 1
        fi
	if $is_cleanup; then
		rm $OUTPUT.read_qual.r -f
		rm $OUTPUT.ATCG_GC_Percentage_plot.r -f
		rm $OUTPUT.GC_hist_plot.r -f
	fi
	echo "Pre-QC Done!"
fi
########################################################################QC#########################################################################
if $is_QC; then
	#----------------------Convert----------------------#
        echo "convert gtf format to bed format"
        grep -w -v 'gene' $ANNOTATION >"$OUTPUT"_annotation_without_gene.gtf
	# Output: "$OUTPUT"_Table_convert.txt
	echo "gtfToGenePred "$OUTPUT"_annotation_without_gene.gtf "$OUTPUT"_Table_convert.txt"
        gtfToGenePred "$OUTPUT"_annotation_without_gene.gtf "$OUTPUT"_Table_convert.txt
	ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running gtfToGenePred, stopping the pipeline!"
                exit 1
        fi

	# Output: "$OUTPUT"_Bed_convert.bed
	echo "genePredToBed "$OUTPUT"_Table_convert.txt >"$OUTPUT"_Bed_convert.bed"
        genePredToBed "$OUTPUT"_Table_convert.txt >"$OUTPUT"_Bed_convert.bed
	ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running genePredToBed, stopping the pipeline!"
                exit 1
        fi
        #--------------------Read Quality-------------------#
        #echo "Calculate reads quality"
        #read_quality.py -i $FASTQ -o $OUTPUT
        #------------------------ATCG-----------------------#
        #echo "Calculate squence(ATCG & GC) content across all bases"
        #fastx_quality_stats -i $FASTQ -o "$OUTPUT".reads_stat.txt
        #read_stat_plot.py -i "$OUTPUT".reads_stat.txt
        #---------------------GC Content--------------------#
        #echo "Calculate GC distribution over all sequences"
        #read_GC.py -i $FASTQ -o $OUTPUT
	#-----------------mapping statistics----------------#
	echo "Calculate mapping statistics"
	# Output: "$OUTPUT".mapping_report.txt
	echo "mapping_stat.py -s $MergeSam -q $FASTQ -o $OUTPUT"
	mapping_stat.py -s $MergeSam -q $FASTQ -o $OUTPUT
	ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running mapping_stat.py, stopping the pipeline!"
                exit 1
        fi
        echo "converting $OUTPUT.mapping_report.txt to pdf file"
        echo "text2pdf $OUTPUT.mapping_report.txt > $OUTPUT.mapping_report.pdf"
        text2pdf $OUTPUT.mapping_report.txt > $OUTPUT.mapping_report.pdf
    
	#-----------------strand specificity----------------#
	echo "Calculate reads strand specificity"
	# Output: "$OUTPUT".strand_stat.txt
	echo "strand_specificity.py -s $MergeSam  -a "$OUTPUT"_Bed_convert.bed -o $OUTPUT"
	strand_specificity.py -s $MergeSam  -a "$OUTPUT"_Bed_convert.bed -o $OUTPUT
	ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running strand_specificity.py, stopping the pipeline!"
                exit 1
        fi
        echo "converting $OUTPUT.strand_stat.txt to pdf file"
        echo "text2pdf $OUTPUT.strand_stat.txt > $OUTPUT.strand_stat.pdf"
        text2pdf $OUTPUT.strand_stat.txt > $OUTPUT.strand_stat.pdf


	#---------------Gene Coverage Analysis--------------#
	echo "Get genes with single transcript for gene body coverage analysis"
	# Output: "$OUTPUT"_singleGene.gtf
	echo "single_gene_forGtf.py -g "$OUTPUT"_annotation_without_gene.gtf -o "$OUTPUT"_singleGene.gtf"
	single_gene_forGtf.py -g "$OUTPUT"_annotation_without_gene.gtf -o "$OUTPUT"_singleGene.gtf
	ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running single_gene_forGtf.py, stopping the pipeline!"
                exit 1
        fi
	# Output: "$OUTPUT"_singleGene_Table_convert.txt
	echo "gtfToGenePred "$OUTPUT"_singleGene.gtf "$OUTPUT"_singleGene_Table_convert.txt"
	gtfToGenePred "$OUTPUT"_singleGene.gtf "$OUTPUT"_singleGene_Table_convert.txt
	ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running gtfToGenePred, stopping the pipeline!"
                exit 1
        fi
	# Output: "$OUTPUT"_singleGene_Bed_convert.bed
	echo "genePredToBed "$OUTPUT"_singleGene_Table_convert.txt >"$OUTPUT"_singleGene_Bed_convert.bed"
	genePredToBed "$OUTPUT"_singleGene_Table_convert.txt >"$OUTPUT"_singleGene_Bed_convert.bed
	ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running genePredToBed, stopping the pipeline!"
                exit 1
        fi
	# Output: "$OUTPUT"_single_transcriptGene.geneBodyCoverage_plot.r, "$OUTPUT"_single_transcriptGene.geneBody_coverage.pdf, "$OUTPUT"_single_transcriptGene.geneBodyCoverage.txt
	echo "geneBody_coverage.py -s $MergeSam -a "$OUTPUT"_singleGene_Bed_convert.bed -o "$OUTPUT"_single_transcriptGene"
	geneBody_coverage.py -s $MergeSam -a "$OUTPUT"_singleGene_Bed_convert.bed -o "$OUTPUT"_single_transcriptGene
	ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running geneBody_coverage.py(for the genes containing single transcript), stopping the pipeline!"
                exit 1
        fi
	echo "Calculate gene coverage for whole genes"
	# Output: "$OUTPUT"_whole_genes.geneBodyCoverage_plot.r, "$OUTPUT"_whole_genes.geneBody_coverage.pdf, "$OUTPUT"_whole_genes.geneBodyCoverage.txt
	echo "geneBody_coverage.py -s $MergeSam -a "$OUTPUT"_Bed_convert.bed -o "$OUTPUT"_whole_genes"
	geneBody_coverage.py -s $MergeSam -a "$OUTPUT"_Bed_convert.bed -o "$OUTPUT"_whole_genes
	ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running geneBody_coverage.py(for whole genes), stopping the pipeline!"
                exit 1
        fi
	#-----------------Reads Distribution----------------#
	echo "Calculate reads distribution"
	# Output: "$OUTPUT".utr_3.txt, "$OUTPUT".utr_5.txt, "$OUTPUT".cds_exon.txt, "$OUTPUT".intron.txt
	echo "Get_UTR_Exon_Intron.py -a "$OUTPUT"_Bed_convert.bed -o $OUTPUT"
	Get_UTR_Exon_Intron.py -a "$OUTPUT"_Bed_convert.bed -o $OUTPUT #get UTR exon intron region
        ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors come out when running Get_UTR_Exon_Intron.py, stopping the pipeline!"
                exit 1
        fi
	# Output: "$OUTPUT".intergenic_up_1kb.txt and "$OUTPUT".intergenic_down_1kb.txt
	echo "Get_Intergenic_1kb.py -a "$OUTPUT"_Bed_convert.bed -p $OUTPUT -o $OUTPUT"
        Get_Intergenic_1kb.py -a "$OUTPUT"_Bed_convert.bed -p $OUTPUT -o $OUTPUT #get Intergenic upstream and downstream 1kb regopm
        ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors come out when running Get_Intergenic_1kb.py, stopping the pipeline!"
                exit 1
        fi
	# Output: "$OUTPUT".intergenic_up_5kb.txt and "$OUTPUT".intergenic_down_5kb.txt
	echo "Get_Intergenic_5kb.py -a "$OUTPUT"_Bed_convert.bed -p $OUTPUT -o $OUTPUT"
        Get_Intergenic_5kb.py -a "$OUTPUT"_Bed_convert.bed -p $OUTPUT -o $OUTPUT #get Intergenic upstream and downstream 5kb regop
        ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors come out when running Get_Intergenic_5kb.py, stopping the pipeline!"
                exit 1
        fi
	# Output: "$OUTPUT".intergenic_up_10kb.txt and "$OUTPUT".intergenic_down_10kb.txt
        echo "Get_Intergenic_10kb.py -a "$OUTPUT"_Bed_convert.bed -p $OUTPUT -o $OUTPUT"
        Get_Intergenic_10kb.py -a "$OUTPUT"_Bed_convert.bed -p $OUTPUT -o $OUTPUT #get Intergenic upstream and downstream 10kb regop
        ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors come out when running Get_Intergenic_10kb.py, stopping the pipeline!"
                exit 1
        fi

	#  Internal Input Files required for read_distribution.py : "$OUTPUT".utr_3.txt, "$OUTPUT".utr_5.txt, "$OUTPUT".cds_exon.txt, "$OUTPUT".intron.txt, "$OUTPUT".intergenic_up_1kb.txt, "$OUTPUT".intergenic_up_5kb.txt, "$OUTPUT".intergenic_up_10kb.txt, "$OUTPUT".intergenic_down_1kb.txt, "$OUTPUT".intergenic_down_5kb.txt, "$OUTPUT".intergenic_down_10kb.txt

	# Output: "$OUTPUT".read_distribution.txt
	echo "read_distribution.py -s $MergeSam -p $OUTPUT -o $OUTPUT"
        read_distribution.py -s $MergeSam -p $OUTPUT -o $OUTPUT
	ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running read_distribution.py, stopping the pipeline!"
                exit 1
        fi
        echo "converting $OUTPUT.read_distribution.txt to pdf file"
        echo "text2pdf $OUTPUT.read_distribution.txt > $OUTPUT.read_distribution.pdf"
        text2pdf $OUTPUT.read_distribution.txt > $OUTPUT.read_distribution.pdf

        if $is_cleanup; then
                rm "$OUTPUT".utr_3.txt -f
                rm "$OUTPUT".utr_5.txt -f
                rm "$OUTPUT".cds_exon.txt -f
                rm "$OUTPUT".intron.txt -f
                rm "$OUTPUT".intergenic_up_1kb.txt -f
                rm "$OUTPUT".intergenic_up_5kb.txt -f
                rm "$OUTPUT".intergenic_up_10kb.txt -f
                rm "$OUTPUT".intergenic_down_1kb.txt -f
                rm "$OUTPUT".intergenic_down_5kb.txt -f
                rm "$OUTPUT".intergenic_down_10kb.txt -f
        fi
	#----------------chrM and rRNA report---------------#
	if [[ -n $READ1 ]]; then
		FASTQ=$READ1
	fi
	if [[ -n $RIBOSOME ]]; then
		# Output: "$OUTPUT"_rRNA.simple_mapping_report.txt
		echo "simple_mapping_report.py -s $RiboSam -q $FASTQ -o "$OUTPUT"_rRNA"
		simple_mapping_report.py -s $RiboSam -q $FASTQ -o "$OUTPUT"_rRNA -t 'Ribosomal RNA mapping report'
		ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running simple_mapping_report.py(for ribosome), stopping the pipeline!"
                exit 1
        fi
        echo "converting $OUTPUT_rRNA.simple_mapping_report.txt to pdf file"
        echo "text2pdf $OUTPUT_rRNA.simple_mapping_report.txt > $OUTPUT_rRNA.simple_mapping_report.pdf"
        text2pdf $OUTPUT_rRNA.simple_mapping_report.txt > $OUTPUT_rRNA.simple_mapping_report.pdf

	fi
	if [[ -n $MITOCHONDRIAL ]]; then
		# Output: "$OUTPUT"_chrM.simple_mapping_report.txt
                echo "simple_mapping_report.py -s $MitoSam -q $FASTQ -o "$OUTPUT"_chrM"
		simple_mapping_report.py -s $MitoSam -q $FASTQ -o "$OUTPUT"_chrM  -t 'Mitochondrial genome mapping report'
		ERR=$?
		if [ $ERR -ne 0 ]; then
            echo "Errors when running simple_mapping_report.py(for mitochondrial), stopping the pipeline!"
            exit 1
        fi
        echo "converting $OUTPUT_chrM.simple_mapping_report.txt to pdf file"            echo "text2pdf $OUTPUT_chrM.simple_mapping_report.txt > $OUTPUT_chrM.simple_mapping_report.pdf"
        text2pdf $OUTPUT_chrM.simple_mapping_report.txt > $OUTPUT_chrM.simple_mapping_report.pdf
	fi

        # If "ImageMagick" installed, then combine all pdfs into a single pdf
        echo "combining QC pdfs into a single pdf using Image Magick"
        echo "identify -verison"
        imageMagick=`identify -version`
        if [[ $imageMagick =~ "Version: ImageMagick" ]]; then
            echo "found ImageMagick, combining QC pdfs into a single pdf: ${OUTPUT}_combined_QC_report.pdf"
            echo "convert *pdf ${OUTPUT}_combined_QC_report.pdf"
            convert ${OUTPUT}*pdf ${OUTPUT}_combined_QC_report.pdf
        else
            echo "ImageMagick not found. If you wish to combine pdfs, please download and install ImageMagick."
        fi
         
	if $is_cleanup; then
		rm "$OUTPUT"_annotation_without_gene.gtf "$OUTPUT"_Table_convert.txt "$OUTPUT"_Bed_convert.bed -f
		rm "$OUTPUT"_singleGene.gtf "$OUTPUT"_singleGene_Table_convert.txt "$OUTPUT"_singleGene_Bed_convert.bed -f
		rm "$OUTPUT"_whole_genes.geneBodyCoverage_plot.r "$OUTPUT"_single_transcriptGene.geneBodyCoverage_plot.r -f
	fi
	echo "QC Done!"
fi
####################################################################SNP calling####################################################################
if $is_SNP; then
	# Output: $GENOME.fai
	echo "samtools faidx $GENOME"
	samtools faidx $GENOME

	# Output: "$OUTPUT"_Unique.bam
	echo "samtools view -bt "$GENOME".fai "$OUTPUT"_Unique.sam > "$OUTPUT"_Unique.bam"
	samtools view -bt "$GENOME".fai "$OUTPUT"_Unique.sam > "$OUTPUT"_Unique.bam
	ERR=$?
        if [ $ERR -ne 0 ]; then
            echo "Errors when running samtools(convert sam to bam), stopping the pipeline!"
            exit 1
        fi

        # Output: "$OUTPUT"_Unique_sort.bam
        echo "samtools sort "$OUTPUT"_Unique.bam "$OUTPUT"_Unique_sort"
        samtools sort "$OUTPUT"_Unique.bam "$OUTPUT"_Unique_sort
        ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running samtools(sort the bam file), stopping the pipeline!"
                exit 1
        fi
        # Output: "$OUTPUT"_SNPsvar.raw.bcf
        echo "samtools mpileup -uf $GENOME "$OUTPUT"_Unique_sort.bam|bcftools view -bvcg - >"$OUTPUT"_SNPsvar.raw.bcf"
        samtools mpileup -uf $GENOME "$OUTPUT"_Unique_sort.bam|bcftools view -bvcg - >"$OUTPUT"_SNPsvar.raw.bcf
        ERR=$?
        if [ $ERR -ne 0 ]; then
                echo "Errors when running samtools and bcftools(do SNP calling), stopping the pipeline!"
                exit 1
        fi


	if $is_cleanup; then
		rm "$GENOME".fai -f
		rm "$OUTPUT"_Unique.sam "$OUTPUT"_Unique.bam -f
		rm "$OUTPUT"_Unique_sort.bam -f
	fi
	echo "SNP Calling Done!"
fi
########################################################################Done#######################################################################
echo "QC_SNP Part Done!"
