#!/bin/bash

USAGE=$'Usage:

FileFormat -i <inputFile> [-r <referenceFile>] -b|-m|-d|-w -o <outputPrefix>

-i/--input <inputFile>			input file that you want to convert, must have .sam, .bam, or .mrf extension
-o/--output <outputPrefix>		output prefix for the converted file. The converted format extension will be added automatically
-r/--reference <referenceFile>		If a file is being converted from SAM to BAM format, the sam file should contain either the appropriate @SQ headers or you need to use the --reference option to supply the reference file that was used to create the sam file 
-b/--toBAM				Convert SAM to BAM
-m/--toMRF				Convert SAM to MRF
-d/--toBED				Convert BAM to BED
-w/--toWIG				Convert SAM to WIG or MRF to WIG'


if [ $# -eq 0 ]; then
        echo "No options or arguments!"
        echo "$USAGE"
        exit 1
fi

declare -a ARGS

ARGS=($@)

OUTPUT="FileFormatConversion_output"

for ((i=0;i<$#;++i))
do
        if [[ ${ARGS[i]} = '-i' || ${ARGS[i]} = '--input' ]]; then # Input files
                INPUT=${ARGS[(i+1)]}
                i=($i+1)
        elif [[ ${ARGS[i]} = '-o' || ${ARGS[i]} = '--output' ]]; then # Prefix of output files
                OUTPUT=${ARGS[(i+1)]}
                i=($i+1)
	elif [[ ${ARGS[i]} = '-r' || ${ARGS[i]} = '--reference' ]]; then
		REFERENCE=${ARGS[(i+1)]}
                i=($i+1)
	elif [[ ${ARGS[i]} = '-b' || ${ARGS[i]} = '--toBAM' ]]; then
		if [[ -z $convertFormat ]]; then
			convertFormat='bam'
			#i=($i+1)
		else
			echo "Can't convert input file to multiple formats in the same command! Please choose only one format at a time for conversion."
			exit 1
		fi
	elif [[ ${ARGS[i]} = '-m' || ${ARGS[i]} = '--toMRF' ]]; then
		if [[ -z $convertFormat ]]; then
                        convertFormat='mrf'
                        #i=($i+1)
                else
                        echo "Can't convert input file to multiple formats in the same command! Please choose only one format at a time for conversion."
                        exit 1
                fi
	elif [[ ${ARGS[i]} = '-d' || ${ARGS[i]} = '--toBED' ]]; then
                if [[ -z $convertFormat ]]; then
                        convertFormat='bed'
                        #i=($i+1)
                else
                        echo "Can't convert input file to multiple formats in the same command! Please choose only one format at a time for conversion."
                        exit 1
                fi
	elif [[ ${ARGS[i]} = '-w' || ${ARGS[i]} = '--toWIG' ]]; then
                if [[ -z $convertFormat ]]; then
                        convertFormat='wig'
                        #i=($i+1)
                else
                        echo "Can't convert input file to multiple formats in the same command! Please choose only one format at a time for conversion."
                        exit 1
                fi
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
if [[ -z $INPUT ]]; then
	echo "Required input file for conversion is missing!"
	echo "Try option: -i < in.file >"
	exit 1
fi

FileExtension=`echo ${INPUT##*.}`
if [[ "$FileExtension"x = 'samx' ]]; then
	
	if [[ $convertFormat = 'bed' ]]; then
		echo "Conversion not provided for $FileExtension to $convertFormat!"
		echo "Current conversions: SAM to BAM, SAM to MRF, SAM to WIG, BAM to BED, MRF to WIG"
		exit 1
	fi
elif [[ "$FileExtension"x = 'bamx' ]]; then
	if [[ $convertFormat != 'bed' ]]; then
		echo "Conversion not provided for $FileExtension to $convertFormat!"
                echo "Current conversions: SAM to BAM, SAM to MRF, SAM to WIG, BAM to BED, MRF to WIG"
                exit 1
        fi
elif [[ "$FileExtension"x = 'mrfx' ]]; then
        if [[ $convertFormat != 'wig' ]]; then
                echo "Conversion not provided for $FileExtension to $convertFormat!"
                echo "Current conversions: SAM to BAM, SAM to MRF, SAM to WIG, BAM to BED, MRF to WIG."
                exit 1
        fi
else 
	echo "Input file can't be converted, only files in SAM, BAM and MRF format can be converted."
	exit 1
fi

if [[ $FileExtension = 'sam' ]]; then
	if [[ $convertFormat = 'bam' ]]; then
		if [[ -z $REFERENCE ]]; then
			echo "If a file is being converted from SAM to BAM format, the sam file should contain either the appropriate @SQ headers or you need to use the --reference option to supply the reference file that was used to create the sam file "
			samtools view -bS $INPUT > "$OUTPUT".bam
		else
			samtools faidx $REFERENCE
			samtools view -bt "$REFERENCE".fai $INPUT > "$OUTPUT".bam
			rm "$REFERENCE".fai -f
		fi
	elif [[ $convertFormat = 'mrf' ]]; then
		sort -r $INPUT | sam2mrf > "$OUTPUT".mrf 
	elif [[ $convertFormat = 'wig' ]]; then
		sort -r $INPUT | sam2mrf > "$OUTPUT".mrf
		mrf2wig $OUTPUT < "$OUTPUT".mrf
		rm "$OUTPUT".mrf -f
	fi
elif [[ $FileExtension = 'bam' ]]; then
	bamToBed -i $INPUT > "$OUTPUT".bed
elif [[ $FileExtension = 'mrf' ]]; then
	mrf2wig $OUTPUT < $INPUT
else
        echo "Input file can't be converted, only files in SAM, BAM and MRF format can be converted."
        exit 1
fi
