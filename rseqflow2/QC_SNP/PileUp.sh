#!/bin/bash

set -e

declare -a ARGS 

ARGS=($@)

# Command line arguments
for ((i=0;i<$#;++i))
do
	if [[ ${ARGS[i]} = '-i' || ${ARGS[i]} = '--bam-input' ]]; then # BAM Input file
		BAMINPUT=${ARGS[(i+1)]}
		i=($i+1)
	elif [[ ${ARGS[i]} = '-g' || ${ARGS[i]} = '--genome' ]]; then # Genome
		GENOME=${ARGS[(i+1)]}
		i=($i+1)
	elif [[ ${ARGS[i]} = '-o' || ${ARGS[i]} = '--output' ]]; then # SNP Output file
		OUTPUT=${ARGS[(i+1)]}
		i=($i+1)
	else
		echo "No switch encountered"
		exit 1
	fi
done

if [[ -z $OUTPUT || -z ${BAMINPUT} || -z $GENOME  ]]; then
	echo "Required arguments not specified."
	exit 2 
fi

if [[ ! -e ${BAMINPUT} || ! -e ${GENOME} ]]; then
	echo "Required input file(s) do not exist."
        exit 3
fi

echo "samtools mpileup -uf $GENOME ${BAMINPUT}.bam | bcftools view -bvcg - > $OUTPUT"
samtools mpileup -uf $GENOME ${BAMINPUT}.bam | bcftools view -bvcg - > $OUTPUT

ERR=$?
if [ ${ERR} -ne 0 ]; then
        echo "Error: samtools mpileup: $ERR"
	exit $ERR
fi

