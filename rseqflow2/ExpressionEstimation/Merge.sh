#!/bin/bash

# Exit Codes
# 0 Successful
# 2 File(s) not found
# 3 Minimum input files not provided
# 4 Output filename not specified

usage ()
{
	echo "Usage: `basename $0` <-o | --output OUTPUT_FILENAME> <INPUT_1> <INPUT_2> [ INPUT_3 .. INPUT_n ]"
}

declare -a ARGS
declare -a INPUT

ARGS=($@)

last=0
for ((i=0;i<$#;++i))
do
	if [[ ${ARGS[i]} = '-o' || ${ARGS[i]} = '--output' ]]; then
		OUTPUT=${ARGS[(i+1)]}
		unset ARGS[$i]
		i=($i+1)
		unset ARGS[$i]
	else
        INPUT[$last]=${ARGS[(i)]}
        last=($last+1)
	fi
done

# Input array must have at least two entries
if [ ${#INPUT[@]} -le 1 ]; then
	usage
	echo "Error: At least two input files are required"
	exit 3
fi

# Output file name is required
if [ -z $OUTPUT ]; then
	usage
	echo "Error: Output filename not specified"
	exit 4
fi

grep 'Strand' ${INPUT[0]} > $OUTPUT
grep --invert-match 'Strand' --no-filename ${INPUT[@]} >> $OUTPUT

# Ignore if no matches are found.
EC=$?
if [ $EC -le 1 ]; then
	exit 0
fi
