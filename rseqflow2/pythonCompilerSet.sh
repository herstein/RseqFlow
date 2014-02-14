#!/bin/bash
CurrentPath=`pwd`

if [ $# -eq 0 ];
then
        echo "No options and arguments!"
        echo "Try option: -h | --help. It will show you the usage."
        exit 1
fi

declare -a ARGS

ARGS=($@)

for ((i=0;i<$#;++i))
do
        if [[ ${ARGS[i]} = '-p' ]]; then 
                NEWPATH=${ARGS[(i+1)]}
                i=($i+2)
	else
		UNKNOWN=${ARGS[i]}
                echo "No switch encountered: $UNKNOWN"
                exit 1
        fi
done	

NEWPATH='\#!'$NEWPATH

cd $CurrentPath/QC_SNP

pyList=`ls *.py`
for pyFile in $pyList
do
	#sed -i '1i\#!/usr/bin/env python' $pyFile
	sed -i '1i\'$NEWPATH'' $pyFile
	sed -i '2d' $pyFile
done

cd $CurrentPath/ExpressionEstimation

pyList=`ls *.py`
for pyFile in $pyList
do
	#sed -i '1i\#!/usr/bin/env python' $pyFile
        sed -i '1i\'$NEWPATH'' $pyFile
        sed -i '2d' $pyFile
	
done

cd $CurrentPath/DE

pyList=`ls *.py`
for pyFile in $pyList
do
	#sed -i '1i\#!/usr/bin/env python' $pyFile
        sed -i '1i\'$NEWPATH'' $pyFile
        sed -i '2d' $pyFile
	
done
