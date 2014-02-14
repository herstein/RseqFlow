#!/bin/bash

DE_USAGE=$'Usage:

DE.sh de --c1 <condition_1_ID> --c2 <condition_2_ID> --f1 <condition_1_file_1,condition_1_file_2,...condition_1_file_n> --f2 <condition_1_file_1,condition_2_file_2,...condition_2_file_n> -o <output prefix>

Arguments:

--c1 < condition_1_ID >			ID for condition 1
--c2 < condition_2_ID >			ID for condition 2
--f1 < file list for condition_1 >	format: --f1 con1_s1.txt,con1_s2.txt,...
--f2 < file list for condition_2 >	format: --f2 con2_s1.txt,con2_s1.txt,...
-o/--output < output prefix >		prefix of output files
-h/--help				print this usage message'

#JSH 2013-08-20   Image command is not supported at this time
#IMAGE_USAGE=$'Usage:

#DE.sh image -f < file1,file2,... > -o < output_prefix >

#-f/--filelist < file list >	format: -f MyPrefix_+con1-con2+_Table.txt, MyPrefix_+con2-con3+_Table.txt,...
#				*Note: file name must contain "+con2-con2+" where con1 and con2 are the condition Ids.

#-o/--output < output prefix >	prefix of output files
#-h/--help			print this usage message'

#USAGE="Usage:

#DE.sh contains two different function options: 

#1) de		calculates differential expression
#2) image	generate images based on the results of differential expression run with de

#Use option \"de\" to calculate differential expression: 
#$DE_USAGE


#Use option \"image\" to generate images based on the results of differential expression run with de:
#$IMAGE_USAGE"

USAGE="Usage:

DE.sh contains option \"de\" to calculate differential expression. 

$DE_USAGE"

echo ""
echo "You are running: $VERSION"
echo ""

if [ $# -eq 0 ];
then
        echo "No arguments or options!"
        echo "$USAGE"
        exit 1
fi

declare -a ARGS

ARGS=($@)

OUTPUT="DE_output"
Function=${ARGS[0]}
if [ $Function = 'de' ]; then
	for ((i=1;i<$#;++i))
	do
		if [[ ${ARGS[i]} = '--c1' ]]; then
			Condition1=${ARGS[(i+1)]}
			i=($i+1)
		elif [[ ${ARGS[i]} = '--c2' ]]; then
			Condition2=${ARGS[(i+1)]}
       	        	i=($i+1)
		elif [[ ${ARGS[i]} = '--f1' ]]; then
			Con1_Sample=${ARGS[(i+1)]}
			i=($i+1)
			#JSH -- remove "|" for Pegasus
			#List1_line=`echo "$Con1_Sample"|sed "s/,/ /g"`
			List1_line=`sed "s/,/ /g" <<<$Con1_Sample`
	                NSAMPLE1=0
			for fileName in $List1_line
	                do
	                        Con1_list[$NSAMPLE1]=$fileName
	                        NSAMPLE1=$(($NSAMPLE1+1))
	                done
		elif [[ ${ARGS[i]} = '--f2' ]]; then
	                Con2_Sample=${ARGS[(i+1)]}
	                i=($i+1)
			#JSH -- remove "|" for Pegasus
	                #List2_line=`echo "$Con2_Sample"|sed "s/,/ /g"`
			List2_line=`sed "s/,/ /g" <<<$Con2_Sample`
	                NSAMPLE2=0
	                for fileName in $List2_line
	                do
	                        Con2_list[$NSAMPLE2]=$fileName
	                        NSAMPLE2=$(($NSAMPLE2+1))
	                done
		elif [[ ${ARGS[i]} = '-o' || ${ARGS[i]} = '--output' ]]; then
	                OUTPUT=${ARGS[(i+1)]}
			i=($i+1)
		elif [[ ${ARGS[i]} = '-h' || ${ARGS[i]} = '--help' ]]; then
                        echo "$DE_USAGE"
	                exit 0
		else
	                UNKNOWN=${ARGS[i]}
	                echo "Unknown switch encountered: $UNKNOWN"
			echo "$DE_USAGE"
	                exit 1
	        fi
	done

	if [[ -z $Con1_Sample || -z $Con2_Sample || -z $Condition1 || -z $Condition2 ]]; then
		echo "required arguments not specified!"
		echo "$DE_USAGE"
        	exit 1
	fi

	if [[ $NSAMPLE1 -gt 1 || $NSAMPLE2 -gt 1 ]]; then
		is_replicate=true
	else
		is_replicate=false
	fi

	if $is_replicate; then
		DE_Output1="$OUTPUT"_DE_all_WithReplicate_+"$Condition1"-"$Condition2"+_Table.txt
		DE_Output2="$OUTPUT"_DE_Significant_WithReplicate_+"$Condition1"-"$Condition2"+_Table.txt
		DE_Table="$OUTPUT"_WithReplicate_Table.txt
		SampleCombine_WithReplicate.py --c1 $Condition1 --c2 $Condition2 --n1 $NSAMPLE1 --n2 $NSAMPLE2 -1 $Con1_Sample -2 $Con2_Sample -o $DE_Table
		R CMD BATCH --args -filename="$DE_Table" -fdr="0.05" -output1="$DE_Output1"  -output2="$DE_Output2" -con1="$Condition1" -con2="$Condition2" $RSEQFLOWPATH/Deseq_WithReplicate_FDR.r "$OUTPUT"_temp.out
	else
		DE_Output1="$OUTPUT"_DE_all_WithoutReplicate_+"$Condition1"-"$Condition2"+_Table.txt
		DE_Output2="$OUTPUT"_DE_Significant_WithoutReplicate_+"$Condition1"-"$Condition2"+_Table.txt
		DE_Table="$OUTPUT"_WithoutReplicate_Table.txt
		SampleCombine_WithoutReplicate.py --c1 $Condition1 --c2 $Condition2 -1 $List1_line -2 $List2_line -o $DE_Table 
        	R CMD BATCH --args -filename="$DE_Table" -fdr="0.05" -output1="$DE_Output1" -output2="$DE_Output2" $RSEQFLOWPATH/Deseq_WithoutReplicate_FDR.r "$OUTPUT"_temp.out
	fi
	rm "$OUTPUT"_temp.out -f

#JSH 2013-08-20   Image command is not supported at this time
#elif [ $Function = 'image' ]; then
#	for ((i=1;i<$#;++i))
#        do
#                if [[ ${ARGS[i]} = '-f' || ${ARGS[i]} = '--filelist' ]]; then
#                        FileList=${ARGS[(i+1)]}
#                        i=($i+1)
#			#JSH -- remove "|" for Pegasus
#                        #List_line=`echo "$FileList"|sed "s/,/ /g"`
#			List_line=`sed "s/,/ /g" <<<$FileList`
#                        FileNumber=0
#                        for fileName in $List_line
#                        do
#                                DE_Output[$FileNumber]=$fileName
#                                FileNumber=$(($FileNumber+1))
#                        done
#                elif [[ ${ARGS[i]} = '-o' || ${ARGS[i]} = '--output' ]]; then
#                        OUTPUT=${ARGS[(i+1)]}
#                        i=($i+1)
#                elif [[ ${ARGS[i]} = '-h' || ${ARGS[i]} = '--help' ]]; then
#			echo "$IMAGE_USAGE"
##                        echo "-f/--filelist < file list > format: -f MyPrefix_+con1-con2+_Table.txt, MyPrefix_+con2-con3+_Table.txt,..." 
##			echo "              *Note: file name must contains the part '+con2-con2+', con1 and con2 are conditionID."
##                        echo "-o/--output < output prefix > Prefix of output files"
#                        exit 0
#                else
#                        UNKNOWN=${ARGS[i]}
#                        echo "Unknown switch encountered: $UNKNOWN"
#			echo "$IMAGE_USAGE"
##                        echo "Try option: -h/--help. It will show you the usage."
#                        exit 1	
#		fi
#	done
#        if [[ -z $FileList ]]; then
#                echo "required arguments not specified!"
#                echo "$IMAGE_USAGE"
#                exit 1
#        fi

#	DE_Image "$OUTPUT"_DE_Image.png "$OUTPUT"_DE_matrix.txt ${DE_Output[@]}
elif [[ ${ARGS[i]} = '-h' || ${ARGS[i]} = '--help' ]]; then
	echo "$USAGE"
        exit 0
else
	echo "Unrecognized function: $Function"
	echo "$USAGE"
        exit 0
fi
