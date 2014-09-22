#!/bin/bash

CurrentPath=`pwd`
VERSION="RseqFlow V2.2   2014 Sep 22"

echo "----setup QC_SNP module----"
cd $CurrentPath/QC_SNP
chmod 755 QC_SNP.sh
make
ERR=$?
if [ $ERR -ne 0 ]; then
	echo "Error: failed to make QC_SNP module"
      	exit 1
fi

echo "----setup ucsc_tools----"
# Get the OS and processor
mach=`uname -a`
ERR=$?
if [ $ERR -ne 0 ]; then
        echo "Error: failed to execute 'uname -a' "
        exit 1
fi

if [[ $mach =~ Darwin ]]; then
    if [[ $mach =~ i386 ]]; then
        echo "Linking MacOS i386 versions of ucsc_tools"
        ln -s $CurrentPath/OtherTools/ucsc_tools/genePredToBed-macOSX.i386 genePredToBed
        ln -s $CurrentPath/OtherTools/ucsc_tools/gtfToGenePred-macOSX.i386 gtfToGenePred
	ln -s $CurrentPath/OtherTools/ucsc_tools/genePredToGtf-macOSX.i386 $CurrentPath/OtherTools/ucsc_tools/genePredToGtf
    elif [[ $mach =~ x86_64 ]]; then
        echo "Linking MacOS x86_64 versions of ucsc_tools"
        ln -s $CurrentPath/OtherTools/ucsc_tools/genePredToBed-macOSX.x86_64 genePredToBed
        ln -s $CurrentPath/OtherTools/ucsc_tools/gtfToGenePred-macOSX.x86_64 gtfToGenePred
        ln -s $CurrentPath/OtherTools/ucsc_tools/genePredToGtf-macOSX.x86_64 $CurrentPath/OtherTools/ucsc_tools/genePredToGtf

    else
        echo "Error: failed to link ucsc_tools scripts. Possibly an unsupported Mac processor. Supported processors are i386 and x86_64"  
        exit 1
    fi
else
   echo "Linking Linux versions of ucsc_tools"
   ln -s $CurrentPath/OtherTools/ucsc_tools/genePredToBed-linux.x86_64 genePredToBed
   ln -s $CurrentPath/OtherTools/ucsc_tools/gtfToGenePred-linux.x86_64 gtfToGenePred
   ln -s $CurrentPath/OtherTools/ucsc_tools/genePredToGtf-linux.x86_64 $CurrentPath/OtherTools/ucsc_tools/genePredToGtf
fi


cd $CurrentPath/QC_SNP/lib
tar -xf bx.tar
ERR=$?
if [ $ERR -ne 0 ]; then
        echo "Error: failed to decompress file QC_SNP/lib/bx.tar"
        exit 1
fi
tar -xf bx_extras.tar
ERR=$?
if [ $ERR -ne 0 ]; then
        echo "Error: failed to decompress file QC_SNP/lib/bx_extras.tar"
        exit 1
fi

cd $CurrentPath/QC_SNP/lib/text2pdf
make
ERR=$?
if [ $ERR -ne 0 ]; then
        echo "Error: failed to make QC_SNP/lib/text2pdf module"
        exit 1
fi



echo "----setup ExpressionEstimation module----"
cd $CurrentPath/ExpressionEstimation
chmod 755 ExpressionEstimation.sh
#make
#ERR=$?
#if [ $ERR -ne 0 ]; then
#        echo "Error: failed to make ExpressionEstimation module"
#        exit 1
#fi

#echo "----setup DE module----"
#cd $CurrentPath/DE
#make
#ERR=$?
#if [ $ERR -ne 0 ]; then
#        echo "Error: failed to make DE module"
#        exit 1
#fi

echo "----setup FileFormatConversion modules----"

cd $CurrentPath/OtherTools

echo "----setup Bowtie2----"
cd $CurrentPath/OtherTools
unzip $CurrentPath/OtherTools/bowtie2-2.1.0-linux-x86_64.zip
ERR=$?
if [ $ERR -ne 0 ]; then
	echo "Error: failed to decompress file bowtie2-2.1.0-linux-x86_64.zip"
        exit 1
fi

echo "----setup samtools----"
cd $CurrentPath/OtherTools
tar -xf $CurrentPath/OtherTools/samtools-0.1.18.tar
ERR=$?
if [ $ERR -ne 0 ]; then
        echo "Error: failed to decompress file samtools-0.1.18.tar"
        exit 1
fi
cd $CurrentPath/OtherTools/samtools-0.1.18
make
ERR=$?
if [ $ERR -ne 0 ]; then
        echo "Error: failed to make samtools"
        exit 1
fi

echo "----setup RSEQtools----"
#tar -xf $CurrentPath/OtherTools/RSEQtools-0.6.tar
#$CurrentPath/OtherTools/RSEQtools/mrf/make
#$CurrentPath/OtherTools/RSEQtools/bios/make

echo "----setup BEDTools----"
cd $CurrentPath/OtherTools
tar -xf $CurrentPath/OtherTools/BEDTools.v2.16.2.tar
ERR=$?
if [ $ERR -ne 0 ]; then
        echo "Error: failed to decompress file BEDTools.v2.16.2.tar"
        exit 1
fi
cd $CurrentPath/OtherTools/BEDTools-Version-2.16.2
make
ERR=$?
if [ $ERR -ne 0 ]; then
        echo "Error: failed to make BEDTools"
        exit 1
fi

cd $CurrentPath/DE
Rscript $CurrentPath/import_DESeq_brainwaver.r
ERR=$?
if [ $ERR -ne 0 ]; then
        echo "Error: failed to import R package DESeq"
        exit 1
fi

echo "----Create configure.sh----"
cd $CurrentPath

touch $CurrentPath/configure.sh

echo "PATH=\$PATH:$CurrentPath/QC_SNP" >$CurrentPath/configure.sh

echo "PATH=\$PATH:$CurrentPath/QC_SNP/lib/text2pdf" >>$CurrentPath/configure.sh

echo "PATH=\$PATH:$CurrentPath/ExpressionEstimation" >>$CurrentPath/configure.sh

echo "PATH=\$PATH:$CurrentPath/DE" >>$CurrentPath/configure.sh

echo "PATH=\$PATH:$CurrentPath/FileFormatConversion" >>$CurrentPath/configure.sh

echo "PATH=\$PATH:$CurrentPath/FileFormatConversion" >>$CurrentPath/configure.sh

echo "PATH=\$PATH:$CurrentPath/OtherTools/bowtie2-2.1.0" >>$CurrentPath/configure.sh

echo "PATH=\$PATH:$CurrentPath/OtherTools/samtools-0.1.18" >>$CurrentPath/configure.sh

echo "PATH=\$PATH:$CurrentPath/OtherTools/samtools-0.1.18/bcftools" >>$CurrentPath/configure.sh

echo "PATH=\$PATH:$CurrentPath/OtherTools/RSEQtools" >>$CurrentPath/configure.sh

echo "PATH=\$PATH:$CurrentPath/OtherTools/BEDTools-Version-2.16.2/bin" >>$CurrentPath/configure.sh

echo "PATH=\$PATH:$CurrentPath/OtherTools/ucsc_tools" >>$CurrentPath/configure.sh

echo "PYTHONPATH=\$PYTHONPATH:$CurrentPath/QC_SNP/lib" >>$CurrentPath/configure.sh

echo "RSEQFLOWPATH=$CurrentPath/DE" >>$CurrentPath/configure.sh

echo "VERSION='$VERSION'" >>$CurrentPath/configure.sh

echo "export PATH" >>$CurrentPath/configure.sh

echo "export PYTHONPATH" >>$CurrentPath/configure.sh

echo "export RSEQFLOWPATH" >>$CurrentPath/configure.sh

echo "export VERSION" >>$CurrentPath/configure.sh

chmod +x $CurrentPath/configure.sh
