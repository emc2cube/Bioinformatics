#!/bin/bash
#
# Usage: sh_mergeFastQ.sh <.fastq(.gz) folder>
#
##############################################################
##                       Description                        ##
##############################################################
#
# Simple script to consolidate fragmented .fastq files from different sequencing lanes.
# Original files will be backed up in a FastQbackup folder.
#
##

# Help!
if [ "${1}" == "--help" ] || [ "${2}" == "--help" ] || [ "${3}" == "--help" ]
then
	echo "Usage: $(basename $0) <.fastq(.gz) folder>"
	echo ""
	echo "Description"
	echo ""
	echo "Simple script to consolidate fragmented .fastq files from different sequencing lanes."
	echo "Original files will be backed up in a FastQbackup folder."
	echo ""
	echo "Options:"
	echo "$(basename $0) --help : Display this help message."
	echo ""
	exit
fi

# Version
if [ "${1}" == "--version" ] || [ "${2}" == "--version" ] || [ "${3}" == "--version" ]
then
	echo "$(basename $0) version 1.0"
	exit
fi

# Get fastq directory
dir="${1}"


# Check paths and trailing / in directories
if [ -z "${dir}" ]
then
	$(echo "${0} --help")
	exit
fi

if [ ${dir: -1} == "/" ]
then
	dir=${dir%?}
fi


# Test if sequence files are .fastq or .fastq.gz
fastqgz=$(ls ${dir}/ | grep .fastq.gz)
fastq=$(ls ${dir}/ --hide=*.gz | grep .fastq)
if [ -z "${fastqgz}" ] && [ -z "${fastq}" ]
then
	echo ""
	echo "No .fastq or .fastq.gz files are present in ${dir}/"
	exit
fi
if [ -n "${fastqgz}" ] && [ -n "${fastq}" ]
then
	echo ""
	echo "Both .fastq and .fastq.gz files are present in ${dir}/"
	echo "Existing .fastq.gz files will now be converted to .fastq files"
	gunzip ${dir}/*.fastq
	fileext=".fastq"
fi
if [ -n "${fastqgz}" ] && [ -z "${fastq}" ]
then
	fileext=".fastq.gz"
fi
if [ -n "${fastq}" ] && [ -z "${fastqgz}" ]
then
	fileext=".fastq"
fi


# Merging files
echo ""
echo "-- Merging files --"
echo ""
mkdir -p ${dir}/FastQbackup/
filestomergeR1=$(ls ${dir}/ | grep _L[0-9][0-9][0-9]_R1)
filestomergeR2=$(ls ${dir}/ | grep _L[0-9][0-9][0-9]_R2)
echo "Processing R1 files"
echo ""
for i in ${filestomergeR1}
do
	sampleoutput=$(echo ${i} | sed 's/_L[0-9][0-9][0-9]_R1/_R1/g')
	echo "Processing ${i} -> ${sampleoutput}"
	cat ${dir}/${i} >> ${dir}/${sampleoutput}
	mv ${dir}/${i} ${dir}/FastQbackup/
done
echo ""
echo "Processing R2 files"
echo ""
for i in ${filestomergeR2}
do
	sampleoutput=$(echo ${i} | sed 's/_L[0-9][0-9][0-9]_R2/_R2/g')
	echo "Processing ${i} -> ${sampleoutput}"
	cat ${dir}/${i} >> ${dir}/${sampleoutput}
	mv ${dir}/${i} ${dir}/FastQbackup/
done
echo ""
echo "-- Files merged --"
