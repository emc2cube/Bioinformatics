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
	echo "Usage: $(basename "$0") <.fastq(.gz) folder>"
	echo ""
	echo "Description"
	echo ""
	echo "Simple script to consolidate fragmented .fastq files from different sequencing lanes."
	echo "Original files will be backed up in a FastQbackup folder."
	echo ""
	echo "Options:"
	echo "$(basename "$0") --help : Display this help message."
	echo ""
	exit
fi

# Version
if [ "${1}" == "--version" ] || [ "${2}" == "--version" ] || [ "${3}" == "--version" ]
then
	echo "$(basename "$0") version 1.0.1"
	exit
fi

# Get fastq directory
dir="${1}"

# Check paths and trailing / in directories
if [ -z "${dir}" ]
then
	${0} --help
	exit
fi

if [ "${dir: -1}" = "/" ]
then
	dir=${dir%?}
fi

# Test if sequence files are .fastq or .fastq.gz
fastqgz=$(find -L "${dir}" -maxdepth 1 -name '*.fastq.gz')
fastq=$(find -L "${dir}" -maxdepth 1 -name '*.fastq')
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
	gunzip "${dir}/*.fastq.gz"
fi

# Concatenate files if split
echo ""
echo "-- Merging files --"
echo ""
mkdir -p "${dir}/FastQbackup/"
filestomergeR1=$(ls "${dir}"/*_L[0-9][0-9][0-9]_R1)
filestomergeR2=$(ls "${dir}"/*_L[0-9][0-9][0-9]_R2)
echo "Processing R1 files"
echo ""
for i in ${filestomergeR1}
do
	sampleoutput="${i//_L[0-9][0-9][0-9]_R1/_R1}"
	echo "Processing ${i}"
	cat "${i}" >> "${sampleoutput}"
	mv "${i}" "${dir}/FastQbackup/"
done
echo ""
echo "Processing R2 files"
echo ""
for i in ${filestomergeR2}
do
	sampleoutput="${i//_L[0-9][0-9][0-9]_R2/_R2}"
	echo "Processing ${i}"
	cat "${i}" >> "${sampleoutput}"
	mv "${i}" "${dir}/FastQbackup/"
done
echo ""
echo "-- Files merged --"
