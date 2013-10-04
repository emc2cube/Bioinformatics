#!/bin/bash
#
# Usage: sh_FastQToTrimmedBam.sh </path/to/fastq(.gz)/folder> [/path/to/destination/folder] [/path/to/config/file.ini]
#
## Description ##
#
# Will call sh_bowtie2_AlignAll.sh to convert fastq to .sam
# then sh_samtools_ProcessSams.sh to convert .sam into sorted and indexed .bam
#
##

if [ -z "$1" -o -z "$2" ]
then
    echo "Usage: sh_FastQToTrimmedBam.sh </path/to/fastq(.gz)/folder> [/path/to/destination/folder] [/path/to/config/file.ini]"
    exit
fi

if [ -n "$3" ] && [ "${3: -4}" != ".ini" ]
then
    echo "Invalid config file detected. Is it an .ini file?"
    echo "Usage: sh_FastQToTrimmedBam.sh </path/to/fastq(.gz)/folder> [/path/to/destination/folder] [/path/to/config/file.ini]"
    exit
fi

sh_bowtie2_AlignAll.sh "$1" "$2" "$3"
sh_samtools_ProcessSams.sh "$2" "$3"
