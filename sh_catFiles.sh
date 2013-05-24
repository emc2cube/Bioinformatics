#!/bin/bash
#
# Usage: sh_catFiles.sh </path/to/original/fastq(.gz)/folder> </path/to/merged/fastq(.gz)/folder>
#
##############################################################
##                      Description                         ##
##############################################################
#
# This script will merge all fastq(.gz) from a folder to fastq(.gz) with the same
# name in a second folder.
#
##

# Get fastq directory
dir="$1"

# Get destination directory
dir2="$2"

# Check paths and trailing / in directories
if [ -z "$dir" ]
then
    echo "usage: sh_catFiles.sh <.fastq(.gz) folder> [destination folder]"
    exit
fi

if [ ${dir: -1} == "/" ]
then
    dir=${dir%?}
fi

if [ -z $dir2 ]
then
    dir2="."
fi

if [ ${dir2: -1} == "/" ]
then
    dir2=${dir2%?}
fi

# Test if sequence files are .fastq or .fastq.gz
fastqgz=`ls $dir/ | grep .fastq.gz`
fastq=`ls $dir/ --hide=*.gz | grep .fastq`
if [ -z "${fastqgz}" ] && [ -z "${fastq}" ]
then
    echo ""
    echo "No .fastq or .fastq.gz files are present in $dir/"
    exit
fi
if [ -n "${fastqgz}" ] && [ -n "${fastq}" ]
then
    echo ""
    echo "Both .fastq and .fastq.gz files are present in $dir/"
    echo "Existing .fastq files will now be converted to .fastq.gz files"
    gzip $dir/*.fastq
    fileext=".fastq.gz"
fi
if [ -n "${fastqgz}" ] && [ -z "${fastq}" ]
then
    fileext=".fastq.gz"
fi
if [ -n "${fastq}" ] && [ -z "${fastqgz}" ]
then
    fileext=".fastq"
fi

fastqfiles=`ls $dir/ | grep $fileext`

for i in $fastqfiles
do
    echo "Appending $dir/$i to $dir2/$i"
    cat $dir/$i >> $dir2/$i
done
