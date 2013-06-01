#!/bin/bash
#
# Usage: sh_samtools_ProcessSams.sh </path/to/sam/files/> [/path/to/config/file.ini]
#
##############################################################
##                      Description                         ##
##############################################################
#
# This script will process .sam files and convert them to .bam files.
# Optional: reads can be filtered on their MapQ score to remove poorly aligned reads.
# These .bam files are then sorted and indexed.
# Intermediate files are deleted to save space.
#
##############################################################
##                  Configurable variables                  ##
##############################################################
#
# Filter by MapQ quality score to remove poor aligned reads.
# samtools defaut is 0 (no filtering), 10<MapQ<20 is advised
mapq="10"
#
## Setup done. You should not need to edit below this point ##

# Check paths and trailing / in directories
dir="$1"

# Get config file location
config="$2"

if [ -z "$dir" ]
then
    echo "Usage: sh_samtools_ProcessSams.sh </path/to/sam/files/> [/path/to/config/file.ini]"
    exit
fi

if [ ${dir: -1} == "/" ]
then
    dir=${dir%?}
fi

if [ -n "$config" ]
then
    if [ ${config: -4} == ".ini" ]
    then
        source "$config"
    else
        echo "Invalid config file detected. Is it an .ini file?"
        echo "Usage: sh_samtools_ProcessSams.sh </path/to/sam/files/> [/path/to/config/file.ini]"
        exit
    fi
fi

# Displaying variables to shell
echo ""
echo "-- Informations --"
echo ""
if [ -z $mapq ]
then
    echo "MapQ score will not be used to filter reads"
else
    if [ $mapq == "0" ]
        then
        echo "MapQ score will not be used to filter reads"
        else
        echo "Reads with a MapQ score <$mapq will be discarded"
    fi
fi
echo ""
if [ -z "$config" ]
then
    echo "You can change these parameters by using a custom config file"
else
    echo "You are using a custom config file: $config"
fi

# Get sam files list
sam_files=`ls $dir | grep .sam`

# Loop through sam files and convert them to filtered, sorted and indexed .bam
echo ""
echo "-- Starting convertion of .sam files to sorted and indexed .bam files --"
echo ""
for i in $sam_files
do
    echo "Processing" $dir/$i

# Generate output files name from input file name
    out1=`echo $i | sed 's/.sam/.bam/'`
    out2=`echo $i | sed 's/.sam/.sorted/'`

# Convert .sam to .bam with optional filter on MapQ quality score
    samtools view -bS -q $mapq $dir/$i -o $dir/$out1

# Remove .sam file
    rm $dir/$i

# Sort .bam file
    samtools sort $dir/$out1 $dir/$out2

# Clean temporary unsorted .bam file
    rm $dir/$out1

# Index sorted .bam
    samtools index $dir/$out2.bam

done

echo ""
echo "-- Conversion done! --"
