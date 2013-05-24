#!/bin/bash
#
# Usage: sh_bowtie2_AlignAll.sh </path/to/fastq(.gz)/folder> </path/to/Aligned(.bam)/destination/folder> [/path/to/config/file.ini]
#
##############################################################
##                      Description                         ##
##############################################################
#
# This script will process fastq files to .sam files
# First it will trim fastqs with Jason script to remove sequences in the adapter region.
# It will also generate a merged file of all individual samples to be aligned.
# Finally trimmed sequences are aligned to a reference genome using bowtie2
#
##############################################################
##                  Configurable variables                  ##
##############################################################
#
# bowtie2 indexed reference genome location
bt_refgenome="/Tools/RefGenomes/gatk_2.3_ucsc_hg19/gatk_2.3_ucsc_hg19"
#
# Run nproc and get the numbers of all installed CPU
#threads=$(nproc --all)
# Use less processors to allow other tasks to run (n-1 here)
threads=$(nproc --all --ignore=1)
#
# log folder
logs="logs"
#
# Process Undetermined files (tag not properly recognized)
# Yes = 1 ; No = 0
underdet="0"
#
# Process BLANK files (Sample name should be "BLANK" in the MiSeq spreadsheet)
# Yes = 1 ; No = 0
blank="0"
#
# If sample files are split into multiple FastQ files, merge them into the first one.
# Yes = 1 ; No = 0
merge="0"
#
# Read group parameters
# Library
# If empty LB will use the destination folder name.
#LB="Library_Name"
#
# Platform/technology used to produce the reads. Valid values: CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT and PACBIO.
PL="ILLUMINA"
#
## Setup done. You should not need to edit below this point ##

# Get fastq directory
dir="$1"

# Get destination directory
dir2="$2"

# Get config file location
config="$3"

# Check paths and trailing / in directories
if [ -z "$dir" -o -z "$dir2" ]
then
    echo "Usage: sh_bowtie2_AlignAll.sh <.fastq(.gz) folder> <destination folder> [/path/to/config/file.ini]"
    exit
fi

if [ ${dir: -1} == "/" ]
then
    dir=${dir%?}
fi

if [ ${dir2: -1} == "/" ]
then
    dir2=${dir2%?}
fi

if [ -n "$config" ]
then
    if [ ${config: -4} == ".ini" ]
    then
        source "$config"
    else
        echo "Invalid config file detected. Is it an .ini file?"
        echo "Usage: sh_gatkSNPcalling.sh <.bam folder> <destination folder> [/path/to/config/file.ini]"
        exit
    fi
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

# Displaying variables to shell
echo ""
echo "-- Informations --"
echo ""
echo "Using $bt_refgenome as reference genome"
echo "Your $fileext files are located in $dir/"
echo "Sorted and indexed .bam will be created into $dir2/"
if [ $underdet -eq "0" ]
then
    echo "Undetermined files will not be processed"
fi
if [ $blank -eq "0" ]
then
    echo "\"BLANK\" files will not be processed"
fi
if [ $merge -eq "1" ]
then
    echo "Individual files will also be merged into a big \"MERGED\" file"
fi
echo "This computer have" $(nproc --all) "CPUs installed, $threads CPUs will be used"
echo ""
if [ -z "$config" ]
then
    echo "You can change these parameters by using a custom config file"
else
    echo "You are using a custom config file: $config"
fi

# Initialize
[ -f $dir2/temp1 ] && rm $dir2/temp1
[ -f $dir2/temp2 ] && rm $dir2/temp2
[ -f $dir2/ReadFastqs ] && rm $dir2/ReadFastqs
[ -f $dir2/trim1 ] && rm $dir2/trim1
[ -f $dir2/trim2 ] && rm $dir2/trim2
[ -f $dir2/TrimFastqs ] && rm $dir2/TrimFastqs
mkdir -p $dir2/
mkdir -p $dir2/$logs

# Concatenate files if split (HiSeq)
if [ $merge -eq "1" ]
then
echo ""
echo "-- Merging files --"
echo ""
mkdir -p $dir/FastQbackup/
filestomergeR1=`ls $dir/ | grep _L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9]`
filestomergeR2=`ls $dir/ | grep _L[0-9][0-9][0-9]_R2_[0-9][0-9][0-9]`
echo "Processing R1 files"
echo ""
for i in $filestomergeR1
do
    sampleoutput=`echo $i | sed 's/_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9]/_L001_R1_001/g'`
    echo "Processing $i"
    cat $dir/$i >> $dir/$sampleoutput
    mv $dir/$i $dir/FastQbackup/
done
echo ""
echo "Processing R2 files"
echo ""
for i in $filestomergeR2
do
    sampleoutput=`echo $i | sed 's/_L[0-9][0-9][0-9]_R2_[0-9][0-9][0-9]/_L001_R2_001/g'`
    echo "Processing $i"
    cat $dir/$i >> $dir/$sampleoutput
    mv $dir/$i $dir/FastQbackup/
done
echo ""
echo "-- Files merged --"
fi

# Trim data
if [ $underdet -eq "0" ]
then
    if [ $blank -eq "0" ]
    then
# Remove all Undetermined_* and BLANK* files
        ls $dir/ --hide=Undetermined_* --hide=BLANK* | grep R1 > $dir2/trim1
        ls $dir/ --hide=Undetermined_* --hide=BLANK* | grep R2 > $dir2/trim2
        paste $dir2/trim1 $dir2/trim2 > $dir2/TrimFastqs
    else
# Remove all Undetermined_* files
        ls $dir/ --hide=Undetermined_* | grep R1 > $dir2/trim1
        ls $dir/ --hide=Undetermined_* | grep R2 > $dir2/trim2
        paste $dir2/trim1 $dir2/trim2 > $dir2/TrimFastqs
    fi
else
    if [ $blank -eq "0" ]
    then
# Remove all BLANK* files
        ls $dir/ --hide=BLANK* | grep R1 > $dir2/trim1
        ls $dir/ --hide=BLANK* | grep R2 > $dir2/trim2
        paste $dir2/trim1 $dir2/trim2 > $dir2/TrimFastqs
    else
# Process all the files!
        ls $dir/ | grep R1 > $dir2/trim1
        ls $dir/ | grep R2 > $dir2/trim2
       paste $dir2/trim1 $dir2/trim2 > $dir2/TrimFastqs
    fi
fi
echo ""
echo "-- Starting Trimming --"
echo ""
while read j;
do
    read1=`echo $j | cut -d" " -f1`
    read2=`echo $j | cut -d" " -f2`
    out2=`basename $read1 | sed "s/$fileext/.trimlog/g" | sed 's/_R1//g'`
    echo "Trimming" $dir/$read1
    pyadapter_trim.py -a $dir/$read1 -b $dir/$read2 -o $dir2/ >$dir2/$logs/$out2
done < $dir2/TrimFastqs
echo ""
echo "-- Trimming done! --"

# Get files for alignment
ls $dir2/ | grep R1_001.trim > $dir2/temp1
ls $dir2/ | grep R2_001.trim > $dir2/temp2
paste $dir2/temp1 $dir2/temp2 > $dir2/ReadFastqs

# Align data
echo ""
echo "-- Starting Alignment --"
echo ""
while read j;
do
    read1=`echo $j | cut -d" " -f1`
    read2=`echo $j | cut -d" " -f2`
    out=`basename $read1 | sed 's/.fastq/.sam/g' | sed 's/_R1//g'`
    out2=`basename $read1 | sed 's/.fastq/.alignlog/g' | sed 's/_R1//g'`
    if [ -z $LB ]
    then
        LB=`basename $dir2`
    fi
    SM=`echo $read1 | awk -F_L001_ '{print $1}'`
    CN=`head -n 1 $dir2/$read1 | awk -F: '{print $1}' | sed 's/@//'`
    PU=`head -n 1 $dir2/$read1 | awk -F: '{print $3}'`
    echo "Aligning" $dir2/$read1
    bowtie2 -p $threads --rg-id "$LB"_"$SM" --rg CN:$CN --rg LB:$LB --rg PL:$PL --rg PU:$PU --rg SM:$SM -x $bt_refgenome -S $dir2/$out -1 $dir2/$read1 -2 $dir2/$read2 2>$dir2/$logs/$out2
done < $dir2/ReadFastqs
echo ""
echo "-- Alignment done! --"

# Clean .trim.fastq files
rm $dir2/*.trim.fastq

# Clean temporary files
[ -f $dir2/temp1 ] && rm $dir2/temp1
[ -f $dir2/temp2 ] && rm $dir2/temp2
[ -f $dir2/ReadFastqs ] && rm $dir2/ReadFastqs
[ -f $dir2/trim1 ] && rm $dir2/trim1
[ -f $dir2/trim2 ] && rm $dir2/trim2
[ -f $dir2/TrimFastqs ] && rm $dir2/TrimFastqs
