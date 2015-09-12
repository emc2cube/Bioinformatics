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
# Trim sequences? (Remove adapters sequence, take a long time!)
# 0 = No ; 1 = Yes, using Trim Galore http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/ ; 2 = Yes, using trimmomatic http://www.usadellab.org/cms/index.php?page=trimmomatic ; 3 = Yes, using Jason script
trim="0"
#
# bowtie2 indexed reference genome location
bt_refgenome="/Tools/RefGenomes/gatk_2.3_ucsc_hg19/gatk_2.3_ucsc_hg19"
#
# Trimmomatic location
# example:
# Trimmomatic="/media/Userland/Applications/Trimmomatic/"
Trimmomatic="/media/Userland/Applications/Trimmomatic/"
#
# Run nproc and get the numbers of all installed CPU
#threads=$(nproc --all)
# Use less processors to allow other tasks to run (n-1 here)
threads=$(nproc --all --ignore=1)
#
# Amount of memory, in GB to allocate to programs
# example:
# mem="24"
mem="24"
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
# Filter by MapQ quality score to remove poor aligned reads.
# samtools defaut is 0 (no filtering), 10<MapQ<20 is advised
mapq="10"
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
        echo "Usage: sh_bowtie2_AlignAll.sh <.fastq folder> <.bam destination folder> [/path/to/config/file.ini]"
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
if [ -z "$config" ]
then
    echo "You can change the following  parameters by using a custom config file"
else
    echo "You are using a custom config file: $config"
fi
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
if [ $trim -ne "0" ]
then
    echo "Sequence will be trimmed"
fi
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
echo "This computer have" $(nproc --all) "CPUs installed, $threads CPUs will be used"
echo "$memGB of memory will be allocated to the programs"

# Initialize
[ -f $dir2/temp1 ] && rm $dir2/temp1
[ -f $dir2/temp2 ] && rm $dir2/temp2
[ -f $dir2/ReadFastqs ] && rm $dir2/ReadFastqs
[ -f $dir2/trim1 ] && rm $dir2/trim1
[ -f $dir2/trim2 ] && rm $dir2/trim2
[ -f $dir2/TrimFastqs ] && rm $dir2/TrimFastqs
mkdir -p $dir2/
mkdir -p $dir2/$logs
mkdir -p $dir2/tmp

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
    sampleoutput=`echo $i | sed 's/_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9]/_L001_R1/g'`
    echo "Processing $i"
    cat $dir/$i >> $dir/$sampleoutput
    mv $dir/$i $dir/FastQbackup/
done
echo ""
echo "Processing R2 files"
echo ""
for i in $filestomergeR2
do
    sampleoutput=`echo $i | sed 's/_L[0-9][0-9][0-9]_R2_[0-9][0-9][0-9]/_L001_R2/g'`
    echo "Processing $i"
    cat $dir/$i >> $dir/$sampleoutput
    mv $dir/$i $dir/FastQbackup/
done
echo ""
echo "-- Files merged --"
fi

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


if [ $trim -eq "0" ]
then
    # Need to copy and rename fastq
    echo ""
	echo "-- Copying $fileext files to $dir2 --"
	echo ""
    while read j;
    do
        read1=`echo $j | cut -d" " -f1`
        read2=`echo $j | cut -d" " -f2`
        echo "Copying" $dir/$read1
        cp $dir/$read1 $dir2/$read1
        echo "Copying" $dir/$read2
        cp  $dir/$read2 $dir2/$read2
    done < $dir2/TrimFastqs
    echo ""
    echo "-- Copy done! --"
elif [ $trim -eq "1" ]
then
    # Trim data with Trim Galore
    echo ""
	echo "-- Starting Trimming --"
	echo ""
    while read j;
    do
        read1=`echo $j | cut -d" " -f1`
        read2=`echo $j | cut -d" " -f2`
        out1=`echo $read1 | sed "s/$fileext/_val_1.fq/g"`
        out2=`echo $read2 | sed "s/$fileext/_val_2.fq/g"`
        outtr1=`echo $read1 | sed "s/$fileext/.trim.fastq/g"`
        outtr2=`echo $read2 | sed "s/$fileext/.trim.fastq/g"`
        trim_galore --no_report_file --dont_gzip --stringency 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o $dir2/ --paired $dir/$read1 $dir/$read2 #--fastqc_args "-o $dir2/$logs/ --noextract" ## Optional FastQC analyzis
    # Fix file name to be compatible with next steps
    mv $dir2/$out1 $dir2/$outtr1
    mv $dir2/$out2 $dir2/$outtr2
    done < $dir2/TrimFastqs
    fileext=".fastq"
    echo ""
    echo "-- Trimming done! --"
elif [ $trim -eq "2" ]
then
    # Trim data with trimmomatic
    echo ""
	echo "-- Starting Trimming --"
	echo ""
    while read j;
    do
        read1=`echo $j | cut -d" " -f1`
        read2=`echo $j | cut -d" " -f2`
        out1=`echo $read1 | sed "s/$fileext/.trim$fileext/g"`
        out2=`echo $read2 | sed "s/$fileext/.trim$fileext/g"`
#        unpaired1=`echo $dir2/$read1 | sed "s/$fileext/.unpaired$fileext/g"`  # Save unpaired reads
#        unpaired2=`echo $dir2/$read2 | sed "s/$fileext/.unpaired$fileext/g"`  # Save unpaired reads
        unpaired1="/dev/null"  # Unpaired reads are discarded
        unpaired2="/dev/null"  # Unpaired reads are discarded
        date
        echo "Trimming" $dir/$read1
        echo ""
        java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $Trimmomatic/trimmomatic.jar PE -threads $threads -phred33 $dir/$read1 $dir/$read2 $dir2/$out1 $unpaired1 $dir2/$out2 $unpaired2 ILLUMINACLIP:$Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
        # Run fastqc to generate quality control files
#        fastqc -o $dir2/$logs/ --noextract $dir2/$out1 $dir2/$out2
    done < $dir2/TrimFastqs
    echo ""
    echo "-- Trimming done! --"
elif [ $trim -eq "3" ]
then
    # Trim data with Jason script
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
    fileext=".fastq"
    echo ""
    echo "-- Trimming done! --"
fi

# Get files for alignment
ls $dir2/ | grep _R1. > $dir2/temp1
ls $dir2/ | grep _R2. > $dir2/temp2
paste $dir2/temp1 $dir2/temp2 > $dir2/ReadFastqs

# Align data
echo ""
echo "-- Starting Alignment --"
echo ""
while read j;
do
    read1=`echo $j | cut -d" " -f1`
    read2=`echo $j | cut -d" " -f2`
    samout=`basename $read1 | sed "s/$fileext/.sam/g" | sed 's/_L001_R1//g'`
    alignlog=`basename $read1 | sed "s/$fileext/.alignlog/g" | sed 's/_L001_R1//g'`
    bamout=`basename $read1 | sed "s/$fileext/.bam/g" | sed 's/_L001_R1//g'`
    bamsortedout=`basename $read1 | sed "s/$fileext/.sorted/g" | sed 's/_L001_R1//g'`
    if [ -z $LB ]
    then
        LB=`basename $dir2`
    fi
    SM=`echo $read1 | awk -F_L001 '{print $1}'`
    if [ "$fileext" = ".fastq.gz" ]
    then
        CN=`gzip -cd $dir2/$read1 | head -n 1 | awk -F: '{print $1}' | sed 's/@//'`
        PU=`gzip -cd $dir2/$read1 | head -n 1 | awk -F: '{print $3}'`
    else
        CN=`head -n 1 $dir2/$read1 | awk -F: '{print $1}' | sed 's/@//'`
        PU=`head -n 1 $dir2/$read1 | awk -F: '{print $3}'`
    fi
    date
    echo "Processing" $dir2/$read1
    # Perform alignment with bowtie
    bowtie2 -p $threads --phred33 --rg-id "$LB"_"$SM" --rg CN:$CN --rg LB:$LB --rg PL:$PL --rg PU:$PU --rg SM:$SM -x $bt_refgenome -S $dir2/$samout -1 $dir2/$read1 -2 $dir2/$read2 2>$dir2/$logs/$alignlog
    # Clean .fastq files from bam folder
    rm $dir2/$read1
    rm $dir2/$read2
    # Convert .sam to .bam with optional filter on MapQ quality score
    samtools view -bS -q $mapq $dir2/$samout -o $dir2/$bamout
    # Remove .sam file
    rm $dir2/$samout
    # Sort .bam file
    samtools sort $dir2/$bamout $dir2/$bamsortedout
    # Clean temporary unsorted .bam file
    rm $dir2/$bamout
    # Index sorted .bam
    samtools index $dir2/$bamsortedout.bam
done < $dir2/ReadFastqs
echo ""
echo "-- Alignment done! --"

# Clean temporary files
[ -f $dir2/temp1 ] && rm $dir2/temp1
[ -f $dir2/temp2 ] && rm $dir2/temp2
[ -f $dir2/ReadFastqs ] && rm $dir2/ReadFastqs
[ -f $dir2/trim1 ] && rm $dir2/trim1
[ -f $dir2/trim2 ] && rm $dir2/trim2
[ -f $dir2/TrimFastqs ] && rm $dir2/TrimFastqs
[ -d $dir2/tmp ] && rm -rf $dir2/tmp ## Temporary folder used by Java
