#!/bin/bash
#
##############################################################
##                      Description                         ##
##############################################################
#
# This script will process .bam files and compute coverage information
#
# usage: sh_gatkCoverage.sh <.bam folder> <destination folder> [/path/to/config/file.ini]
#
##############################################################
##                  Configurable variables                  ##
##############################################################
#
#
# GATK .jar file location
# example:
# gatk=/media/Userland/Applications/GATK/GenomeAnalysisTK.jar
gatk="/media/Userland/Applications/GATK/GenomeAnalysisTK.jar"
#
# Reference genome (fasta) file.
# It must be the same one that was indexed by bowtie2 for alignment
# example:
# fasta_refgenome="/media/Tools/RefGenomes/gatk_2.3_ucsc_hg19/ucsc.hg19.fasta"
fasta_refgenome="/media/Tools/RefGenomes/gatk_2.3_ucsc_hg19/ucsc.hg19.fasta"
#
# List of the targeted intervals we wanted to sequence. 
# This is necessary as only about 60-70% of all the reads will end up in exonic regions and the rest may align anywhere else in the genome.
# To restrict the output to exonic sequences, generate a file containing for example all the exons plus 50bp at each end for getting splice site information as well.
# This can be done using the UCSC Table Browser (http://genome.ucsc.edu/cgi-bin/hgTables?command=start).
# Choose the hg19 assembly of the human genome and set the track to RefSeq genes and the table to refGene.
# Use BED format as output format and assign the file an appropriate name.
# By clicking on get output, several more options can be made: Choose "create one bed record per Exon plus 50bp at each end" and save the file.
# If we have it, we can also use a .bed file specific to the library preparation kit targets.
# example:
# regions="/media/Tools/GATK/hg19_50bp_RefSeq_refGene.bed"
regions="/media/Tools/GATK/hg19_50bp_RefSeq_refGene.bed"
#
# Run nproc and get the numbers of all installed CPU
#threads=$(nproc --all)
# Use less processors to allow other tasks to run (n-1 here)
threads=$(nproc --all --ignore=1)
#
# Amount of memory, in GB to allocate to java
# example:
# mem="24"
mem="24"
#
# log folder
logs="logs"
#
##############################################################
## Setup done. You should not need to edit below this point ##
##############################################################

# Get fastq directory
dir="$1"

# Get destination directory
dir2="$2"

# Get config file location
config="$3"

# Check paths and trailing / in directories
if [ -z $dir -o -z "$dir2" ]
then
    echo "usage: sh_gatkCoverage.sh <.bam folder> <destination folder> [/path/to/config/file.ini]"
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

if [ ! -z "$config" ]
then
    source "$config"
fi

# Displaying variables to shell
echo ""
echo "-- Informations --"
echo ""
echo "PROGRAMS"
echo "GATK is installed in $gatk"
echo "$fasta_refgenome will be used as reference genome"
echo ""
echo "DATABASES"
echo "Regions of interest are defined in $regions"
echo ""
echo "FILES"
echo "Your .bam files are located in $dir/"
echo "Coverage files will be created into $dir2/"
echo "Log files will be created in $dir2/$logs"
echo ""
echo "HARDWARE"
echo "This computer have" `nproc --all` "CPUs installed, $threads CPUs will be used"
echo "$mem GB of memory will be allocated to java"
echo ""
echo "You can change these parameters by using a custom config file"

# Starting script.
echo ""
echo "                       \|/"
echo "                      (@ @)"
echo ".---------------oOO----(_)----OOo---------------."
echo "|               GATK SNPs CALLING               |"
echo ".-----------------------------------------------."

# Initialize
mkdir -p $dir2/
mkdir -p $dir2/$logs
mkdir -p $dir2/tmp

# Coverage computing
echo ""
echo "-- Computing coverage --"
echo ""
# Get files for coverage
covfiles=`ls $dir2/ --hide=*.bai | grep .realigned.fixed.recal.bam`

for i in $covfiles
do
# Generate file names from sample name
    echo "Processing" $i
    snpsfolder=`echo $i | sed 's/_L001_001.realigned.fixed.recal.bam//g'`
    coverage=`echo $i | sed 's/_L001_001.realigned.fixed.recal.bam/.coverage/g'`

    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T DepthOfCoverage -R $fasta_refgenome -I $dir2/$i -o $dir2/$coverage -L $regions -ct 10 -ct 15 -ct 30
    # Copy the coverage summary file to SNP folder as it will be archived later for easy download
    cp $dir2/$coverage".sample_summary" $dir2/$snpsfolder/

echo ""
echo "-- Coverage computed --"

echo ""
echo "-- Sample $snpsfolder is done, results are in $dir2/$snpsfolder --"
echo "----------------"

# Remove indermediate files
    rm $dir2/$coverage ## Huge file! not sure if we should keep it.
    rm -rf $dir2/tmp ## Temporary folder used by Java

done

# Final processing of result files
echo ""
echo "-- Final processing of all SNPs files --"

[ -d $dir2/ALL_SNPs ] && rm -rf $dir2/ALL_SNPs
[ -f $dir2/ALL_SNPs.tar.gz ] && rm $dir2/ALL_SNPs.tar.gz
[ -d $dir2/ALL_CSVs ] && rm -rf $dir2/ALL_CSVs

# Listing all SNPs folders

archive=`find $dir2/* -maxdepth 0 -mindepth 0 -type d -not -name $logs`

mkdir -p $dir2/ALL_SNPs/
mkdir -p $dir2/ALL_CSVs/

# Merging all exome_summary files into a All_SNPs_merged.csv file
cp $dir2/*/*.snps.exome_summary.csv $dir2/ALL_CSVs/
sh_csvmerge.sh $dir2/ALL_CSVs/ $dir2/ALL_SNPs/
[ -d $dir2/ALL_CSVs ] && rm -rf $dir2/ALL_CSVs

# Creating and archive of all SNPs folders for easy download
for i in $archive
do
    cp -rf $i $dir2/ALL_SNPs/`basename $i`
done

tar --remove-files -C $dir2 -pczf $dir2/ALL_SNPs.tar.gz ALL_SNPs

echo ""
echo "-- Archive ready at $dir2/ALL_SNPs.tar.gz --"

# That's all folks
echo ""
echo "                       \|/"
echo "                      (@ @)"
echo ".---------------oOO----(_)----OOo---------------."
echo "|             ALL SAMPLES PROCESSED             |"
echo ".-----------------------------------------------."
