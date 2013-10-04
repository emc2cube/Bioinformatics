#!/bin/bash
#
##############################################################
##                      Description                         ##
##############################################################
#
# This script will process .bam files
# Optional: Will first remove duplicate reads (edit options).
# This script will, for all samples:
# - perform a local realignment around known indels.
# - perform a quality score recalibration.
# - call SNPs.
# - annotate using annovar
# - merge csv files and do some cleaning for an easy downloadable file
#
# usage: sh_gatkSNPcalling.sh <.bam folder> <destination folder> [/path/to/config/file.ini]
#
##############################################################
##                  Configurable variables                  ##
##############################################################
#
# Picard-tools location
# example:
# picard="/media/Userland/Applications/picard-tools"
picard="/media/Userland/Applications/picard-tools"
#
# GATK .jar file location
# example:
# gatk=/media/Userland/Applications/GATK/GenomeAnalysisTK.jar
gatk="/media/Userland/Applications/GATK/GenomeAnalysisTK.jar"
#
# Annovar folder
# example:
# annovar="/media/Userland/Applications/annovar"
annovar="/media/Userland/Applications/annovar"
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
# Indel reference files (from GATK bundle)
# Mills and 1000 genomes gold standard vcf file location
# example:
# mills="/media/Tools/GATK/Mills_and_1000G_gold_standard.indels.hg19.vcf"
millsgold="/media/Tools/GATK/Mills_and_1000G_gold_standard.indels.hg19.vcf"
#
# 1000 genomes phase 1 indels vcf file location
# example:
# onekGph1="/media/Tools/GATK/1000G_phase1.indels.hg19.vcf"
onekGph1="/media/Tools/GATK/1000G_phase1.indels.hg19.vcf"
#
# dbSNP known SNPs .vcf file location
# example:
# dbSNP="/Tools/GATK/dbsnp_137.hg19.vcf"
dbSNP="/Tools/GATK/dbsnp_137.hg19.vcf"
#
# Mark for duplicates
# In case of exome sequencing (not targetted sequencing) we want to mark and remove duplicates.
# Yes = 1 ; No = 0
markduplicates="0"
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
    echo "usage: sh_gatkSNPcalling.sh <.bam folder> <destination folder> [/path/to/config/file.ini]"
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

# Displaying variables to shell
echo ""
echo "-- Informations --"
echo ""
echo "PROGRAMS"
echo "Picard-tools are installed in $picard"
echo "GATK is installed in $gatk"
echo "ANNOVAR is installed in $annovar"
echo "$fasta_refgenome will be used as reference genome"
echo ""
echo "DATABASES"
echo "Local realigment around known indels will be done using:"
echo "Mills and 1000 genomes gold standard, located at $millsgold"
echo "1000 genomes phase 1 database, located at $onekGph1"
echo ""
echo "SNPs will be filtered using:"
echo "dbSNP database, located at $dbSNP"
echo ""
echo "Regions of interest are defined in $regions"
echo ""
echo "FILES"
echo "Your .bam files are located in $dir/"
echo "SNP calling files will be created into $dir2/"
echo "Log files will be created in $dir2/$logs"
if [ $markduplicates -eq "1" ]
then
echo "Duplicates will be removed"
else
echo "Duplicates won't be removed"
fi
echo ""
echo "HARDWARE"
echo "This computer have" `nproc --all` "CPUs installed, $threads CPUs will be used"
echo "$mem GB of memory will be allocated to java"
echo ""
if [ -z "$config" ]
then
    echo "You can change these parameters by using a custom config file"
else
    echo "You are using a custom config file: $config"
fi

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

# Get files for filtering
filterfiles=`ls $dir/ --hide=*.idx | grep .rawsnp.vcf`

echo ""
echo "-- Filtering SNPs --"
echo ""
for i in $filterfiles
do
# Generate file names from sample name
    echo "Processing" $i
    snpmetrics=`echo $i | sed 's/.rawsnp.vcf/.snpmetrics/g'`
    filteredSNP=`echo $i | sed 's/.rawsnp.vcf/.filteredsnps.vcf/g'`
    annovarfile=`echo $i | sed 's/.rawsnp.vcf/.annovar/g'`
    snpsfolder=`echo $i | sed 's/.rawsnp.vcf//g'`
    snpssummary=`echo $i | sed 's/.rawsnp.vcf/.snps/g'`

    # Filter SNPs
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T VariantFiltration -R $fasta_refgenome -V $dir/$i -o $dir2/$filteredSNP --clusterWindowSize 10 --filterExpression "MQ >= 30.0 && MQ < 40.0 || FS > 50.0 && FS <= 60 || HaplotypeScore >= 13.0 && HaplotypeScore < 60.0" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 5.0" --filterName "LowCoverage" --filterExpression "DP >= 5.0 && DP < 10.0" --filterName "ShallowCoverage" --filterExpression "QUAL < 30.0" --filterName "VeryLowQual" --filterExpression "QUAL >= 30.0 && QUAL < 50.0" --filterName "LowQual" --filterExpression "QD < 2.0" --filterName "LowQD" --filterExpression "FS > 60.0" --filterName "HighFSScore" --filterExpression "FS > 50.0 && FS <= 60" --filterName "MedFSScore" --filterExpression "HaplotypeScore > 13.0" --filterName "HighHaplotypeScore" --filterExpression "HaplotypeScore >= 13.0 && HaplotypeScore < 60.0" --filterName "MedHighHaplotypeScore" --filterExpression "MQ < 40.0" --filterName "LowMQ" --filterExpression "MQ >= 30.0 && MQ < 40.0" --filterName "MedMQ" --filterExpression "MQRankSum < -12.5" --filterName "LowMQRankSum" --filterExpression "ReadPosRankSum < -8.0" --filterName "LowReadPosRankSum"

echo ""
echo "-- SNPs filtered --"

echo ""
echo "-- Annotate SNPs using annovar --"
echo ""

    echo "Processing" $filteredSNP

    # Annotate using annovar
    # Convert to annovar format from GATK .vcf file
    convert2annovar.pl --format vcf4 --includeinfo $dir2/$filteredSNP --outfile $dir2/$annovarfile

    # Annotate using annovar
    mkdir -p $dir2/$snpsfolder
    summarize_annovar.pl --buildver hg19 $dir2/$annovarfile $annovar/humandb -outfile $dir2/$snpsfolder/$snpssummary -verdbsnp 137 -ver1000g 1000g2012apr -veresp 6500

echo ""
echo "-- Annotation done. --"

echo ""
echo "-- Sample $snpsfolder is done, results are in $dir2/$snpsfolder --"
echo "----------------"

# Remove indermediate files
#    rm $dir2/$recal ## Last .bam file, after al cleanup steps. Will be used to call SNPs, better to keep it
#    rm $dir2/$rawSNP
#    rm $dir2/"$rawSNP".idx
#    rm $dir2/$filteredSNP ## SNPs file after filter step. Will be used for annotations, better to keep it
#    rm $dir2/"$filteredSNP".idx SNPs index file of $filteredSNP. To be removed or kept depending of $filteredSNP
#    rm $dir2/$annovarfile ## Annovar file. Will be used for annotations, better to keep it
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
