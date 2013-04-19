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

if [ ! -z "$config" ]
then
    source "$config"
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

# Look for duplicates reads with picard-tools
if [ $markduplicates -eq "1" ]
then
echo ""
echo "-- Marking duplicates --"
echo ""
# Get files for duplicates marking
dupfiles=`ls $dir/ --hide=*.bai | grep .trim.sorted.bam`

for i in $dupfiles
do
    echo "Processing" $i
    out=`echo $i | sed 's/.bam/.nodup.bam/g'`
    metrics=`echo $i | sed "s/.trim.sorted.bam/.nodup.metrics/g"`
    stdout=`echo $i | sed "s/.trim.sorted.bam/.noduplog/g"`

    # Use picard tools MarkDuplicates with removal of duplicates and index creation options.
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $picard/MarkDuplicates.jar I=$dir/$i O=$dir2/$out METRICS_FILE=$dir2/$logs/$metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true 2>$dir2/$logs/$stdout
done
echo ""
echo "-- Duplicates removed --"
else
echo ""
echo "-- Copying files (no duplication analysis) --"
echo ""
# Move files without duplicates marking
mvfiles=`ls $dir/ | grep .trim.sorted.bam`
for i in $mvfiles
do
    out=`echo $i | sed 's/.bam/.nodup.bam/g'`

    # Copy files
    cp $dir/$i $dir2/$out
done
echo ""
echo "-- Files copied --"
fi

# Perform local realignment around known indels
echo ""
echo "-- Local realignment around indels --"
echo ""
# Get files for local realignment
nodupfiles=`ls $dir2/ --hide=*.bai | grep .nodup.bam`

for i in $nodupfiles
do
# Generate file names from sample name
    echo "Processing" $i
    nodupbai=`echo $i | sed 's/.bam/.bai/g'`
    intervals=`echo $i | sed 's/.trim.sorted.nodup.bam/.intervals/g'`
    realigned=`echo $i | sed 's/.trim.sorted.nodup.bam/.realigned.bam/g'`
    realignedbai=`echo $i | sed 's/.trim.sorted.nodup.bam/.realigned.bai/g'`
    matefixed=`echo $i | sed 's/.trim.sorted.nodup.bam/.realigned.fixed.bam/g'`
    matefixedbai=`echo $i | sed 's/.trim.sorted.nodup.bam/.realigned.fixed.bai/g'`
    recal_data=`echo $i | sed 's/.trim.sorted.nodup.bam/.recal.grp/g'`
    recal=`echo $i | sed 's/.trim.sorted.nodup.bam/.realigned.fixed.recal.bam/g'`
    rawSNP=`echo $i | sed 's/.trim.sorted.nodup.bam/.rawsnp.vcf/g'`
    snpmetrics=`echo $i | sed 's/.trim.sorted.nodup.bam/.snpmetrics/g'`
    filteredSNP=`echo $i | sed 's/.trim.sorted.nodup.bam/.filteredsnps.vcf/g'`
    annovarfile=`echo $i | sed 's/.trim.sorted.nodup.bam/.annovar/g'`
    snpsfolder=`echo $i | sed 's/_L001_001.trim.sorted.nodup.bam//g'`
    snpssummary=`echo $i | sed 's/_L001_001.trim.sorted.nodup.bam/.snps/g'`
    coverage=`echo $i | sed 's/.trim.sorted.nodup.bam/.coverage/g'`

    # Determining (small) suspicious intervals which are likely in need of realignment
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T RealignerTargetCreator -R $fasta_refgenome -I $dir2/$i -o $dir2/$intervals -known $millsgold -known $onekGph1 -L $regions -nt $threads

    # Running the realigner over those intervals
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T IndelRealigner -R $fasta_refgenome -I $dir2/$i -o $dir2/$realigned -targetIntervals $dir2/$intervals

    # When using paired end data, the mate information must be fixed, as alignments may change during the realignment process
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $picard/FixMateInformation.jar I=$dir2/$realigned O=$dir2/$matefixed SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

echo ""
echo "-- Local realignment done --"

echo ""
echo "-- Quality score recalibration --"
echo ""

    echo "Processing" $matefixed
    # Quality score recalibration
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T BaseRecalibrator -R $fasta_refgenome -I $dir2/$matefixed -knownSites $dbSNP -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o $dir2/$recal_data -L $regions
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T PrintReads -R $fasta_refgenome -I $dir2/$matefixed -BQSR $dir2/$recal_data -o $dir2/$recal

echo ""
echo "-- Recalibration done --"

echo ""
echo "-- Calling SNPs --"
echo ""

    echo "Processing" $recal
    # Produce raw SNP calls
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T UnifiedGenotyper -R $fasta_refgenome -I $dir2/$recal -D $dbSNP -o $dir2/$rawSNP -glm BOTH -nt $threads -L $regions -metrics $dir2/$logs/$snpmetrics

    # Filter SNPs
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T VariantFiltration -R $fasta_refgenome -V $dir2/$rawSNP -o $dir2/$filteredSNP --clusterWindowSize 10 --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 5 " --filterName "LowCoverage" --filterExpression "QUAL < 30.0 " --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0 " --filterName "LowQual" --filterExpression "QD < 1.5 " --filterName "LowQD" --filterExpression "FS > 150 " --filterName "StrandBias"

echo ""
echo "-- SNPs called --"

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
echo "-- Computing coverage --"
echo ""

    echo "Processing" $recal
    # Will compute all coverage informations needed
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T DepthOfCoverage -R $fasta_refgenome -I $dir2/$recal -o $dir2/$coverage -nt $threads -L $regions
    # Copy the coverage summary file to SNP folder as it will be archived later for easy download
    cp $dir2/$coverage".sample_summary" $dir2/$snpsfolder/

echo ""
echo "-- Coverage computed --"

echo ""
echo "-- Sample $snpsfolder is done, results are in $dir2/$snpsfolder --"
echo "----------------"

# Remove indermediate files
    rm $dir2/$i
    [ -f $dir2/$nodupbai ] && rm $dir2/$nodupbai
    [ -f $dir2/$i.bai ] && rm $dir2/$i.bai
    rm $dir2/$intervals
    rm $dir2/$realigned
    rm $dir2/$realignedbai
    rm $dir2/$recal_data
    rm $dir2/$matefixed
    rm $dir2/$matefixedbai
#    rm $dir2/$recal ## Last .bam file, after al cleanup steps. Will be used to call SNPs, better to keep it
    rm $dir2/$rawSNP
    rm $dir2/"$rawSNP".idx
#    rm $dir2/$filteredSNP ## SNPs file after filter step. Will be used for annotations, better to keep it
#    rm $dir2/"$filteredSNP".idx SNPs index file of $filteredSNP. To be removed or kept depending of $filteredSNP
#    rm $dir2/$annovarfile ## Annovar file. Will be used for annotations, better to keep it
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
