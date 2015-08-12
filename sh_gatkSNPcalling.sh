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
# OR
# - create a Family.vcf file with all the family members
#
# usage: sh_gatkSNPcalling.sh <.bam folder> <destination folder> [/path/to/config/file.ini]
#
##############################################################
##                  Configurable variables                  ##
##############################################################
#
# Are these samples part of a same family?
# Yes = 1 ; No = 0
family="0"
#
# PED file to specify family relations if available
ped=""
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
# Send an email to this address when a job is done. Separate by commas for multiple recipients.
# default no email sent
email="/dev/null"
#
# Send this email as if it was from a custom address
# Change it to use it with IFTTT.com for example
# default is user@host
fromemail="`id -un`@`hostname -A`"
#
##############################################################
## Setup done. You should not need to edit below this point ##
##############################################################

# Get bam directory
dir="$1"

# Get destination directory
dir2="$2"

# Get config file location
config="$3"

# Check paths and trailing / in directories
if [ -z "$dir" -o -z "$dir2" ]
then
    echo "Usage: sh_gatkSNPcalling.sh <.bam folder> <destination folder> [/path/to/config/file.ini]"
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
echo "SNPs filtering options: $filterparameters"
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
echo "${mem}GB of memory will be allocated to java"
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

# Look for duplicates reads with picard-tools
if [ $markduplicates -eq "1" ]
then
echo ""
echo "-- Marking duplicates --"
echo ""
# Get files for duplicates marking
dupfiles=`ls $dir/ --hide=*.bai | grep .sorted.bam`

for i in $dupfiles
do
    echo "Processing" $i
    out=`echo $i | sed 's/.bam/.nodup.bam/g'`
    metrics=`echo $i | sed "s/.sorted.bam/.nodup.metrics/g"`
    stdout=`echo $i | sed "s/.sorted.bam/.noduplog/g"`

    # Use picard tools MarkDuplicates with removal of duplicates and index creation options.
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $picard/MarkDuplicates.jar I=$dir/$i O=$dir2/$out METRICS_FILE=$dir2/$logs/$metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true 2>$dir2/$logs/$stdout
done
echo ""
echo "-- Duplicates removed --"
echo ""
else
echo ""
echo "-- Copying files (no duplication analysis) --"
echo ""
# Move files without duplicates marking
mvfiles=`ls $dir/ | grep .sorted.bam`
for i in $mvfiles
do
    out=`echo $i | sed 's/.bam/.nodup.bam/g'`

    # Copy files
    cp $dir/$i $dir2/$out
done
echo ""
echo "-- Files copied --"
echo ""
fi

# Start GATK best practice sample processing

# List bam files to process
nodupfiles=`ls $dir2/ --hide=*.bai | grep .nodup.bam`

for i in $nodupfiles
do
# Generate file names from sample name
    
    samplename=`echo $i | awk -F. '{print $1}'`
    nodupbai=`echo $i | sed 's/.bam/.bai/g'`
    intervals="$samplename.intervals"
    realigned="$samplename.realigned.bam"
    realignedbai="$samplename.realigned.bai"
    matefixed="$samplename.realigned.fixed.bam"
    matefixedbai="$samplename.realigned.fixed.bai"
    recal_data="$samplename.recal.grp"
    recal="$samplename.realigned.fixed.recal.bam"
    recalbai="$samplename.realigned.fixed.recal.bai"
    rawSNP="$samplename.rawsnp.vcf"
    snpmetrics="$samplename.snpmetrics"
    filteredSNP="$samplename.filteredsnps.vcf"
    annovarfile="$samplename.annovar"
    snpssummary="$samplename.snps"
    coverage="$samplename.coverage"

echo ""
echo "-- Local realignment around indels --"
echo ""

echo "Processing" $i

    # Determining (small) suspicious intervals which are likely in need of realignment
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T RealignerTargetCreator -R $fasta_refgenome -I $dir2/$i -o $dir2/$intervals -known $millsgold -known $onekGph1 -L $regions -nt $threads -ped $ped

    # Running the realigner over those intervals
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T IndelRealigner -R $fasta_refgenome -I $dir2/$i -o $dir2/$realigned -targetIntervals $dir2/$intervals -ped $ped

    # When using paired end data, the mate information must be fixed, as alignments may change during the realignment process
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $picard/FixMateInformation.jar I=$dir2/$realigned O=$dir2/$matefixed SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

echo ""
echo "-- Local realignment done --"
echo ""
echo "-- Quality score recalibration --"
echo ""

    echo "Processing" $matefixed
    # Quality score recalibration
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T BaseRecalibrator -R $fasta_refgenome -I $dir2/$matefixed -knownSites $dbSNP -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o $dir2/$recal_data -L $regions -ped $ped
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T PrintReads -R $fasta_refgenome -I $dir2/$matefixed -BQSR $dir2/$recal_data -o $dir2/$recal -ped $ped

echo ""
echo "-- Recalibration done --"

echo ""
echo "-- Calling SNPs --"
echo ""

    date
    echo "Processing" $recal
    # Produce raw SNP calls
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T UnifiedGenotyper -R $fasta_refgenome -I $dir2/$recal -D $dbSNP -o $dir2/$rawSNP -glm BOTH -nt $threads -L $regions -metrics $dir2/$logs/$snpmetrics -ped $ped

    # Filter SNPs
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T VariantFiltration -R $fasta_refgenome -V $dir2/$rawSNP -o $dir2/$filteredSNP --clusterWindowSize 10 --filterExpression "MQ >= 30.0 && MQ < 40.0 || FS > 50.0 && FS <= 60 || HaplotypeScore >= 13.0 && HaplotypeScore < 60.0" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 5.0" --filterName "LowCoverage" --filterExpression "DP >= 5.0 && DP < 10.0" --filterName "ShallowCoverage" --filterExpression "QUAL < 30.0" --filterName "VeryLowQual" --filterExpression "QUAL >= 30.0 && QUAL < 50.0" --filterName "LowQual" --filterExpression "QD < 2.0" --filterName "LowQD" --filterExpression "FS > 60.0" --filterName "HighFSScore" --filterExpression "FS > 50.0 && FS <= 60" --filterName "MedFSScore" --filterExpression "HaplotypeScore > 13.0" --filterName "HighHaplotypeScore" --filterExpression "HaplotypeScore >= 13.0 && HaplotypeScore < 60.0" --filterName "MedHighHaplotypeScore" --filterExpression "MQ < 40.0" --filterName "LowMQ" --filterExpression "MQ >= 30.0 && MQ < 40.0" --filterName "MedMQ" --filterExpression "MQRankSum < -12.5" --filterName "LowMQRankSum" --filterExpression "ReadPosRankSum < -8.0" --filterName "LowReadPosRankSum" -ped $ped

echo ""
echo "-- SNPs called --"

echo ""
echo "-- Annotate SNPs using annovar --"
echo ""

    echo "Processing" $filteredSNP

    # Convert to annovar format from GATK .vcf file
    $annovar/convert2annovar.pl --format vcf4 --includeinfo $dir2/$filteredSNP --outfile $dir2/$annovarfile

    # Annotate using annovar
    $annovar/table_annovar.pl --buildver hg19 $dir2/$annovarfile $annovar/humandb/ --protocol refGene,phastConsElements46way,genomicSuperDups,gwasCatalog,esp6500si_all,1000g2014oct_all,snp138,ljb26_all,clinvar_20150330 --operation g,r,r,r,f,f,f,f,f --otherinfo --outfile $dir2/$snpssummary --remove #No csv output as annovar result file is full of commas by itself.

    # Fixing headers and cleaning files
    sed -i "1s/Otherinfo/`cat $dir2/$filteredSNP | grep CHROM | sed 's/#//g'`/g" $dir2/$snpssummary.hg19_multianno.txt		# Fixing headers to add back sample names in annovar csv output files
	sed -i "s/,/;/g" $dir2/$snpssummary.hg19_multianno.txt		# Remove any potential existing commas and replace them by semi columns
	sed -i "s/\t/,/g" $dir2/$snpssummary.hg19_multianno.txt		# Convert tabs to commas
	sed -i 's/\\x2c//g' $dir2/$snpssummary.hg19_multianno.txt	# ClinVar database is full of \x2c (comma in hexadecimal), clean it up!
	mv $dir2/$snpssummary.hg19_multianno.txt $dir2/$snpssummary.hg19_multianno.csv

echo ""
echo "-- Annotation done. --"

echo ""
echo "-- Computing coverage --"
echo ""

    echo "Processing" $recal
    # Will compute all coverage informations needed
    java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T DepthOfCoverage -R $fasta_refgenome -I $dir2/$recal -o $dir2/$coverage -L $regions -ct 10 -ct 15 -ct 30 -ped $ped

echo ""
echo "-- Coverage computed --"

# Generate list of mutations in ACMG genes
#sh_ACMGfilter.sh $dir2/$snpsfolder $dir2/$snpsfolder

echo ""
echo "-- Sample $samplename is done, results are in $dir2/ --"
echo "----------------"
echo "----------------"
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
#    rm $dir2/$recal ## Last .bam file, after all cleanup steps. Will be used to call SNPs, better to keep it
#    rm $dir2/$recalbai ## .bai for $recal. To be removed or kept depending of $recal
#    rm $dir2/$rawSNP ## Raw SNPs file. Better to keep it if we want to change filter parameters.
#    rm $dir2/"$rawSNP".idx ## .idx SNPs index file of $rawSNP. To be removed or kept depending of $rawSNP
#    rm $dir2/$filteredSNP ## SNPs file after filter step. Used for annotations, better to keep it
#    rm $dir2/"$filteredSNP".idx ## .idx SNPs index file of $filteredSNP. To be removed or kept depending of $filteredSNP
#    rm $dir2/$annovarfile ## Annovar file. Will be used for annotations, better to keep it
    rm $dir2/$coverage ## Huge file! not sure if we should keep it.
    rm $dir2/$coverage.sample_cumulative_coverage_counts ## Other files from coverage calc.
    rm $dir2/$coverage.sample_cumulative_coverage_proportions
    rm $dir2/$coverage.sample_interval_statistics
    rm $dir2/$coverage.sample_interval_summary
    rm $dir2/$coverage.sample_statistics
    rm -rf $dir2/tmp ## Temporary folder used by Java

done

if [ $family -eq "1" ]

then
	echo "-- Processing a family"
	echo ""

	# Initialize
	familyvcffiles=""
	mkdir -p $dir2/$logs
	mkdir -p $dir2/tmp

	#Getting files names
	vcffiles=`ls $dir2/*.filteredsnps.vcf`

	for f in $vcffiles
	do
		familyvcffiles="$familyvcffiles --variant:`cat $f | grep CHROM | awk -F'\t' '{print $10}'` $f"
	done

	# Create an unique vcf file for the family
	java -Xmx"$mem"g -Djava.io.tmpdir=$dir2/tmp -jar $gatk -T CombineVariants -R $fasta_refgenome $familyvcffiles -o $dir2/Family.vcf -L $regions -ped $ped

	echo ""
	echo "-- Annotating Familial SNPs using annovar --"
	echo ""
	annovarfile="Family.annovar"
	snpssummary="Family.snps"
	filteredSNP="Family.vcf"

	# Annotate using annovar
	# Convert to annovar format from GATK .vcf file
	$annovar/convert2annovar.pl --format vcf4old --includeinfo $dir2/$filteredSNP --outfile $dir2/$annovarfile
	
	# Annotate using annovar
	$annovar/table_annovar.pl --buildver hg19 $dir2/$annovarfile $annovar/humandb/ --protocol refGene,phastConsElements46way,genomicSuperDups,gwasCatalog,esp6500si_all,1000g2014oct_all,snp138,ljb26_all,clinvar_20150330 --operation g,r,r,r,f,f,f,f,f --otherinfo --outfile $dir2/$snpssummary --remove #No csv output as annovar would add the otherinfos data as an unique text field, delimited by "", which confuse an import to excel.

	# Fixing headers to add back sample names in annovar txt output file, and convert it to a proper csv
	sed -i "1s/Otherinfo/`cat $dir2/$filteredSNP | grep CHROM | sed 's/#//g'`/g" $dir2/$snpssummary.hg19_multianno.txt		# Add back sample names in annovar output file
	sed -i "s/,/;/g" $dir2/$snpssummary.hg19_multianno.txt		# Remove any potential existing commas and replace them by semi columns
	sed -i "s/\t/,/g" $dir2/$snpssummary.hg19_multianno.txt		# Convert tabs to commas
	sed -i 's/\\x2c//g' $dir2/$snpssummary.hg19_multianno.txt	# ClinVar database is full of \x2c (comma in hexadecimal), clean it up!
	mv $dir2/$snpssummary.hg19_multianno.txt $dir2/$snpssummary.hg19_multianno.csv

	# Remove indermediate files
	#rm $dir2/$annovarfile ## Annovar file. Will be used for annotations, better to keep it
	rm -rf $dir2/tmp ## Temporary folder used by Java

	echo ""
	echo "-- Family annotation done. --"

	# Creating folders for archive
	[ -d $dir2/ALL_SNPs ] && rm -rf $dir2/ALL_SNPs
	[ -f $dir2/ALL_SNPs.tar.gz ] && rm $dir2/ALL_SNPs.tar.gz
	mkdir -p $dir2/ALL_SNPs/
	
else # This is not a family, consider multiple sporadic samples, create an unique csv file for all of them.

	[ -d $dir2/ALL_SNPs ] && rm -rf $dir2/ALL_SNPs
	mkdir -p $dir2/ALL_SNPs/

	echo ""
	echo "-- Merging csv files. --"

	[ -f $dir2/ALL_SNPs.tar.gz ] && rm $dir2/ALL_SNPs.tar.gz
	[ -f $dir2/All_SNPs_merged.csv ] && rm $dir2/All_SNPs_merged.csv
	[ -d $dir2/ALL_SNPs ] && rm -rf $dir2/ALL_SNPs
	mkdir -p $dir2/ALL_SNPs/

	# Merge all .csv files into an All_SNPs_merged.csv file
    csvfiles=`ls $dir2/ | grep .csv`
    echo "Sample,`head -1 $dir2/\`ls $dir2 | grep .csv | head -1\` | sed s/,CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,.*$//g`,CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GENOTYPE" > $dir2/All_SNPs_merged.csv

    for i in $csvfiles
    do
        filename=`echo $i | awk -F. '{print $1}'` 
        tail -n +2 $dir2/$i > $dir2/filecontent
        while read line
        do 
            echo "\"$filename\",$line" >> $dir2/All_SNPs_merged.csv
            done < $dir2/filecontent
    done

    rm $dir2/filecontent

fi

# Final processing of result files
echo ""
echo "-- Final processing of all SNPs files --"

# Listing all csv files

archive=`ls $dir2/*.csv`

# Creating and archive of all csv files for easy download
for i in $archive
do
    cp -rf $i $dir2/ALL_SNPs/`basename $i`
done

tar --remove-files -C $dir2 -pczf $dir2/ALL_SNPs.tar.gz ALL_SNPs

echo ""
echo "-- Archive ready at $dir2/ALL_SNPs.tar.gz --"

# Send you an email when it's ready, nice isn't it?
sendemail -f $fromemail -t $email -u "`basename $(dirname $dir2)` #SNPCalling done" -m "`date`: `basename $(dirname $dir2)` SNPs calling job is completed." ## -a $dir2/ALL_SNPs.tar.gz ## Optional -a parameter: add a file as attachment (/!\ size)

# That's all folks!
echo ""
echo "                       \|/"
echo "                      (@ @)"
echo ".---------------oOO----(_)----OOo---------------."
echo "|             ALL SAMPLES PROCESSED             |"
echo ".-----------------------------------------------."
