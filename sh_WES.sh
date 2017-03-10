#!/bin/bash
#
# Usage: sh_WES.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]
#
##############################################################
##					  Description					 	    ##
##############################################################
#
# This script will process fastq(.gz) files and align them to
# a reference genome using bowtie2. It will then use Picard
# and GATK following GATK according to June 2016 best
# practices workflow.
# SNPs will then be annotated with ANNOVAR.
#
##############################################################
##				  Configurable variables				    ##
##############################################################
#
## Hardware options
#
# Maximum number of threads (or CPUs) to request and allocate to programs.
# In some case less than this value may automatically be allowed.
threads=`nproc --all --ignore=1`
#
# Maximum amount of memory (in GB) to request and allocate to programs.
# In some case less than this value may automatically be allowed.
mem="24"
#
# Log folder, will be created in your destination folder.
logs="logs"
#
# Advanced: Path to temporary folder. Useful if your cluster system have a fast local I/O disk.
# Leave empty to use the destination folder. In all cases temporary files will be purged at the end of the script.
tmp=""
#
## SLURM options
#
# email to use to receive SLURM notifications.
SLURMemail=""
#
# SLURM account, partition and qos settings if required.
SLURMaccount=""
SLURMpartition=""
SLURMqos=""
#
# Custom commands to run before any other program.
# Use it to load modules for example
# customcmd="module load R"
customcmd=""
#
## FastQ options
#
# If sample files are split into multiple FastQ files, merge them.
# 0 = No ; 1 = Yes
merge="0"
#
# Perform fastqc on fastq files?
# (If sequences are trimmed, fastqc will be performed on trimmed sequences)
# 0 = No ; 1 = Yes
fastqc="0"
#
# Trim sequences with Trimmomatic?
# 0 = No ; 1 = Yes
trim="0"
#
# Trimmomatic location
Trimmomatic="/bin/Trimmomatic/"
#
# Keep unpaired trimmed sequences?
# 0 = No ; 1 = Yes
unpaired="0"
#
#
## bowtie2 options
#
# bowtie2 indexed reference genome location
bt_refgenome="/Tools/RefGenomes/gatk_2.3_ucsc_hg19/gatk_2.3_ucsc_hg19"
#
# Process Undetermined files (tag not properly recognized)
# 0 = No ; 1 = Yes
underdet="0"
#
# Process BLANK files (Sample name should be "BLANK")
# 0 = No ; 1 = Yes
blank="0"
#
# Read group parameters
# Library: if empty LB will be the destination folder name.
LB=""
#
# Platform used to produce the reads. Valid values: CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT and PACBIO.
PL="ILLUMINA"
#
# Filter by MapQ quality score to remove poor aligned reads.
# samtools defaut is 0 (no filtering), 10<MapQ<20 is advised
mapq="10"
#
#
## GATK options
#
# Picard-tools picard.jar location
picard="/media/Userland/Applications/picard.jar"
#
# GATK GenomeAnalysisTK.jar file location
gatk="/media/Userland/Applications/GenomeAnalysisTK.jar"
#
# Reference genome (fasta) file.
# It must be the same one that was indexed by bowtie2 for alignment
fasta_refgenome="/media/Tools/RefGenomes/gatk_2.3_ucsc_hg19/ucsc.hg19.fasta"
#
# List of the targeted intervals we wanted to sequence. 
# This is necessary as only about 60-70% of all the reads will end up in exonic regions and the rest may align anywhere else in the genome.
# To restrict the output to exonic sequences, generate a file containing for example all the exons plus 50bp at each end for getting splice site information as well.
# This can be done using the UCSC Table Browser (http://genome.ucsc.edu/cgi-bin/hgTables?command=start).
# Choose the hg19 assembly of the human genome and set the track to RefSeq genes and the table to refGene.
# Use BED format as output format and assign the file an appropriate name.
# By clicking on get output, several more options can be made: Choose "create one bed record per Exon plus 50bp at each end" and save the file.
# If you have it, you can also use a .bed file corresponding to the library preparation kit used.
regions="/Tools/AgilentBedFiles/SureSelect_Human_All_Exon_V6+UTR_r2/S07604624_Covered.bed"
#
# VQSR reference files (from GATK bundle):
# HapMap vcf file location (SNP)
hapmap="/Tools/GATK/hapmap_3.3.hg19.sites.vcf"
# 1000 genome omni vcf file location (SNP)
omni="/Tools/GATK/1000G_omni2.5.hg19.sites.vcf"
# 1000 genomes phase 1 indels vcf file location
onekGph1="/Tools/GATK/1000G_phase1.snps.high_confidence.vcf"
# dbSNP known SNPs .vcf file location
dbSNP="/Tools/GATK/dbsnp_137.hg19.vcf"
# Mills and 1000 genomes gold standard vcf file location
millsgold="/Tools/GATK/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
#
# Indel realignment is generally no longer necessary for variant discovery
# unless in specific non common cases.
# 0 = No ; 1 = Yes
realign="0"
#
# Should the script stop after .g.vcf creation? (No filtering, no annotation)
# 0 = No ; 1 = Yes
gvcf="0"
#
# PED file to specify family relations if available
ped=""
#
# Folder containing gVCFs to be used as a learning control for VQSR calibration
# Leave empty if you have a ready to use .vcf file and use popgvcf below.
popfolder=""
# VCF file to be used as a learning control for VQSR calibration.
popgvcf=""
#
# Should coverage be computed? (slow)
# 0 = No ; 1 = Yes
coverage="0"
#
#
## Annovar options
#
# Annovar folder (containing humandb folder)
annovar="/bin/annovar"
#
# Genome build version
# buildver="hg19"
buildver="hg19"
#
# Matching protocols to be used for annotation
protocol="refGene,phastConsElements46way,genomicSuperDups,gwasCatalog,esp6500siv2_all,1000g2015aug_all,exac03,gme,kaviar_20150923,avsnp147,dbnsfp33a,dbnsfp31a_interpro,dbscsnv11,clinvar_20170130,intervar_20170202,hrcr1,revel,mcap"
#
# Matching operations to be used for annotation
operation="g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f"
#
## IFTTT options
#
# Trigger IFTTT when script is done.
# You must register the "Maker channel" on https://ifttt.com/maker
# Copy your private key here. Leave blank to disable this function.
# iftttkey="AbCd_15CdhUIvbsFJTHGMcfgjsdHRTgcyjt" # Not a real key, you have to use your own private key
iftttkey=""
#
# Event name used in your IFTTT recipes.
# The maker channel will look for the combination private key + event name to then trigger your recipe
# You can create a recipe to send an email, a text message or a push notification.
# iftttevent="SNPCalling"
iftttevent="SNPCalling"
#
#
## Setup done. You should not need to edit below this point ##

# Get fastq directory
dir="$1"

# Get destination directory
dir2="$2"

# Get config file location
config="$3"

# Check paths and trailing / in directories
if [ -z "${dir}" -o -z "${dir2}" ]
then
	echo "Usage: sh_WES.sh </path/to/fastq(.gz)/folder> </path/to/SNPsCalled/folder> [/path/to/config/file.ini]"
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

if [ -n "${config}" ]
then
	if [ ${config: -4} == ".ini" ]
	then
		source "${config}"
	else
		echo "Invalid config file detected. Is it an .ini file?"
		echo "Usage: sh_WES.sh </path/to/fastq(.gz)/folder> </path/to/SNPsCalled/folder> [/path/to/config/file.ini]"
		exit
	fi
fi

if [ -n "${popgvcf}" ]
then
	if [ ${popgvcf: -4} != ".vcf"  ]
	then
		echo "Invalid population vcf file detected. Is ${popgvcf} a .vcf file?"
		exit
	fi
fi

if [ -n "${popgvcf}" ] && [ -n "${popfolder}" ]
then
	echo "Error, either \"popgvcf\" or \"popfolder\" should be defined"
	exit
fi

# Test if sequence files are .fastq or .fastq.gz
fastqgz=`ls ${dir}/ | grep .fastq.gz`
fastq=`ls ${dir}/ --hide=*.gz | grep .fastq`
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
	gunzip ${dir}/*.fastq
	fileext=".fastq"
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
if [ -z "${config}" ]
then
	echo "You can change the following  parameters by using a custom config file"
else
	echo "You are using a custom config file: ${config}"
fi
echo "Your ${fileext} files are located in ${dir}/"
echo "Final results will be located in ${dir2}/"
echo "Log files will be created in ${dir2}/${logs}"

echo ""
echo "-- Alignment options"
echo ""
echo "Using ${bt_refgenome} as reference genome"
if [ ${underdet} -eq "0" ]
then
	echo "Undetermined files will not be processed"
fi
if [ ${blank} -eq "1" ]
then
	echo "\"BLANK\" files will be processed"
fi
if [ ${merge} -eq "1" ]
then
	echo "Individual files will also be merged into a big \"MERGED\" file"
fi
if [ ${trim} -eq "1" ]
then
	echo "Sequence will be trimmed"
	echo "Trimmomatic is installed in ${Trimmomatic}"
fi
if [ ${mapq} == "0" ] || [ -z ${mapq} ]
then
	echo "MapQ score will not be used to filter reads"
else
	echo "Reads with a MapQ score <${mapq} will be discarded"
fi

echo ""
echo "-- GATK options"
echo ""
echo "Picard-tools are installed in ${picard}"
echo "GATK is installed in ${gatk}"
echo "ANNOVAR is installed in ${annovar}"
echo "Using ${fasta_refgenome} as reference genome"
echo ""
if [ ${realign} -eq "1" ]
then
	echo "Local realigment around known indels will be performed"
fi
echo "Joint Genotyping will be performed using:"
echo "Mills and 1000 genomes gold standard, located at ${millsgold}"
echo "1000 genomes phase 1 database, located at ${onekGph1}"
echo ""
echo "VQSR calibration will be performed using:"
if [ -n "${popgvcf}" ]
then
	echo "A previously generated file, ${popgvcf}"
else
	echo "A newly generated file using g.vcf files present in ${popfolder}"
fi
echo ""
echo "Regions of interest are defined in ${regions}"
echo ""
echo "-- Hardware"
echo ""
echo "Up to ${threads} CPUs will be used"
echo "Up to ${mem}GB of memory will be allocated to the programs"

# Initialize
[ -f ${dir2}/files1 ] && rm ${dir2}/files1
[ -f ${dir2}/files2 ] && rm ${dir2}/files2
[ -f ${dir2}/Fastqs ] && rm ${dir2}/Fastqs
[ -d ${dir2}/logs ] && rm -rf ${dir2}/logs
mkdir -p ${dir2}/
mkdir -p ${dir2}/${logs}

if [ -z ${tmp} ]
then
	tmp="${dir2}/tmp"
fi

# Concatenate files if split
if [ ${merge} -eq "1" ]
then

	echo ""
	echo "-- Merging files --"
	echo ""
	mkdir -p ${dir}/FastQbackup/
	filestomergeR1=`ls ${dir}/ | grep _L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9]`
	filestomergeR2=`ls ${dir}/ | grep _L[0-9][0-9][0-9]_R2_[0-9][0-9][0-9]`
	echo "Processing R1 files"
	echo ""
	for i in ${filestomergeR1}
	do
		sampleoutput=`echo ${i} | sed 's/_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9]/_R1/g'`
		echo "Processing ${i}"
		cat ${dir}/${i} >> ${dir}/${sampleoutput}
		mv ${dir}/${i} ${dir}/FastQbackup/
	done
	echo ""
	echo "Processing R2 files"
	echo ""
	for i in ${filestomergeR2}
	do
		sampleoutput=`echo ${i} | sed 's/_L[0-9][0-9][0-9]_R2_[0-9][0-9][0-9]/_R2/g'`
		echo "Processing ${i}"
		cat ${dir}/${i} >> ${dir}/${sampleoutput}
		mv ${dir}/${i} ${dir}/FastQbackup/
	done
	echo ""
	echo "-- Files merged --"
fi

if [ ${underdet} -eq "0" ]
then
	if [ ${blank} -eq "0" ]
	then
		# Remove all Undetermined_* and BLANK* files
		ls ${dir}/ --hide=Undetermined_* --hide=BLANK* | grep R1 > ${dir2}/files1
		ls ${dir}/ --hide=Undetermined_* --hide=BLANK* | grep R2 > ${dir2}/files2
		paste ${dir2}/files1 ${dir2}/files2 > ${dir2}/Fastqs
	else
		# Remove all Undetermined_* files
		ls ${dir}/ --hide=Undetermined_* | grep R1 > ${dir2}/files1
		ls ${dir}/ --hide=Undetermined_* | grep R2 > ${dir2}/files2
		paste ${dir2}/files1 ${dir2}/files2 > ${dir2}/Fastqs
	fi
else
	if [ ${blank} -eq "0" ]
	then
		# Remove all BLANK* files
		ls ${dir}/ --hide=BLANK* | grep R1 > ${dir2}/files1
		ls ${dir}/ --hide=BLANK* | grep R2 > ${dir2}/files2
		paste ${dir2}/files1 ${dir2}/files2 > ${dir2}/Fastqs
	else
		# Process all the files!
		ls ${dir}/ | grep R1 > ${dir2}/files1
		ls ${dir}/ | grep R2 > ${dir2}/files2
		paste ${dir2}/files1 ${dir2}/files2 > ${dir2}/Fastqs
	fi
fi

######################
## Processing samples
######################

echo ""
echo "-- Queuing sample jobs --"
echo ""

while read line;
do

	# General variables
	read1=`echo ${line} | cut -d" " -f1`
	read2=`echo ${line} | cut -d" " -f2`
	samplename=`echo ${read1} | awk -F_R1 '{print $1}'`

	echo "Processing" ${samplename}
	


	###################### 
	## Prepare and queue trimming job
	job="trim"

	# Trimmomatic variables
	trimout1="`basename ${read1} ${fileext}`.trim.fastq"
	trimout2="`basename ${read2} ${fileext}`.trim.fastq"
	if [ ! -z ${unpaired} ] && [ ${unpaired} -eq "1" ]
	then
		unpaired1="${dir2}/`basename ${read1} ${fileext}`.unpaired.fastq"  # Save unpaired reads
		unpaired2="${dir2}/`basename ${read2} ${fileext}`.unpaired.fastq"  # Save unpaired reads
	else
		unpaired1="/dev/null"  # Unpaired reads are discarded
		unpaired2="/dev/null"  # Unpaired reads are discarded
	fi
	
	if [ ! -z ${trim} ] && [ ${trim} -eq "1" ]
	then
		
		# General SLURM parameters
		echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH `if [ ! -z ${mem} ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi` `if [ ! -z ${threads} ] && [ ${threads} -gt "8" ]; then echo "--cpus-per-task=8"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH --time=1:00:00" >> ${dir2}/${samplename}_${job}.sbatch
		if [ -n "${SLURMemail}" ]
		then
			echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
		if [ -n "${SLURMpartition}" ]
		then
			echo "#SBATCH --partition=${SLURMpartition}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
		if [ -n "${SLURMqos}" ]
		then
			echo "#SBATCH --qos=${SLURMqos}" >> ${dir2}/${samplename}_${job}.sbatch
		fi

		# General commands
		echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
		
		# Job specific commands
		# Trim fastq files with trimmomatic
		echo "java -Xmx`if [ ! -z ${mem} ] && [ ${mem} -gt "8" ]; then echo "8"; else echo "${mem}"; fi`g -Djava.io.tmpdir=${tmp} -jar ${Trimmomatic}/trimmomatic.jar PE -threads `if [ ! -z ${threads} ] && [ ${threads} -gt "8" ]; then echo "2"; else echo "${threads}"; fi` -phred33 ${dir}/${read1} ${dir}/${read2} ${dir2}/${trimout1} ${unpaired1} ${dir2}/${trimout2} ${unpaired2} ILLUMINACLIP:${Trimmomatic}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch

		# Cleaning commands
		# remove .sbatch
		echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
		
		# Queue job
		SBtrim=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
		echo -e "\t Trimming job queued"
	else
		# Need to rename some variables
		trimout1="${read1}"
		trimout2="${read2}"
	fi
	
	

	######################
	## Prepare and queue alignment job
	job="align"
	
	# bowtie2 variables
	samout=`basename ${read1} | sed "s/_R1${fileext}/.sam/g"`
	bamout=`basename ${read1} | sed "s/_R1${fileext}/.bam/g"`
	bamsortedout=`basename ${read1} | sed "s/_R1${fileext}/.sorted.bam/g"`
	if [ -z $LB ]
	then
		LB=`basename ${dir2}`
	fi
	SM=`echo ${read1} | awk -F_ '{print $1}'`
	if [ "${fileext}" = ".fastq.gz" ]
	then
		CN=`gzip -cd ${dir}/${read1} | head -n 1 | awk -F: '{print $1}' | sed 's/@//'`
		PU=`gzip -cd ${dir}/${read1} | head -n 1 | awk -F: '{print $3}'`
	else
		CN=`head -n 1 ${dir}/${read1} | awk -F: '{print $1}' | sed 's/@//'`
		PU=`head -n 1 ${dir}/${read1} | awk -F: '{print $3}'`
	fi
	
	# General SLURM parameters
	echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --time=8:00:00" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${SLURMemail}" ]
	then
		echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMaccount}" ]
	then
		echo "#SBATCH --account=${SLURMaccount}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMpartition}" ]
	then
		echo "#SBATCH --partition=${SLURMpartition}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMqos}" ]
	then
		echo "#SBATCH --qos=${SLURMqos}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# Require previous job successful completion
	if [ -n "${SBtrim}" ]
	then
		echo "#SBATCH --dependency=afterok:${SBtrim##* }" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# General commands
	echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# Job specific commands
	# Perform alignment with bowtie
	echo "bowtie2 -p ${threads} --phred33 --rg-id ${LB}_${SM} --rg CN:${CN} --rg LB:${LB} --rg PL:${PL} --rg PU:${PU} --rg SM:${SM} -x ${bt_refgenome} -S ${dir2}/${samout} -1 ${dir2}/${trimout1} -2 ${dir2}/${trimout2} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch
	# Convert .sam to .bam with optional filter on MapQ quality score
	echo "samtools view -bS -q ${mapq} -o ${dir2}/${bamout} ${dir2}/${samout} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch
	# Sort .bam file
	echo "samtools sort -@ ${threads} -o ${dir2}/${bamsortedout} -O bam -T ${tmp} ${dir2}/${bamout} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch
	# Index sorted .bam
	echo "samtools index ${dir2}/${bamsortedout} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch

	# Cleaning commands
	# remove trim.fastq files from destination folder
	if ([ ! -z ${trim} ] && [ ${trim} -eq "1" ]) && ([ -z ${fastqc} ] || [ ${fastqc} -eq "0" ])
	then
		echo "rm ${dir2}/${trimout1} ${dir2}/${trimout2}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	# remove .sam file
	echo "rm ${dir2}/${samout}" >> ${dir2}/${samplename}_${job}.sbatch
	# remove unsorted .bam file
	echo "rm ${dir2}/${bamout}" >> ${dir2}/${samplename}_${job}.sbatch
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBalign=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	echo -e "\t Alignment job queued"



	###################### 
	## Run fastqc to generate quality control files
	if [ ! -z ${fastqc} ] && [ ${fastqc} -eq "1" ]
	then	
		job="fqc"
	
		# General SLURM parameters
		echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH `if [ ! -z ${mem} ] && [ ${mem} -gt "4" ]; then echo "--mem=4000"; else echo "--mem=${mem}000"; fi` `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH --time=1:00:00" >> ${dir2}/${samplename}_${job}.sbatch
		if [ -n "${SLURMemail}" ]
		then
			echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
		if [ -n "${SLURMpartition}" ]
		then
			echo "#SBATCH --partition=${SLURMpartition}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
		if [ -n "${SLURMqos}" ]
		then
			echo "#SBATCH --qos=${SLURMqos}" >> ${dir2}/${samplename}_${job}.sbatch
		fi

		# Require previous job successful completion
		if [ -n "${SBtrim}" ]
		then
				echo "#SBATCH --dependency=afterok:${SBalign##* }" >> ${dir2}/${samplename}_${job}.sbatch
		fi

		# General commands
		echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
		
		# Job specific commands
		echo "fastqc -o ${dir2}/ --noextract ${dir2}/${trimout1} ${dir2}/${trimout2} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch
		
		# Cleaning commands
		# remove trim.fastq files from destination folder
		if [ ! -z ${trim} ] && [ ${trim} -eq "1" ]
		then
			echo "rm ${dir2}/${trimout1} ${dir2}/${trimout2}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
		# remove .sbatch
		echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
		
		# Queue job
		SBfqc=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
		echo -e "\t FastQC job queued"
	fi



	######################
	## Prepare and queue duplicates marking job
	job="dup"

	# Picard variables
	dupout="`basename ${bamsortedout} .bam`.nodup.bam"
	dupbai="`basename ${bamsortedout} .bam`.nodup.bai"
	dupmetrics="`basename ${read1} _R1${fileext}`.dupmetrics"
	dupstdout="`basename ${read1} _R1${fileext}`.duplog"


	# General SLURM parameters
	echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH `if [ ! -z ${mem} ] && [ ${mem} -gt "64" ]; then echo "--mem=64000"; else echo "--mem=${mem}000"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --time=4:00:00" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${SLURMemail}" ]
	then
		echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMaccount}" ]
	then
		echo "#SBATCH --account=${SLURMaccount}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMpartition}" ]
	then
		echo "#SBATCH --partition=${SLURMpartition}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMqos}" ]
	then
		echo "#SBATCH --qos=${SLURMqos}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# Require previous job successful completion
	echo "#SBATCH --dependency=afterok:${SBalign##* }" >> ${dir2}/${samplename}_${job}.sbatch
	
	# General commands
	echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# Job specific commands
	# Use picard tools MarkDuplicates with removal of duplicates and index creation options.
	echo "java -Xmx`if [ ! -z ${mem} ] && [ ${mem} -gt "64" ]; then echo "64"; else echo "${mem}"; fi`g -Djava.io.tmpdir=${tmp} -jar ${picard} MarkDuplicates I=${dir2}/${bamsortedout} O=${dir2}/${dupout} METRICS_FILE=${dir2}/${logs}/${dupmetrics} REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true 2>${dir2}/${logs}/${dupstdout} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch

	# Cleaning commands
	# remove sorted files
	echo "rm ${dir2}/${bamsortedout}" >> ${dir2}/${samplename}_${job}.sbatch
	echo "rm ${dir2}/${bamsortedout}.bai" >> ${dir2}/${samplename}_${job}.sbatch
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBdup=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	echo -e "\t Duplicates marking job queued"



	## Start GATK best practice (June 2016) sample processing
	
	
	
	######################
	## Local realignment around indels
	# Note that indel realignment is no longer necessary for variant discovery if you plan to use a variant
	# caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2. However it is still
	# required when using legacy callers such as UnifiedGenotyper or the original MuTect.
	if [ ! -z ${realign} ] && [ ${realign} -eq "1" ]
	then
		job="realign"

		# GATK variables
		intervals="${samplename}.intervals"
		realigned="`basename ${dupout} .bam`.realigned.bam"
		realignedbai="`basename ${realigned} .bam`.bai"
		matefixed="`basename ${realigned} .bam`.fixed.bam"
		matefixedbai="`basename ${matefixed} .bam`.bai"

		# General SLURM parameters
		echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}" >> ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH --time=8:00:00" >> ${dir2}/${samplename}_${job}.sbatch
		if [ -n "${SLURMemail}" ]
		then
			echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
		if [ -n "${SLURMaccount}" ]
		then
			echo "#SBATCH --account=${SLURMaccount}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
		if [ -n "${SLURMpartition}" ]
		then
			echo "#SBATCH --partition=${SLURMpartition}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
		if [ -n "${SLURMqos}" ]
		then
			echo "#SBATCH --qos=${SLURMqos}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
	
		# Require previous job successful completion
		echo "#SBATCH --dependency=afterok:${SBdup##* }" >> ${dir2}/${samplename}_${job}.sbatch
	
		# General commands
		echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
	
		# Job specific commands
		# Determining (small) suspicious intervals which are likely in need of realignment (technically not required on GATK > 3.6)
		echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T RealignerTargetCreator -R ${fasta_refgenome} -I ${dir2}/${dupout} -o ${dir2}/${intervals} -known ${millsgold} -known ${onekGph1} -L ${regions} -nt ${threads} `if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch
		# Running the realigner over those intervals (technically not required on GATK > 3.6)
		echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T IndelRealigner -R ${fasta_refgenome} -I ${dir2}/${dupout} -o ${dir2}/${realigned} -targetIntervals ${dir2}/${intervals} `if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch
		# When using paired end data, the mate information must be fixed, as alignments may change during the realignment process
		echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${picard} FixMateInformation I=${dir2}/${realigned} O=${dir2}/${matefixed} SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch

		# Cleaning commands
		# remove intermediate files
		echo "rm ${dir2}/${intervals}" >> ${dir2}/${samplename}_${job}.sbatch
		echo "rm ${dir2}/${realigned}" >> ${dir2}/${samplename}_${job}.sbatch
		echo "rm ${dir2}/${realigned}bai" >> ${dir2}/${samplename}_${job}.sbatch
		# remove .sbatch
		echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

		# Queue job
		SBrealign=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
		echo -e "\t Realignment job queued"
	
		# Update variable names for next step
		dupout="${matefixed}"
		SBdup="${SBrealign}"
	fi

	######################
	## Quality score recalibration
	job="recal"

	# GATK variables
	recal_data="${samplename}.recal.grp"
	recal="`basename ${dupout} .bam`.recal.bam"
	recalbai="`basename ${recal} .bam`.bai"	

	# General SLURM parameters
	echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --time=8:00:00" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${SLURMemail}" ]
	then
		echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMaccount}" ]
	then
		echo "#SBATCH --account=${SLURMaccount}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMpartition}" ]
	then
		echo "#SBATCH --partition=${SLURMpartition}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMqos}" ]
	then
		echo "#SBATCH --qos=${SLURMqos}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# Require previous job successful completion
	echo "#SBATCH --dependency=afterok:${SBdup##* }" >> ${dir2}/${samplename}_${job}.sbatch
	
	# General commands
	echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# Job specific commands
	echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T BaseRecalibrator -R ${fasta_refgenome} -I ${dir2}/${dupout} -knownSites ${dbSNP} -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${dir2}/${recal_data} -L ${regions} -nct ${threads} `if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch
	echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T PrintReads -R ${fasta_refgenome} -I ${dir2}/${dupout} -BQSR ${dir2}/${recal_data} -o ${dir2}/${recal} -nct ${threads} `if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch

	# Cleaning commands
	# remove deduplicated files
	echo "rm ${dir2}/${dupout}" >> ${dir2}/${samplename}_${job}.sbatch
	echo "rm ${dir2}/${dupbai}" >> ${dir2}/${samplename}_${job}.sbatch
	# remove intermediate files
	echo "rm ${dir2}/${recal_data}" >> ${dir2}/${samplename}_${job}.sbatch
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBrecal=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	echo -e "\t Quality score recalibration job queued"



	######################
	## Calling SNPs
	job="gvcf"

	# GATK variables
	rawgvcf="${samplename}.raw.g.vcf"

	# General SLURM parameters
	echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --time=8:00:00" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${SLURMemail}" ]
	then
		echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMaccount}" ]
	then
		echo "#SBATCH --account=${SLURMaccount}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMpartition}" ]
	then
		echo "#SBATCH --partition=${SLURMpartition}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMqos}" ]
	then
		echo "#SBATCH --qos=${SLURMqos}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# Require previous job successful completion
	echo "#SBATCH --dependency=afterok:${SBrecal##* }" >> ${dir2}/${samplename}_${job}.sbatch
	
	# General commands
	echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
	fi

	# Job specific commands
	# Produce raw SNP calls
	echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T HaplotypeCaller --emitRefConfidence GVCF -R ${fasta_refgenome} -I ${dir2}/${recal} -D ${dbSNP} -o ${dir2}/${rawgvcf} -L ${regions} -nct ${threads} `if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch

	# Cleaning commands
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBgvcf=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	SBgvcfIDs=${SBgvcfIDs}:${SBgvcf##* }
	rawgvcfs="${rawgvcfs} ${dir2}/${rawgvcf}"
	echo -e "\t GVCF creation job queued"



	######################
	## Compute coverage
	if [ ! -z ${coverage} ] && [ ${coverage} -eq "1" ]
	then
		job="coverage"

		# GATK variables
		coverageout="${samplename}.coverage"

		# General SLURM parameters
		echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}" >> ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH --time=8:00:00" >> ${dir2}/${samplename}_${job}.sbatch
		if [ -n "${SLURMemail}" ]
		then
			echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
		if [ -n "${SLURMaccount}" ]
		then
			echo "#SBATCH --account=${SLURMaccount}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
		if [ -n "${SLURMpartition}" ]
		then
			echo "#SBATCH --partition=${SLURMpartition}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
		if [ -n "${SLURMqos}" ]
		then
			echo "#SBATCH --qos=${SLURMqos}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
	
		# Require previous job successful completion
		echo "#SBATCH --dependency=afterok:${SBrecal##* }" >> ${dir2}/${samplename}_${job}.sbatch
	
		# General commands
		echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
		fi

		# Job specific commands
		# Will compute all coverage informations needed
		echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T DepthOfCoverage -R ${fasta_refgenome} -I ${dir2}/${recal} -o ${dir2}/${coverageout} -L ${regions} -ct 10 -ct 15 -ct 30 `if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch

		# Cleaning commands
		# Remove indermediate files
		echo "rm ${dir2}/${coverageout}" >> ${dir2}/${samplename}_${job}.sbatch
		echo "rm ${dir2}/${coverageout}.sample_cumulative_coverage_counts" >> ${dir2}/${samplename}_${job}.sbatch
		echo "rm ${dir2}/${coverageout}.sample_cumulative_coverage_proportions" >> ${dir2}/${samplename}_${job}.sbatch
		echo "rm ${dir2}/${coverageout}.sample_interval_statistics" >> ${dir2}/${samplename}_${job}.sbatch
		echo "rm ${dir2}/${coverageout}.sample_interval_summary" >> ${dir2}/${samplename}_${job}.sbatch
		echo "rm ${dir2}/${coverageout}.sample_statistics" >> ${dir2}/${samplename}_${job}.sbatch
		# remove .sbatch
		echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

		# Queue job
		SBcoverage=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
		SBcoverageIDs=${SBcoverageIDs}:${SBcoverage##* }
		echo -e "\t Coverage job queued"
	fi
	
done < ${dir2}/Fastqs



# We are no longer processing sample files, but batch of GVCFs
samplename="GVCFs"



echo ""
echo "-- All sample jobs queued! --"
echo ""

	
	
######################
## GVCFs ready notification
if [ ${gvcf} -eq "1" ]
then
	job="notif"

	# General SLURM parameters
	echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH `if [ ! -z ${mem} ] && [ ${mem} -gt "1" ]; then echo "--mem=1000"; else echo "--mem=${mem}000"; fi` `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --time=2:00" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${SLURMemail}" ]
	then
	echo "#SBATCH --mail-type=END --mail-user=${SLURMemail}" >> ${dir2}/${samplename}_${job}.sbatch
	fi	
	if [ -n "${SLURMaccount}" ]
	then
		echo "#SBATCH --account=${SLURMaccount}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMpartition}" ]
	then
		echo "#SBATCH --partition=${SLURMpartition}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMqos}" ]
	then
		echo "#SBATCH --qos=${SLURMqos}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# Require previous job successful completion
	echo "#SBATCH --dependency=afterok${SBsnpscallIDs}`if [ -n \"${SBcoverageIDs}\" ]; then echo \"${SBcoverageIDs}\"; fi`" >> ${dir2}/${samplename}_${job}.sbatch	  

	# General commands
	echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
	fi

	# Job specific commands
	if [ -n "${iftttkey}" ]
	then
		# Trigger IFTTT maker channel event when it's ready, nice isn't it?
		echo "curl –silent -X POST https://maker.ifttt.com/trigger/${iftttevent}/with/key/${iftttkey}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	# Cleaning commands
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBnotif=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	echo "-- GVCFs ready notification job queued --"
	echo ""

	# Clean temporary files
	[ -f ${dir2}/files1 ] && rm ${dir2}/files1
	[ -f ${dir2}/files2 ] && rm ${dir2}/files2
	[ -f ${dir2}/Fastqs ] && rm ${dir2}/Fastqs
	
	# That's all folks!
	exit
fi



# Can be launched right away -> wait for it before raw calling
if [ -n "${popfolder}" ]
then
	job="popvcf"
	
	if [ ${popfolder: -1} == "/" ]
	then
		popfolder=${popfolder%?}
	fi
	 	
	allgvcf=`find ${popfolder} -name '*.g.vcf' | sed ':a;N;$!ba;s/\n/ -V /g'`

	# General SLURM parameters
	echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --time=8:00:00" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${SLURMemail}" ]
	then
		echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMaccount}" ]
	then
		echo "#SBATCH --account=${SLURMaccount}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMpartition}" ]
	then
		echo "#SBATCH --partition=${SLURMpartition}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ -n "${SLURMqos}" ]
	then
		echo "#SBATCH --qos=${SLURMqos}" >> ${dir2}/${samplename}_${job}.sbatch
	fi

	# General commands
	echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
	fi

	# Job specific commands
	# Produce raw SNP calls from all .g.vcf files
	echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T GenotypeGVCFs -R ${fasta_refgenome} -D ${dbSNP} -V `find ${popfolder} -name '*.g.vcf' | sed ':a;N;$!ba;s/\n/ -V /g'` -o ${popfolder}/aggregate.vcf -nt ${threads} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch
	
	# Cleaning commands
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBpopvcf=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	echo "-- Creation of .vcf file for VQSR agreggation using samples in ${popfolder}/ job queued --"

	# Update variable names for next step   
	popgvcf="${popfolder}/aggregate.vcf"
fi


######################
## Joint genotyping on gVCF files
job="jointgeno"

# GATK variables
rawgvcfs=`echo ${rawgvcfs} | sed 's/ / -V /g'`
rawSNP="${LB}-JointGenotypeOutput.vcf"
recalSNP="${LB}-JointGenotypeOutput.recalibratedsnps.vcf"
filteredSNP="${LB}-JointGenotypeOutput.filtered.vcf"
	
# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=8:00:00" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${SLURMemail}" ]
then
	echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}" >> ${dir2}/${samplename}_${job}.sbatch
fi
if [ -n "${SLURMaccount}" ]
then
	echo "#SBATCH --account=${SLURMaccount}" >> ${dir2}/${samplename}_${job}.sbatch
fi
if [ -n "${SLURMpartition}" ]
then
	echo "#SBATCH --partition=${SLURMpartition}" >> ${dir2}/${samplename}_${job}.sbatch
fi
if [ -n "${SLURMqos}" ]
then
	echo "#SBATCH --qos=${SLURMqos}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Require previous job successful completion
echo "#SBATCH --dependency=afterok${SBgvcfIDs}`if [ -n "${popfolder}" ]; then echo :${SBpopvcf##* }; fi`" >> ${dir2}/${samplename}_${job}.sbatch

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Job specific commands
# Produce raw SNP calls from all .g.vcf files
echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T GenotypeGVCFs -R ${fasta_refgenome} -D ${dbSNP} -V ${rawgvcfs} -o ${dir2}/${rawSNP} -nt ${threads} `if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch
# Variants recalibration
# Pass #1 for SNPs
echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T VariantRecalibrator -R ${fasta_refgenome} -input ${dir2}/${rawSNP} `if [ -n "${popgvcf}" ]; then echo "-aggregate ${popgvcf}"; fi` -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} -resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${onekGph1} -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 ${dbSNP} -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 90.0 -recalFile ${dir2}/recalibrate_SNP.recal -tranchesFile ${dir2}/recalibrate_SNP.tranches -rscriptFile ${dir2}/recalibrate_SNP_plots.R -nt ${threads} `if [ -n "${ped}" ]; then echo "-ped ${ped} -pedValidationType SILENT"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch
# Pass #2 for SNPs ApplyRecalibration
echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T ApplyRecalibration -R ${fasta_refgenome} -input ${dir2}/${rawSNP} -tranchesFile ${dir2}/recalibrate_SNP.tranches -recalFile ${dir2}/recalibrate_SNP.recal -o ${dir2}/${recalSNP} --ts_filter_level 99.5 -mode SNP -nt ${threads} `if [ -n "${ped}" ]; then echo "-ped ${ped} -pedValidationType SILENT"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch
# Pass #3 for Indels
echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T VariantRecalibrator -R ${fasta_refgenome} -input ${dir2}/${recalSNP} `if [ -n "${popgvcf}" ]; then echo "-aggregate ${popgvcf}"; fi` --maxGaussians 4 -resource:mills,known=false,training=true,truth=true,prior=12.0 ${millsgold} -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbSNP} -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 90.0 -recalFile ${dir2}/recalibrate_INDEL.recal -tranchesFile ${dir2}/recalibrate_INDEL.tranches -rscriptFile ${dir2}/recalibrate_INDEL_plots.R -nt ${threads} `if [ -n "${ped}" ]; then echo "-ped ${ped} -pedValidationType SILENT"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch
# Pass #4 for Indels ApplyRecalibration
echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T ApplyRecalibration -R ${fasta_refgenome} -input ${dir2}/${recalSNP} -tranchesFile ${dir2}/recalibrate_INDEL.tranches -recalFile ${dir2}/recalibrate_INDEL.recal -o ${dir2}/$filteredSNP --ts_filter_level 99.0 -mode INDEL -nt ${threads} `if [ -n "${ped}" ]; then echo "-ped ${ped} -pedValidationType SILENT"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch
	
# Cleaning commands
# Remove indermediate files
echo "rm ${dir2}/${recalSNP}" >> ${dir2}/${samplename}_${job}.sbatch
echo "rm ${dir2}/${recalSNP}.idx" >> ${dir2}/${samplename}_${job}.sbatch
echo "rm ${dir2}/recalibrate_SNP.tranches" >> ${dir2}/${samplename}_${job}.sbatch
echo "rm ${dir2}/recalibrate_SNP.recal" >> ${dir2}/${samplename}_${job}.sbatch
echo "rm ${dir2}/recalibrate_SNP.recal.idx" >> ${dir2}/${samplename}_${job}.sbatch
echo "rm ${dir2}/recalibrate_INDEL.tranches" >> ${dir2}/${samplename}_${job}.sbatch
echo "rm ${dir2}/recalibrate_INDEL.recal" >> ${dir2}/${samplename}_${job}.sbatch
echo "rm ${dir2}/recalibrate_INDEL.recal.idx" >> ${dir2}/${samplename}_${job}.sbatch
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBjointgeno=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo "-- Joint genotyping job queued --"



# We are no longer processing GVCFs
samplename="SNPs"



######################
## Annotate SNPs using ANNOVAR
job="annovar"

# ANNOVAR variables
annovarfile="${LB}-JointGenotypeOutput.annovar"
snpssummary="${LB}-JointGenotypeOutput.snps"

# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=4:00:00" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${SLURMemail}" ]
then
	echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}" >> ${dir2}/${samplename}_${job}.sbatch
fi
if [ -n "${SLURMaccount}" ]
then
	echo "#SBATCH --account=${SLURMaccount}" >> ${dir2}/${samplename}_${job}.sbatch
fi
if [ -n "${SLURMpartition}" ]
then
	echo "#SBATCH --partition=${SLURMpartition}" >> ${dir2}/${samplename}_${job}.sbatch
fi
if [ -n "${SLURMqos}" ]
then
	echo "#SBATCH --qos=${SLURMqos}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Require previous job successful completion
echo "#SBATCH --dependency=afterok:${SBjointgeno##* }" >> ${dir2}/${samplename}_${job}.sbatch

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Job specific commands
# Convert to ANNOVAR format from GATK .vcf file
echo "${annovar}/convert2annovar.pl --format vcf4 -allsample -withfreq --includeinfo ${dir2}/${filteredSNP} --outfile ${dir2}/${annovarfile} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch

# Annotate using ANNOVAR
#echo "${annovar}/table_annovar.pl --buildver hg19 ${dir2}/${annovarfile} ${annovar}/humandb/ --protocol refGene,phastConsElements46way,genomicSuperDups,gwasCatalog,esp6500siv2_all,1000g2015aug_all,exac03,gme,kaviar_20150923,avsnp147,dbnsfp33a,dbnsfp31a_interpro,dbscsnv11,clinvar_20170130,intervar_20170202,hrcr1,revel,mcap --operation g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f --otherinfo --outfile ${dir2}/${snpssummary} --remove -thread ${threads}" >> ${dir2}/${samplename}_${job}.sbatch
echo "${annovar}/table_annovar.pl --buildver ${buildver} ${dir2}/${annovarfile} ${annovar}/humandb/ --protocol ${protocol} --operation ${operation} --otherinfo --outfile ${dir2}/${snpssummary} --remove -thread ${threads} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then exit 1; else touch ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID}; fi" >> ${dir2}/${samplename}_${job}.sbatch

# Fixing headers and cleaning files
echo "sed -i \"1s/Otherinfo/Otherinfo\t\t\t\`cat ${dir2}/${filteredSNP} | grep CHROM | sed 's/#//g'\`/g\" ${dir2}/${snpssummary}.hg19_multianno.txt" >> ${dir2}/${samplename}_${job}.sbatch		# Fixing headers to add back sample names in annovar csv output files
echo "sed -i 's/,/;/g' ${dir2}/${snpssummary}.hg19_multianno.txt" >> ${dir2}/${samplename}_${job}.sbatch		# Remove any potential existing commas and replace them by semi columns
echo "sed -i 's/\\x2c//g' ${dir2}/${snpssummary}.hg19_multianno.txt" >> ${dir2}/${samplename}_${job}.sbatch		# ClinVar database is full of \x2c (comma in hexadecimal), clean it up!
echo "sed -i 's/\t/,/g' ${dir2}/${snpssummary}.hg19_multianno.txt" >> ${dir2}/${samplename}_${job}.sbatch		# Convert tabs to commas
echo "mv ${dir2}/${snpssummary}.hg19_multianno.txt ${dir2}/${snpssummary}.hg19_multianno.csv" >> ${dir2}/${samplename}_${job}.sbatch

# Cleaning commands
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBannovar=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo "-- Annovar job queued --"



######################
## Final processing of result files
job="notif"
samplename="`basename ${dir2}`-WES"

# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH `if [ ! -z ${mem} ] && [ ${mem} -gt "4" ]; then echo "--mem=4000"; else echo "--mem=${mem}000"; fi` `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=1:00:00" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${SLURMemail}" ]
then
echo "#SBATCH --mail-type=FAIL,END --mail-user=${SLURMemail}" >> ${dir2}/${samplename}_${job}.sbatch
fi	
if [ -n "${SLURMaccount}" ]
then
	echo "#SBATCH --account=${SLURMaccount}" >> ${dir2}/${samplename}_${job}.sbatch
fi
if [ -n "${SLURMpartition}" ]
then
	echo "#SBATCH --partition=${SLURMpartition}" >> ${dir2}/${samplename}_${job}.sbatch
fi
if [ -n "${SLURMqos}" ]
then
	echo "#SBATCH --qos=${SLURMqos}" >> ${dir2}/${samplename}_${job}.sbatch
fi
	
# Require previous job successful completion
echo "#SBATCH --dependency=afterok:${SBannovar##* }`if [ -n \"${SBcoverageIDs}\" ]; then echo \"${SBcoverageIDs}\"; fi`" >> ${dir2}/${samplename}_${job}.sbatch

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Job specific commands
# Creating folders for archive
echo "[ -f ${dir2}/Results.tar.gz ] && rm ${dir2}/Results.tar.gz" >> ${dir2}/${samplename}_${job}.sbatch
echo "[ -d ${dir2}/Results ] && rm -rf ${dir2}/Results" >> ${dir2}/${samplename}_${job}.sbatch
echo "mkdir -p ${dir2}/Results/" >> ${dir2}/${samplename}_${job}.sbatch
# List all csv, html, pdf and zip files
echo "csvarchive=\`ls ${dir2}/*.csv\`" >> ${dir2}/${samplename}_${job}.sbatch
echo "pdfarchive=\`ls ${dir2}/*.pdf\`" >> ${dir2}/${samplename}_${job}.sbatch
echo "htmlarchive=\`ls ${dir2}/*.html\`" >> ${dir2}/${samplename}_${job}.sbatch
echo "ziparchive=\`ls ${dir2}/*.zip\`" >> ${dir2}/${samplename}_${job}.sbatch
# Create an archive of all csv, html, pdf and zip files for easy download
echo "for i in \${csvarchive}; do cp -rf \$i ${dir2}/Results/\`basename \$i\`; done" >> ${dir2}/${samplename}_${job}.sbatch
echo "for i in \${htmlarchive}; do cp -rf \$i ${dir2}/Results/\`basename \$i\`; done" >> ${dir2}/${samplename}_${job}.sbatch
echo "for i in \${pdfarchive}; do cp -rf \$i ${dir2}/Results/\`basename \$i\`; done" >> ${dir2}/${samplename}_${job}.sbatch
echo "for i in \${ziparchive}; do cp -rf \$i ${dir2}/Results/\`basename \$i\`; done" >> ${dir2}/${samplename}_${job}.sbatch
echo "tar --remove-files -C ${dir2} -pczf ${dir2}/Results.tar.gz Results" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${iftttkey}" ] && [ -n "${iftttevent}" ]
then
	# Trigger IFTTT maker channel event when it's ready, nice isn't it?
	echo "curl -X POST -H \"Content-Type: application/json\" -d '{ \"value1\" : \"`basename ${dir2}`\" }' https://maker.ifttt.com/trigger/${iftttevent}/with/key/${iftttkey}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Cleaning commands
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBnotif=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo "-- Results ready notification job queued --"
echo ""

# Clean temporary files
[ -f ${dir2}/files1 ] && rm ${dir2}/files1
[ -f ${dir2}/files2 ] && rm ${dir2}/files2
[ -f ${dir2}/Fastqs ] && rm ${dir2}/Fastqs

# That's all folks!
exit
