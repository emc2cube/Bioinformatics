#!/bin/bash
#
# Usage: sh_WES.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]
#
##############################################################
##                       Description                        ##
##############################################################
#
# This script will process fastq(.gz) files and align them to
# a reference genome using bowtie2. It will then use Picard
# and GATK following the June 2016 best practices workflow.
# SNPs will then be annotated using ANNOVAR.
#
## Options:
#
# --help : Display help message.
# --version : Display version number.
#
##############################################################
##                  Configurable variables                  ##
##############################################################
#
#
## General options
#
# Maximum number of threads (or CPUs) to request and allocate to programs.
# In some case less than this value may automatically be allowed.
threads=$(nproc --all --ignore=1)
#
# Maximum amount of memory (in GB) to request and allocate to programs.
# In some case less than this value may automatically be allowed.
mem="128"
#
# Log folder, will be created in your destination folder.
logs="logs"
#
# Advanced: Path to temporary folder. Useful if your cluster system have a fast local I/O disk.
# Leave empty to use the destination folder. In all cases temporary files will be purged at the end of the script.
tmp="\${L_SCRATCH_JOB}"
#
# Debug mode.
# For troubleshooting only. Will keep all intermediate files,
# sbatch files and logs.
# 0 = No ; 1 = Yes
debug="0"
#
#
## SLURM options
#
# email to use to receive SLURM notifications.
# By default only send FAIL notifications for all jobs and FAIL,END for the last one.
SLURMemail=""
#
# SLURM account, partition and qos settings if required.
SLURMaccount=""
SLURMpartition=""
SLURMqos=""
#
# Custom commands to run before any other program.
# Use it to load modules for example:
# customcmd="module load R"
customcmd=""
#
#
## FastQ options
#
# If sample files are split into multiple FastQ files, merge them.
# 0 = No ; 1 = Yes
merge="0"
#
# Process undetermined/unmatched files (tag not properly recognized).
# 0 = No ; 1 = Yes
underdet="0"
#
# Process BLANK files (Sample name should be "BLANK").
# 0 = No ; 1 = Yes
blank="0"
#
# Perform fastqc on fastq files.
# If sequences are trimmed, fastqc will be performed on trimmed sequences.
# 0 = No ; 1 = Yes
fastqc="0"
#
# Trim sequences with Trimmomatic.
# 0 = No ; 1 = Yes
trim="0"
#
# Trimmomatic folder, should contain a trimmomatic.jar file
Trimmomatic="/Tools/Trimmomatic"
#
# Keep unpaired trimmed sequences.
# 0 = No ; 1 = Yes
unpaired="0"
#
#
## Alignment options
#
# Filter by MapQ quality score to remove poor aligned reads.
# samtools defaut is 0 (no filtering), 10<MapQ<20 is advised.
mapq="10"
#
# bowtie2 indexed reference genome location.
bt_refgenome="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
#
# Read group parameters.
# Library: if empty LB will be set to the destination folder name.
LB=""
#
# Platform used to produce the reads.
# Valid values: CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT and PACBIO.
PL="ILLUMINA"
#
#
## GATK options
#
# Picard-tools picard.jar location.
picard="/Tools/Picard/picard.jar"
#
# GATK GenomeAnalysisTK.jar file location.
gatk="/Tools/GATK/GenomeAnalysisTK.jar"
#
# Reference genome (fasta) file.
# It must be the same one that was indexed by bowtie2 for alignment
fasta_refgenome="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
#
# List of the targeted intervals sequenced. 
# This is necessary as only about 60-70% of all the reads will end up in exonic regions and the rest may align anywhere else in the genome.
# To restrict the output to exonic sequences, generate a file containing for example all the exons plus 50bp at each end for getting splice site information as well.
# This can be done using the UCSC Table Browser (http://genome.ucsc.edu/cgi-bin/hgTables?command=start).
# Choose the hg19 assembly of the human genome and set the track to RefSeq genes and the table to refGene.
# Use BED format as output format and assign the file an appropriate name.
# By clicking on get output, several more options can be made: Choose "create one bed record per Exon plus 50bp at each end" and save the file.
# A .bed file specific to the library preparation kit targets can be used instead if available (ex: https://earray.chem.agilent.com/suredesign/search.htm ).
regions="/Tools/AgilentBedFiles/SureSelect_Human_All_Exon_V6+UTR_r2/S07604624_Covered.bed"
#
# VQSR reference files (from GATK bundle).
# HapMap vcf file location (SNP).
hapmap="/Tools/GATK/hapmap_3.3.hg19.sites.vcf"
# 1000 genome omni vcf file location (SNP).
omni="/Tools/GATK/1000G_omni2.5.hg19.sites.vcf"
# 1000 genomes phase 1 indels vcf file location.
onekGph1="/Tools/GATK/1000G_phase1.snps.high_confidence.hg19.sites.vcf"
# dbSNP known SNPs .vcf file location.
dbSNP="/Tools/GATK/dbsnp_138.hg19.vcf"
# Mills and 1000 genomes gold standard vcf file location.
millsgold="/Tools/GATK/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
#
# Indel realignment is generally no longer necessary for variant discovery,
# unless in specific uncommon use cases.
# 0 = No ; 1 = Yes
realign="0"
#
# Stop after .g.vcf creation (no filtering, no annotation).
# Yes = 1 ; No = 0
gvcf="0"
#
# PED file to specify family relations if available.
ped=$([ -f "${dir}/$(basename ${dir}).ped" ] && echo "${dir}/$(basename ${dir}).ped")
#
# Folder containing gVCFs to be used as a learning control for VQSR calibration.
# Leave empty if you have a ready to use .vcf file and use popgvcf below.
popfolder=""
#
# VCF file to be used as a learning control for VQSR calibration.
popgvcf=""
#
# Compute coverage (slow).
# Yes = 1 ; No = 0
coverage="1"
#
#
## Annovar options
#
# Annovar folder (containing humandb folder).
annovar="/Tools/annovar"
#
# Genome build version.
# Name of the genome assembly used. Can also be determined automatically from your directory structure (see examples).
# examples:
# buildver="hg19"
# buildver=$(echo ${fasta_refgenome} | sed 's/\/Sequence.*//g' | awk -F/ '{print $NF}')
buildver=$(echo ${fasta_refgenome} | sed 's/\/Sequence.*//g' | awk -F/ '{print $NF}')
#
# Protocols to be used for annotation.
protocol="refGene,phastConsElements46way,genomicSuperDups,gwasCatalog,esp6500siv2_all,1000g2015aug_all,exac03,gme,kaviar_20150923,avsnp150,dbnsfp35a,dbnsfp31a_interpro,dbscsnv11,clinvar_20200316,intervar_20180118,hrcr1,revel,mcap,gnomad211_exome"
#
# Matching operations to be used for annotation.
operation="g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f"
#
#
## IFTTT options
#
# Trigger IFTTT when script is done.
# You must register the "Maker channel" on https://ifttt.com/maker
# Copy your private key here. Leave blank to disable this function.
# iftttkey="AbCd_15CdhUIvbsFJTHGMcfgjsdHRTgcyjt" # Not a real key, you have to use your own private key.
iftttkey=""
#
# Event name used in your IFTTT recipes.
# The maker channel will look for the combination private key + event name to then trigger your recipe.
# You can create a recipe to send an email, a text message or a push notification.
iftttevent="SNPCalling"
#
#
## Setup done. You should not need to edit below this point ##

# Help!
if [ "${1}" = "--help" ] || [ "${2}" = "--help" ] || [ "${3}" = "--help" ] || [ "${4}" = "--help" ]
then
	echo "Usage: $(basename "${0}") </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]"
	echo ""
	echo "Description"
	echo ""
	echo "This script will process fastq(.gz) files and align them to"
	echo "a reference genome using bowtie2. It will then use Picard"
	echo "and GATK following the June 2016 best practices workflow."
	echo "SNPs will then be annotated using ANNOVAR."
	echo ""
	echo "Options:"
	echo "$(basename "${0}") --help : Display this help message."
	echo "$(basename "${0}") --version : Display version number."
	echo ""
	exit
fi

# Version
if [ "${1}" = "--version" ] || [ "${2}" = "--version" ] || [ "${3}" = "--version" ] || [ "${4}" = "--version" ]
then
	echo "$(basename "${0}") version 2.0.1"
	echo "Only for GATK 3.X ( 3.8.1 tested, download from https://software.broadinstitute.org/gatk/download/archive )"
	echo "Major code cleaning (2.0.1)"
	echo ".g.vcf support (2.0)"
	exit
fi

# Get fastq directory
dir="${1}"

# Get destination directory
dir2="${2}"

# Get config file location
config="${3}"

# Check paths and trailing / in directories
if [ -z "${dir}" ] || [ -z "${dir2}" ]
then
	${0} --help
	exit
fi

if [ "${dir: -1}" = "/" ]
then
	dir=${dir%?}
fi

if [ "${dir2: -1}" = "/" ]
then
	dir2=${dir2%?}
fi

if [ -n "${config}" ]
then
	if [ "${config: -4}" = ".ini" ]
	then
		# shellcheck disable=SC1090
		source "${config}"
	else
		echo "Invalid config file detected. Is it an .ini file?"
		echo ""
		${0} --help
		exit
	fi
fi

if [ -n "${popgvcf}" ]
then
	if [ "${popgvcf: -4}" != ".vcf"  ]
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

# GATK version check
if [ "$(basename ${gatk})" = "GenomeAnalysisTK.jar" ]
then
	vgatk="3"
elif [ "$(basename ${gatk})" = "gatk" ]
then
	vgatk="4"
	echo "GATK 4.X detected, pipeline not yet compatible."
	echo "Only for GATK 3.X ( 3.8.1 tested, download from https://software.broadinstitute.org/gatk/download/archive )"
	exit
fi

# Test if sequence files are .fastq or .fastq.gz
fastqgz=$(find -L "${dir}" -maxdepth 1 -name '*.fastq.gz')
fastq=$(find -L "${dir}" -maxdepth 1 -name '*.fastq')
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
	gunzip "${dir}/*.fastq.gz"
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
if [ ${mapq} = "0" ] || [ -z ${mapq} ]
then
	echo "MapQ score will not be used to filter reads"
else
	echo "Reads with a MapQ score <${mapq} will be discarded"
fi

echo ""
echo "-- GATK options"
echo ""
echo "Picard-tools are installed in ${picard}"
echo "GATK version ${vgatk} is installed in ${gatk}"
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

echo ""
echo "---------------"

# Initialize
[ -f "${dir2}/files1" ] && rm "${dir2}/files1"
[ -f "${dir2}/files2" ] && rm "${dir2}/files2"
[ -f "${dir2}/Fastqs" ] && rm "${dir2}/Fastqs"
[ -d "${dir2}/logs" ] && rm -rf "${dir2}/logs"
mkdir -p "${dir2}/"
mkdir -p "${dir2}/${logs}"

if [ -z "${tmp}" ]
then
	tmp="${dir2}/tmp"
fi

# Concatenate files if split
if [ "${merge}" = "1" ]
then
	echo ""
	echo "-- Merging files --"
	echo ""
	mkdir -p "${dir}/FastQbackup/"
	filestomergeR1=$(ls "${dir}"/*_L[0-9][0-9][0-9]_R1)
	filestomergeR2=$(ls "${dir}"/*_L[0-9][0-9][0-9]_R2)
	echo "Processing R1 files"
	echo ""
	for i in ${filestomergeR1}
	do
		sampleoutput="${i//_L[0-9][0-9][0-9]_R1/_R1}"
		echo "Processing ${i}"
		cat "${i}" >> "${sampleoutput}"
		mv "${i}" "${dir}/FastQbackup/"
	done
	echo ""
	echo "Processing R2 files"
	echo ""
	for i in ${filestomergeR2}
	do
		sampleoutput="${i//_L[0-9][0-9][0-9]_R2/_R2}"
		echo "Processing ${i}"
		cat "${i}" >> "${sampleoutput}"
		mv "${i}" "${dir}/FastQbackup/"
	done
	echo ""
	echo "-- Files merged --"
fi

if [ "${underdet}" = "0" ]
then
	if [ "${blank}" = "0" ]
	then
		# Remove all Undetermined_* and BLANK* files
		find -L "${dir}" -name '*_R1*' -not -name '*ndetermined*' -not -name '*nmatched*' -not -name 'BLANK*' | sed 's#.*/##' | sort -n > "${dir2}/files1"
		find -L "${dir}" -name '*_R2*' -not -name '*ndetermined*' -not -name '*nmatched*' -not -name 'BLANK*' | sed 's#.*/##' | sort -n > "${dir2}/files2"
		paste "${dir2}/files1" "${dir2}/files2" > "${dir2}/Fastqs"
	else
		# Remove all Undetermined_* files
		find -L "${dir}" -name '*_R1*' -not -name '*ndetermined*' -not -name '*nmatched*' | sed 's#.*/##' | sort -n > "${dir2}/files1"
		find -L "${dir}" -name '*_R2*' -not -name '*ndetermined*' -not -name '*nmatched*' | sed 's#.*/##' | sort -n > "${dir2}/files2"
		paste "${dir2}/files1" "${dir2}/files2" > "${dir2}/Fastqs"
	fi
else
	if [ "${blank}" = "0" ]
	then
		# Remove all BLANK* files
		find -L "${dir}" -name '*_R1*' -not -name 'BLANK*' | sort -n | sed 's#.*/##' > "${dir2}/files1"
		find -L "${dir}" -name '*_R2*' -not -name 'BLANK*' | sort -n | sed 's#.*/##' > "${dir2}/files2"
		paste "${dir2}/files1" "${dir2}/files2" > "${dir2}/Fastqs"
	else
		# Process all the files!
		find -L "${dir}" -name '*_R1*' | sed 's#.*/##' | sort -n > "${dir2}/files1"
		find -L "${dir}" -name '*_R2*' | sed 's#.*/##' | sort -n > "${dir2}/files2"
		paste "${dir2}/files1" "${dir2}/files2" > "${dir2}/Fastqs"
	fi
fi



######################
## Processing samples
######################

echo ""
echo "-- Queuing sample jobs --"

while read -r line;
do

	# General variables
	read1=$(echo "${line}" | cut -f1)
	read2=$(echo "${line}" | cut -f2)
	samplename=$(echo "${read1}" | awk -F_R1 '{print $1}')

	echo ""
	echo "-- Queuing sample jobs"



	######################
	## Prepare and queue trimming job
	job="trim"

	# Trimmomatic variables
	trimout1="${dir2}/$(basename "${read1}" ${fileext}).trim.fastq"
	trimout2="${dir2}/$(basename "${read2}" ${fileext}).trim.fastq"

	if [ -n "${unpaired}" ] && [ "${unpaired}" = "1" ]
	then
		unpaired1="${dir2}/$(basename "${read1}" ${fileext}).unpaired.fastq"  # Save unpaired reads
		unpaired2="${dir2}/$(basename "${read2}" ${fileext}).unpaired.fastq"  # Save unpaired reads
	else
		unpaired1="/dev/null"  # Unpaired reads are discarded
		unpaired2="/dev/null"  # Unpaired reads are discarded
	fi

	if [ -n "${trim}" ] && [ ${trim} -eq "1" ]
	then
		{
			# General SLURM parameters
			echo '#!/bin/bash'
			echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append"
			echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi) $(if [ -n "${threads}" ] && [ "${threads}" -gt "2" ]; then echo "--cpus-per-task=2"; else echo "--cpus-per-task=${threads}"; fi)"
			echo "#SBATCH --time=3:00:00"
			if [ -n "${SLURMemail}" ]
			then
				echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}"
			fi
			if [ -n "${SLURMpartition}" ]
			then
				echo "#SBATCH --partition=${SLURMpartition}"
			fi
			if [ -n "${SLURMqos}" ]
			then
				echo "#SBATCH --qos=${SLURMqos}"
			fi

			# General commands
			echo "mkdir -p ${tmp}"
			if [ -n "${customcmd}" ]
			then
				echo "${customcmd}"
			fi

			# Job specific commands
			# Trim fastq files with trimmomatic
			echo "java -Xmx$(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "8"; else echo "${mem}"; fi)g -Djava.io.tmpdir=${tmp} -jar ${Trimmomatic}/trimmomatic.jar PE -threads $(if [ -n "${threads}" ] && [ "${threads}" -gt "2" ]; then echo "2"; else echo "${threads}"; fi) -phred33 ${dir}/${read1} ${dir}/${read2} ${trimout1} ${unpaired1} ${trimout2} ${unpaired2} ILLUMINACLIP:${Trimmomatic}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

			# Cleaning commands
			if [ "${debug}" != "1" ]
			then
				# Remove .sbatch
				echo "rm ${dir2}/${samplename}_${job}.sbatch"
			fi

		} > "${dir2}/${samplename}_${job}.sbatch"

		# Queue job
		SBtrim=$(until sbatch "${dir2}/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
		echo -e "\t Trimming job queued"
	else
		# Need to rename some variables
		trimout1="${dir}/${read1}"
		trimout2="${dir}/${read2}"
	fi



	######################
	## Prepare and queue alignment job
	job="align"

	# bowtie2 variables
	samout=$(basename "${read1}" | sed "s/_R1${fileext}/.sam/g")
	bamout=$(basename "${read1}" | sed "s/_R1${fileext}/.bam/g")
	bamsortedout=$(basename "${read1}" | sed "s/_R1${fileext}/.sorted.bam/g")
	if [ -z "$LB" ]
	then
		LB=$(basename "${dir2}")
	fi
	SM=$(echo "${read1}" | awk -F_ '{print $1}')
	if [ "${fileext}" = ".fastq.gz" ]
	then
		CN=$(gzip -cd "${dir}"/"${read1}" | head -n 1 | awk -F: '{print $1}' | sed 's/@//')
		PU=$(gzip -cd "${dir}"/"${read1}" | head -n 1 | awk -F: '{print $3}')
	else
		CN=$(head -n 1 "${dir}"/"${read1}" | awk -F: '{print $1}' | sed 's/@//')
		PU=$(head -n 1 "${dir}"/"${read1}" | awk -F: '{print $3}')
	fi

	{
		# General SLURM parameters
		echo '#!/bin/bash'
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append"
		echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
		echo "#SBATCH --time=8:00:00"
		if [ -n "${SLURMemail}" ]
		then
			echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}"
		fi
		if [ -n "${SLURMaccount}" ]
		then
			echo "#SBATCH --account=${SLURMaccount}"
		fi
		if [ -n "${SLURMpartition}" ]
		then
			echo "#SBATCH --partition=${SLURMpartition}"
		fi
		if [ -n "${SLURMqos}" ]
		then
			echo "#SBATCH --qos=${SLURMqos}"
		fi

		# Require previous job successful completion
		if [ -n "${SBtrim}" ]
		then
			echo "#SBATCH --dependency=afterok:${SBtrim##* }"
		fi

		# General commands
		echo "mkdir -p ${tmp}"
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}"
		fi

		# Job specific commands
		# Perform alignment with bowtie
		echo "bowtie2 -p ${threads} --phred33 --rg-id ${LB}_${SM} --rg CN:${CN} --rg LB:${LB} --rg PL:${PL} --rg PU:${PU} --rg SM:${SM} -x ${bt_refgenome} -S ${dir2}/${samout} -1 ${trimout1} -2 ${trimout2} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
		# Convert .sam to .bam with optional filter on MapQ quality score
		echo "samtools view -bS -q ${mapq} -o ${dir2}/${bamout} ${dir2}/${samout} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
		# Sort .bam file
		echo "samtools sort -@ ${threads} -o ${dir2}/${bamsortedout} -O bam -T ${tmp} ${dir2}/${bamout} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
		# Index sorted .bam
		echo "samtools index ${dir2}/${bamsortedout} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

		# Cleaning commands
		if [ "${debug}" != "1" ]
		then
			# Remove trim.fastq files from destination folder
			if [ "${trim}" = "1" ] && [ "${fastqc}" = "0" ]
			then
				echo "rm ${trimout1} ${trimout2}"
			fi
			# Remove .sam file
			echo "rm ${dir2}/${samout}"
			# Remove unsorted .bam file
			echo "rm ${dir2}/${bamout}"
			# Remove .sbatch
			echo "rm ${dir2}/${samplename}_${job}.sbatch"
		fi

	} > "${dir2}/${samplename}_${job}.sbatch"

	# Queue job
	SBalign=$(until sbatch "${dir2}/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
	if [ -n "${SBalign##* }" ]
	then
		SBalignIDs=${SBalignIDs}:${SBalign##* }
	fi
	echo -e "\t Alignment job queued"



	######################
	## Run fastqc to generate quality control files
	if [ -n "${fastqc}" ] && [ ${fastqc} -eq "1" ]
	then
		job="fqc"

		{
			# General SLURM parameters
			echo '#!/bin/bash'
			echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append"
			echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "4" ]; then echo "--mem=4000"; else echo "--mem=${mem}000"; fi) $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi)"
			echo "#SBATCH --time=1:00:00"
			if [ -n "${SLURMemail}" ]
			then
				echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}"
			fi
			if [ -n "${SLURMpartition}" ]
			then
				echo "#SBATCH --partition=${SLURMpartition}"
			fi
			if [ -n "${SLURMqos}" ]
			then
				echo "#SBATCH --qos=${SLURMqos}"
			fi

			# Require previous job successful completion
			if [ -n "${SBalign}" ]
			then
				echo "#SBATCH --dependency=afterok:${SBalign##* }"
			fi

			# General commands
			echo "mkdir -p ${tmp}"
			if [ -n "${customcmd}" ]
			then
				echo "${customcmd}"
			fi

			# Job specific commands
			echo "fastqc -o ${dir2}/ --noextract ${dir2}/${trimout1} ${dir2}/${trimout2} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

			# Cleaning commands
			if [ "${debug}" != "1" ]
			then
				# Remove trim.fastq files from destination folder
				if [ "${trim}" = "1" ]
				then
					echo "rm ${trimout1} ${trimout2}"
				fi
				# Remove .sbatch
				echo "rm ${dir2}/${samplename}_${job}.sbatch"
			fi

		} > "${dir2}/${samplename}_${job}.sbatch"

		# Queue job
		SBfqc=$(until sbatch "${dir2}/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
		if [ -n "${SBfqc##* }" ]
		then
			SBfqcIDs=${SBfqcIDs}:${SBfqc##* }
		fi
		echo -e "\t FastQC job queued"
	fi



	######################
	## Prepare and queue duplicates marking job
	job="dup"

	# Picard variables
	dupout="$(basename "${bamsortedout}" .bam).nodup.bam"
	dupbai="$(basename "${bamsortedout}" .bam).nodup.bai"
	dupmetrics="$(basename "${read1}" _R1${fileext}).dupmetrics"

	{
		# General SLURM parameters
		echo '#!/bin/bash'
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append"
		echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "64" ]; then echo "--mem=64000"; else echo "--mem=${mem}000"; fi)"
		echo "#SBATCH --time=4:00:00"
		if [ -n "${SLURMemail}" ]
		then
			echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}"
		fi
		if [ -n "${SLURMaccount}" ]
		then
			echo "#SBATCH --account=${SLURMaccount}"
		fi
		if [ -n "${SLURMpartition}" ]
		then
			echo "#SBATCH --partition=${SLURMpartition}"
		fi
		if [ -n "${SLURMqos}" ]
		then
			echo "#SBATCH --qos=${SLURMqos}"
		fi

		# Require previous job successful completion
		echo "#SBATCH --dependency=afterok:${SBalign##* }"

		# General commands
		echo "mkdir -p ${tmp}"
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}"
		fi

		# Job specific commands
		# Use picard tools MarkDuplicates with removal of duplicates and index creation options.
		echo "java -Xmx$(if [ -n "${mem}" ] && [ ${mem} -gt "64" ]; then echo "64"; else echo "${mem}"; fi)g -Djava.io.tmpdir=${tmp} -jar ${picard} MarkDuplicates I=${dir2}/${bamsortedout} O=${dir2}/${dupout} METRICS_FILE=${dir2}/${logs}/${dupmetrics} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true SORTING_COLLECTION_SIZE_RATIO=0.20 || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

		# Cleaning commands
		if [ "${debug}" != "1" ]
		then
			# Remove sorted files
			echo "rm ${dir2}/${bamsortedout}"
			echo "rm ${dir2}/${bamsortedout}.bai"
			# Remove .sbatch
			echo "rm ${dir2}/${samplename}_${job}.sbatch"
		fi

	} > "${dir2}/${samplename}_${job}.sbatch"

	# Queue job
	SBdup=$(until sbatch "${dir2}/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
	if [ -n "${SBdup##* }" ]
	then
		SBdupIDs=${SBdupIDs}:${SBdup##* }
	fi
	echo -e "\t Duplicates marking job queued"



	## Start GATK best practice (June 2016) Somatic short variant discovery (SNVs + Indels)



	######################
	## Local realignment around indels
	# Note that indel realignment is no longer necessary for variant discovery if you plan to use a variant
	# caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2. However it is still
	# required when using legacy callers such as UnifiedGenotyper or the original MuTect.
	if [ -n "${realign}" ] && [ ${realign} -eq "1" ]
	then
		job="realign"

		# GATK variables
		intervals="${samplename}.intervals"
		realigned="$(basename "${dupout}" .bam).realigned.bam"
		matefixed="$(basename "${realigned}" .bam).fixed.bam"

		{
			# General SLURM parameters
			echo '#!/bin/bash'
			echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append"
			echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
			echo "#SBATCH --time=8:00:00"
			if [ -n "${SLURMemail}" ]
			then
				echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}"
			fi
			if [ -n "${SLURMaccount}" ]
			then
				echo "#SBATCH --account=${SLURMaccount}"
			fi
			if [ -n "${SLURMpartition}" ]
			then
				echo "#SBATCH --partition=${SLURMpartition}"
			fi
			if [ -n "${SLURMqos}" ]
			then
				echo "#SBATCH --qos=${SLURMqos}"
			fi

			# Require previous job successful completion
			echo "#SBATCH --dependency=afterok:${SBdup##* }"

			# General commands
			echo "mkdir -p ${tmp}"
			if [ -n "${customcmd}" ]
			then
				echo "${customcmd}"
			fi

			# Job specific commands
			if [ "${vgatk}" = 3 ]
			then
				# Determining (small) suspicious intervals which are likely in need of realignment (technically not required on GATK > 3.6)
				echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T RealignerTargetCreator -R ${fasta_refgenome} -I ${dir2}/${dupout} -o ${dir2}/${intervals} -known ${millsgold} -known ${onekGph1} -L ${regions} -nt ${threads} $(if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi) || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
				# Running the realigner over those intervals (technically not required on GATK > 3.6)
				echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T IndelRealigner -R ${fasta_refgenome} -I ${dir2}/${dupout} -o ${dir2}/${realigned} -targetIntervals ${dir2}/${intervals} $(if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi) || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
				# When using paired end data, the mate information must be fixed, as alignments may change during the realignment process
				echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${picard} FixMateInformation I=${dir2}/${realigned} O=${dir2}/${matefixed} SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
			else
				echo "GATK4 command; not implemented yet"
				exit
				"${gatk}" --java-options "-Xmx${mem}g -Djava.io.tmpdir=${tmp}"
			fi

			# Cleaning commands
			if [ "${debug}" != "1" ]
			then
				# Remove intermediate files
				echo "rm ${dir2}/${intervals}"
				echo "rm ${dir2}/${realigned}"
				echo "rm ${dir2}/${realigned}bai"
				# Remove .sbatch
				echo "rm ${dir2}/${samplename}_${job}.sbatch"
			fi

		} > "${dir2}/${samplename}_${job}.sbatch"

		# Queue job
		SBrealign=$(until sbatch "${dir2}/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
		if [ -n "${SBrealign##* }" ]
		then
			SBrealignIDs=${SBrealignIDs}:${SBrealign##* }
		fi
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
	recal="$(basename "${dupout}" .bam).recal.bam"

	{
		# General SLURM parameters
		echo '#!/bin/bash'
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append"
		echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
		echo "#SBATCH --time=8:00:00"
		if [ -n "${SLURMemail}" ]
		then
			echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}"
		fi
		if [ -n "${SLURMaccount}" ]
		then
			echo "#SBATCH --account=${SLURMaccount}"
		fi
		if [ -n "${SLURMpartition}" ]
		then
			echo "#SBATCH --partition=${SLURMpartition}"
		fi
		if [ -n "${SLURMqos}" ]
		then
			echo "#SBATCH --qos=${SLURMqos}"
		fi

		# Require previous job successful completion
		echo "#SBATCH --dependency=afterok:${SBdup##* }"

		# General commands
		echo "mkdir -p ${tmp}"
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}"
		fi

		# Job specific commands
		if [ "${vgatk}" = 3 ]
		then
			echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T BaseRecalibrator -R ${fasta_refgenome} -I ${dir2}/${dupout} -knownSites ${dbSNP} -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${dir2}/${recal_data} -L ${regions} -nct ${threads} $(if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi) || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
			#echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T PrintReads -R ${fasta_refgenome} -I ${dir2}/${dupout} -BQSR ${dir2}/${recal_data} -o ${dir2}/${recal} -nct ${threads} $(if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi) || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
			## Dirty fix for Java fatal error while using PrintReads, max memory set to 16g https://gatkforums.broadinstitute.org/gatk/discussion/10353/gatk-3-8-0-printreads-fatal-error
			echo "java -Xmx$(if [ -n "${mem}" ] && [ ${mem} -gt "16" ]; then echo "16"; else echo "${mem}"; fi)g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T PrintReads -R ${fasta_refgenome} -I ${dir2}/${dupout} -BQSR ${dir2}/${recal_data} -o ${dir2}/${recal} -nct ${threads} $(if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi) || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
		else
			echo "GATK4 command; not implemented yet"
			"${gatk}" --java-options "-Xmx${mem}g -Djava.io.tmpdir=${tmp}"
			exit
		fi

		# Cleaning commands
		if [ "${debug}" != "1" ]
		then
			# Remove deduplicated files
			echo "rm ${dir2}/${dupout}"
			echo "rm ${dir2}/${dupbai}"
			# Remove intermediate files
			echo "rm ${dir2}/${recal_data}"
			# Remove .sbatch
			echo "rm ${dir2}/${samplename}_${job}.sbatch"
		fi

	} > "${dir2}/${samplename}_${job}.sbatch"

	# Queue job
	SBrecal=$(until sbatch "${dir2}/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
	if [ -n "${SBrecal##* }" ]
	then
		SBrecalIDs=${SBrecalIDs}:${SBrecal##* }
	fi
	echo -e "\t Quality score recalibration job queued"



	######################
	## Calling SNPs
	job="gvcf"

	# GATK variables
	rawgvcf="${samplename}.raw.g.vcf"

	{
		# General SLURM parameters
		echo '#!/bin/bash'
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append"
		echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
		echo "#SBATCH --time=8:00:00"
		if [ -n "${SLURMemail}" ]
		then
			echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}"
		fi
		if [ -n "${SLURMaccount}" ]
		then
			echo "#SBATCH --account=${SLURMaccount}"
		fi
		if [ -n "${SLURMpartition}" ]
		then
			echo "#SBATCH --partition=${SLURMpartition}"
		fi
		if [ -n "${SLURMqos}" ]
		then
			echo "#SBATCH --qos=${SLURMqos}"
		fi

		# Require previous job successful completion
		echo "#SBATCH --dependency=afterok:${SBrecal##* }"

		# General commands
		echo "mkdir -p ${tmp}"
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}"
		fi

		# Job specific commands
		if [ "${vgatk}" = 3 ]
		then
			# Produce raw SNP calls
			echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T HaplotypeCaller --emitRefConfidence GVCF -R ${fasta_refgenome} -I ${dir2}/${recal} -D ${dbSNP} -o ${dir2}/${rawgvcf} -L ${regions} -nct ${threads} $(if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi) || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
		else
			echo "GATK4 command; not implemented yet"
			exit
			"${gatk}" --java-options "-Xmx${mem}g -Djava.io.tmpdir=${tmp}"
		fi

		# Cleaning commands
		if [ "${debug}" != "1" ]
		then
			# Remove .sbatch
			echo "rm ${dir2}/${samplename}_${job}.sbatch"
		fi

	} > "${dir2}/${samplename}_${job}.sbatch"

	# Queue job
	SBgvcf=$(until sbatch "${dir2}/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
	if [ -n "${SBgvcf##* }" ]
	then
		SBgvcfIDs=${SBgvcfIDs}:${SBgvcf##* }
	fi

	rawgvcfs="${dir2}/${rawgvcf}${rawgvcfs:+ ${rawgvcfs}}"
	echo -e "\t GVCF creation job queued"



	######################
	## Compute coverage
	if [ -n "${coverage}" ] && [ ${coverage} -eq "1" ]
	then
		job="coverage"

		# GATK variables
		coverageout="${samplename}.coverage"

		{
			# General SLURM parameters
			echo '#!/bin/bash'
			echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append"
			echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
			echo "#SBATCH --time=8:00:00"
			if [ -n "${SLURMemail}" ]
			then
				echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}"
			fi
			if [ -n "${SLURMaccount}" ]
			then
				echo "#SBATCH --account=${SLURMaccount}"
			fi
			if [ -n "${SLURMpartition}" ]
			then
				echo "#SBATCH --partition=${SLURMpartition}"
			fi
			if [ -n "${SLURMqos}" ]
			then
				echo "#SBATCH --qos=${SLURMqos}"
			fi

			# Require previous job successful completion
			echo "#SBATCH --dependency=afterok:${SBrecal##* }"

			# General commands
			echo "mkdir -p ${tmp}"
			if [ -n "${customcmd}" ]
			then
				echo "${customcmd}"
			fi

			# Job specific commands
			if [ "${vgatk}" = 3 ]
			then
				# Will compute all coverage informations needed
				echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T DepthOfCoverage -R ${fasta_refgenome} -I ${dir2}/${recal} -o ${dir2}/${coverageout} -L ${regions} -ct 10 -ct 15 -ct 30 $(if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi) || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
			else
				echo "GATK4 command; not implemented yet"
				exit
				"${gatk}" --java-options "-Xmx${mem}g -Djava.io.tmpdir=${tmp}"
			fi

			# Cleaning commands
			if [ "${debug}" != "1" ]
			then
				# Remove indermediate files
				echo "rm ${dir2}/${coverageout}"
				echo "rm ${dir2}/${coverageout}.sample_cumulative_coverage_counts"
				echo "rm ${dir2}/${coverageout}.sample_cumulative_coverage_proportions"
				echo "rm ${dir2}/${coverageout}.sample_interval_statistics"
				echo "rm ${dir2}/${coverageout}.sample_interval_summary"
				echo "rm ${dir2}/${coverageout}.sample_statistics"
				# Remove .sbatch
				echo "rm ${dir2}/${samplename}_${job}.sbatch"
			fi

		} > "${dir2}/${samplename}_${job}.sbatch"

		# Queue job
		SBcoverage=$(until sbatch "${dir2}/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
		if [ -n "${SBcoverage##* }" ]
		then
			SBcoverageIDs=${SBcoverageIDs}:${SBcoverage##* }
		fi

		echo -e "\t Coverage job queued"
	fi

done < "${dir2}/Fastqs"

echo ""
echo "-- All sample jobs queued! --"



######################
## GVCFs ready notification
if [ ${gvcf} -eq "1" ]
then
	job="notif"
	# We are no longer processing sample files, but batch of GVCFs
	samplename="GVCFs"

	{
		# General SLURM parameters
		echo '#!/bin/bash' > "${dir2}"/${samplename}_"${job}".sbatch
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append"
		echo "#SBATCH $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi)"
		echo "#SBATCH --time=2:00"
		if [ -n "${SLURMemail}" ]
		then
			echo "#SBATCH --mail-type=END --mail-user=${SLURMemail}"
		fi
		if [ -n "${SLURMaccount}" ]
		then
			echo "#SBATCH --account=${SLURMaccount}"
		fi
		if [ -n "${SLURMpartition}" ]
		then
			echo "#SBATCH --partition=${SLURMpartition}"
		fi
		if [ -n "${SLURMqos}" ]
		then
			echo "#SBATCH --qos=${SLURMqos}"
		fi

		# Require previous job successful completion
		echo "#SBATCH --dependency=afterok$(if [ -n "${SBfqcIDs}" ]; then echo "${SBfqcIDs}"; fi)$(if [ -n "${SBcoverageIDs}" ]; then echo "${SBcoverageIDs}"; fi)"

		# General commands
		echo "mkdir -p ${tmp}"
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}"
		fi

		# Job specific commands
		if [ -n "${iftttkey}" ]
		then
			# Trigger IFTTT maker channel event when it's ready, nice isn't it?
			echo "curl -X POST -H \"Content-Type: application/json\" -d '{ \"value1\" : \"$(basename "${dir2}")\" , \"value2\" : \"$(whoami)\"}' https://maker.ifttt.com/trigger/${iftttevent}/with/key/${iftttkey}"
		fi

		# Cleaning commands
		if [ "${debug}" != "1" ]
		then
			# Remove .sbatch
			echo "rm ${dir2}/${samplename}_${job}.sbatch"
		fi

	} > "${dir2}/${samplename}_${job}.sbatch"

	# Queue job
	until sbatch "${dir2}/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done  >/dev/null 2>&1
	echo ""
	echo "-- GVCFs ready notification job queued --"
	echo ""

	# Clean temporary files
	[ -f "${dir2}"/files1 ] && rm "${dir2}"/files1
	[ -f "${dir2}"/files2 ] && rm "${dir2}"/files2
	[ -f "${dir2}"/Fastqs ] && rm "${dir2}"/Fastqs

	# That's all folks!
	exit
fi



# Can be launched right away -> wait for it before raw calling
if [ -n "${popfolder}" ]
then
	job="popvcf"
	# We are no longer processing sample files, but batch of GVCFs
	samplename="GVCFs"

	if [ "${popfolder: -1}" = "/" ]
	then
		popfolder=${popfolder%?}
	fi

	{
		# General SLURM parameters
		echo '#!/bin/bash'
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append"
		echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
		echo "#SBATCH --time=8:00:00"
		if [ -n "${SLURMemail}" ]
		then
			echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}"
		fi
		if [ -n "${SLURMaccount}" ]
		then
			echo "#SBATCH --account=${SLURMaccount}"
		fi
		if [ -n "${SLURMpartition}" ]
		then
			echo "#SBATCH --partition=${SLURMpartition}"
		fi
		if [ -n "${SLURMqos}" ]
		then
			echo "#SBATCH --qos=${SLURMqos}"
		fi

		# General commands
		echo "mkdir -p ${tmp}"
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}"
		fi

		# Job specific commands
		if [ "${vgatk}" = 3 ]
		then
			# Produce raw SNP calls from all .g.vcf files
			echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T GenotypeGVCFs -R ${fasta_refgenome} -D ${dbSNP} -V $(find "${popfolder}" -name '*.g.vcf' | sed ':a;N;$!ba;s/\n/ -V /g') -o ${popfolder}/aggregate.vcf -nt ${threads} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
		else
			echo "GATK4 command; not implemented yet"
			exit
			"${gatk}" --java-options "-Xmx${mem}g -Djava.io.tmpdir=${tmp}"
		fi

		# Cleaning commands
		if [ "${debug}" != "1" ]
		then
			# Remove .sbatch
			echo "rm ${dir2}/${samplename}_${job}.sbatch"
		fi

	} > "${dir2}/${samplename}_${job}.sbatch"

	# Queue job
	SBpopvcf=$(until sbatch "${dir2}/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
	echo ""
	echo "-- Creation of .vcf file for VQSR agreggation using samples in ${popfolder}/ job queued --"

	# Update variable names for next step
	popgvcf="${popfolder}/aggregate.vcf"
fi


######################
## Joint genotyping on gVCF files
job="jointgeno"
# We are no longer processing sample files, but batch of GVCFs
samplename="GVCFs"

# GATK variables
rawgvcfs="${rawgvcfs// / -V }"
rawSNP="${LB}-JointGenotypeOutput.vcf"
recalSNP="${LB}-JointGenotypeOutput.recalibratedsnps.vcf"
filteredSNP="${LB}-JointGenotypeOutput.filtered.vcf"

{
	# General SLURM parameters
	echo '#!/bin/bash'
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append"
	echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
	echo "#SBATCH --time=8:00:00"
	if [ -n "${SLURMemail}" ]
	then
		echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}"
	fi
	if [ -n "${SLURMaccount}" ]
	then
		echo "#SBATCH --account=${SLURMaccount}"
	fi
	if [ -n "${SLURMpartition}" ]
	then
		echo "#SBATCH --partition=${SLURMpartition}"
	fi
	if [ -n "${SLURMqos}" ]
	then
		echo "#SBATCH --qos=${SLURMqos}"
	fi

	# Require previous job successful completion
	echo "#SBATCH --dependency=afterok${SBgvcfIDs}$(if [ -n "${popfolder}" ]; then echo :"${SBpopvcf##* }"; fi)"

	# General commands
	echo "mkdir -p ${tmp}"
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}"
	fi

	# Job specific commands
	if [ "${vgatk}" = 3 ]
	then
		# Produce raw SNP calls from all .g.vcf files
		echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T GenotypeGVCFs -R ${fasta_refgenome} -D ${dbSNP} -V ${rawgvcfs} -o ${dir2}/${rawSNP} -nt ${threads} $(if [ -n "${ped}" ]; then echo "-ped ${ped}"; fi) || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
		# Variants recalibration
		# Pass #1 for SNPs
		echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T VariantRecalibrator -R ${fasta_refgenome} -input ${dir2}/${rawSNP} $(if [ -n "${popgvcf}" ]; then echo "-aggregate ${popgvcf}"; fi) -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} -resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${onekGph1} -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 ${dbSNP} -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 90.0 -recalFile ${dir2}/recalibrate_SNP.recal -tranchesFile ${dir2}/recalibrate_SNP.tranches -rscriptFile ${dir2}/recalibrate_SNP_plots.R -nt ${threads} $(if [ -n "${ped}" ]; then echo "-ped ${ped} -pedValidationType SILENT"; fi) || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
		# Pass #2 for SNPs ApplyRecalibration
		echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T ApplyRecalibration -R ${fasta_refgenome} -input ${dir2}/${rawSNP} -tranchesFile ${dir2}/recalibrate_SNP.tranches -recalFile ${dir2}/recalibrate_SNP.recal -o ${dir2}/${recalSNP} --ts_filter_level 99.5 -mode SNP -nt ${threads} $(if [ -n "${ped}" ]; then echo "-ped ${ped} -pedValidationType SILENT"; fi) || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
		# Pass #3 for Indels
		echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T VariantRecalibrator -R ${fasta_refgenome} -input ${dir2}/${recalSNP} $(if [ -n "${popgvcf}" ]; then echo "-aggregate ${popgvcf}"; fi) --maxGaussians 4 -resource:mills,known=false,training=true,truth=true,prior=12.0 ${millsgold} -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbSNP} -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 90.0 -recalFile ${dir2}/recalibrate_INDEL.recal -tranchesFile ${dir2}/recalibrate_INDEL.tranches -rscriptFile ${dir2}/recalibrate_INDEL_plots.R -nt ${threads} $(if [ -n "${ped}" ]; then echo "-ped ${ped} -pedValidationType SILENT"; fi) || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
		# Pass #4 for Indels ApplyRecalibration
		echo "java -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${gatk} -T ApplyRecalibration -R ${fasta_refgenome} -input ${dir2}/${recalSNP} -tranchesFile ${dir2}/recalibrate_INDEL.tranches -recalFile ${dir2}/recalibrate_INDEL.recal -o ${dir2}/$filteredSNP --ts_filter_level 99.0 -mode INDEL -nt ${threads} $(if [ -n "${ped}" ]; then echo "-ped ${ped} -pedValidationType SILENT"; fi) || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
	else
		echo "GATK4 command; not implemented yet"
		exit
		"${gatk}" --java-options "-Xmx${mem}g -Djava.io.tmpdir=${tmp}"
	fi

	# Cleaning commands
	if [ "${debug}" != "1" ]
	then
	# Remove indermediate files
		echo "rm ${dir2}/${recalSNP}"
		echo "rm ${dir2}/${recalSNP}.idx"
		echo "rm ${dir2}/recalibrate_SNP.tranches"
		echo "rm ${dir2}/recalibrate_SNP.recal"
		echo "rm ${dir2}/recalibrate_SNP.recal.idx"
		echo "rm ${dir2}/recalibrate_INDEL.tranches"
		echo "rm ${dir2}/recalibrate_INDEL.recal"
		echo "rm ${dir2}/recalibrate_INDEL.recal.idx"
	# Remove .sbatch
		echo "rm ${dir2}/${samplename}_${job}.sbatch"
	fi

} > "${dir2}/${samplename}_${job}.sbatch"

# Queue job
SBjointgeno=$(until sbatch "${dir2}/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
echo ""
echo "-- Joint genotyping job queued --"



######################
## Annotate SNPs using ANNOVAR
job="annovar"
# We are no longer processing GVCFs
samplename="SNPs"

# ANNOVAR variables
annovarfile="${LB}-JointGenotypeOutput.annovar"
snpssummary="${LB}-JointGenotypeOutput.snps"

{
	# General SLURM parameters
	echo '#!/bin/bash'
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append"
	echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
	echo "#SBATCH --time=4:00:00"
	if [ -n "${SLURMemail}" ]
	then
		echo "#SBATCH --mail-type=FAIL --mail-user=${SLURMemail}"
	fi
	if [ -n "${SLURMaccount}" ]
	then
		echo "#SBATCH --account=${SLURMaccount}"
	fi
	if [ -n "${SLURMpartition}" ]
	then
		echo "#SBATCH --partition=${SLURMpartition}"
	fi
	if [ -n "${SLURMqos}" ]
	then
		echo "#SBATCH --qos=${SLURMqos}"
	fi

	# Require previous job successful completion
	echo "#SBATCH --dependency=afterok:${SBjointgeno##* }"

	# General commands
	echo "mkdir -p ${tmp}"
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}"
	fi

	# Job specific commands
	# Convert to ANNOVAR format from GATK .vcf file
	echo "${annovar}/convert2annovar.pl --format vcf4 -allsample -withfreq --includeinfo ${dir2}/${filteredSNP} --outfile ${dir2}/${annovarfile} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

	# Annotate using ANNOVAR
	#echo "${annovar}/table_annovar.pl --buildver hg19 ${dir2}/${annovarfile} ${annovar}/humandb/ --protocol refGene,phastConsElements46way,genomicSuperDups,gwasCatalog,esp6500siv2_all,1000g2015aug_all,exac03,gme,kaviar_20150923,avsnp147,dbnsfp33a,dbnsfp31a_interpro,dbscsnv11,clinvar_20170130,intervar_20170202,hrcr1,revel,mcap --operation g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f --otherinfo --outfile ${dir2}/${snpssummary} --remove -thread ${threads}" >> ${dir2}/${samplename}_"${job}".sbatch
	echo "${annovar}/table_annovar.pl --buildver ${buildver} ${dir2}/${annovarfile} ${annovar}/humandb/ --protocol ${protocol} --operation ${operation} --otherinfo --outfile ${dir2}/${snpssummary} --remove -thread ${threads} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

	# Fixing headers and cleaning files
	echo "sed -i \"1s/Otherinfo.*$/Otherinfo1\tOtherinfo2\tOtherinfo3\t\$(grep CHROM ${dir2}/${filteredSNP} | sed 's/#//g')/g\" ${dir2}/${snpssummary}.hg19_multianno.txt"		# Fixing headers to add back sample names in annovar csv output files
	echo "sed -i 's/\\x2c/,/g' ${dir2}/${snpssummary}.hg19_multianno.txt"		# ClinVar databases are sometimes full of \x2c (comma in hexadecimal), clean it up!
	# May be overkill, converting .txt to .csv
	echo "cp ${dir2}/${snpssummary}.hg19_multianno.txt ${dir2}/${snpssummary}.hg19_multianno.csv"
	echo "sed -i 's/,/;/g' ${dir2}/${snpssummary}.hg19_multianno.csv"			# Remove any potential existing commas and replace them by semi columns
	echo "sed -i 's/\t/,/g' ${dir2}/${snpssummary}.hg19_multianno.csv"			# Convert tabs to commas

	# Cleaning commands
	if [ "${debug}" != "1" ]
	then
	# Remove .sbatch
		echo "rm ${dir2}/${samplename}_${job}.sbatch"
	fi

} > "${dir2}/${samplename}_${job}.sbatch"

# Queue job
SBannovar=$(until sbatch "${dir2}/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
echo ""
echo "-- Annovar job queued --"



######################
## Final processing of result files
job="notif"
samplename="$(basename "${dir2}")-WES"

{
# General SLURM parameters
	echo '#!/bin/bash'
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append"
	echo "#SBATCH $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi)"
	echo "#SBATCH --time=10:00"
	if [ -n "${SLURMemail}" ]
	then
		echo "#SBATCH --mail-type=FAIL,END --mail-user=${SLURMemail}"
	fi
	if [ -n "${SLURMaccount}" ]
	then
		echo "#SBATCH --account=${SLURMaccount}"
	fi
	if [ -n "${SLURMpartition}" ]
	then
		echo "#SBATCH --partition=${SLURMpartition}"
	fi
	if [ -n "${SLURMqos}" ]
	then
		echo "#SBATCH --qos=${SLURMqos}"
	fi

	# Require previous job successful completion
	echo "#SBATCH --dependency=afterok:${SBannovar##* }$(if [ -n "${SBcoverageIDs}" ]; then echo "${SBcoverageIDs}"; fi)"

	# General commands
	echo "mkdir -p ${tmp}"
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}"
	fi

	# Job specific commands
	# Creating folders for archive
	echo "[ -f ${dir2}/Results.tar.gz ] && rm ${dir2}/Results.tar.gz"
	echo "[ -d ${dir2}/Results_$(basename "${dir2}") ] && rm -rf ${dir2}/Results"
	echo "mkdir -p ${dir2}/Results_$(basename "${dir2}")/"
	# List all csv, html, pdf and zip files
	echo "csvarchive=\`ls ${dir2}/*.csv\`"
	echo "pdfarchive=\`ls ${dir2}/*.pdf\`"
	echo "htmlarchive=\`ls ${dir2}/*.html\`"
	echo "ziparchive=\`ls ${dir2}/*.zip\`"
	# Create an archive of all csv, html, pdf and zip files for easy download
	echo "for i in \${csvarchive}; do cp -rf \$i ${dir2}/Results_$(basename "${dir2}")/\`basename \$i\`; done"
	echo "for i in \${htmlarchive}; do cp -rf \$i ${dir2}/Results_$(basename "${dir2}")/\`basename \$i\`; done"
	echo "for i in \${pdfarchive}; do cp -rf \$i ${dir2}/Results_$(basename "${dir2}")/\`basename \$i\`; done"
	echo "for i in \${ziparchive}; do cp -rf \$i ${dir2}/Results_$(basename "${dir2}")/\`basename \$i\`; done"
	echo "tar --remove-files -C ${dir2} -pczf ${dir2}/Results.tar.gz Results_$(basename "${dir2}")"
	if [ -n "${iftttkey}" ] && [ -n "${iftttevent}" ]
	then
		# Trigger IFTTT maker channel event when it's ready, nice isn't it?
		echo "curl -X POST -H \"Content-Type: application/json\" -d '{ \"value1\" : \"$(basename "${dir2}")\" , \"value2\" : \"$(whoami)\"}' https://maker.ifttt.com/trigger/${iftttevent}/with/key/${iftttkey}"
	fi

	# Cleaning commands
	if [ "${debug}" != "1" ]
	then
		# Remove error files upon successful completion. Comment to disable.
		echo "rm ${dir2}/*.err"
		# Remove logs folder upon successfull completion. Comment to disable.
		echo "rm -rf ${dir2}/logs"
		# Remove Temporary directory
		echo "rm -rf ${tmp}"
		# Remove .sbatch
		echo "rm ${dir2}/${samplename}_${job}.sbatch"
	fi

} > "${dir2}/${samplename}_${job}.sbatch"

# Queue job
until sbatch "${dir2}/${samplename}_${job}.sbatch" >/dev/null 2>&1; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done
echo ""
echo "-- Results ready notification job queued --"
echo ""

# Clean temporary files
if [ "${debug}" != "1" ]
then
	[ -f "${dir2}/files1" ] && rm "${dir2}/files1"
	[ -f "${dir2}/files2" ] && rm "${dir2}/files2"
	[ -f "${dir2}/Fastqs" ] && rm "${dir2}/Fastqs"
fi

# That's all folks!
exit 0
