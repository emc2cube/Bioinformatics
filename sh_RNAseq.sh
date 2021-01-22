#!/bin/bash
# shellcheck disable=SC2028,SC2030,SC2031
#
# Usage: sh_RNAseq.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini] [OutputDirName]
#
##############################################################
##                       Description                        ##
##############################################################
#
# This script will process fastq(.gz) files and align them to
# a reference genome using either STAR (recommended), hishat2
# or tophat2.
# If STAR is used then RSEM will also be used and differential
# expression will be analysed using DESeq2.
# Differential expression can also be computed using cufflinks
# (cufflinks is pretty much deprecated, should be avoided
# unless trying to reproduce old results).
# Local Splicing Variation can now be computed using MAJIQ
# and/or LeafCutter.
# If a 4th '[OutputDirName]' argument is provided only the
# secondary analyses selected in the config file will be queued
# using the already aligned and processed files from a previous
# run, and results will be saved in a '_OutputDirName' directory.
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
# Select which software to use for alignment jobs.
# This should be installed in your $PATH and corresponding indexes need to be created.
# Valid values: star, hisat2, tophat2
align="star"
#
# tophat2 (bowtie2) indexed reference genome location.
th_refgenome="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
#
# hisat2 indexed reference genome location.
ha_refgenome="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Sequence/Hisat2Index/genome"
#
# STAR indexed reference genome location.
st_refgenome="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Sequence/StarIndex/"
#
# RSEM indexed reference genome location.
rs_refgenome="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Sequence/RsemIndex/genome"
#
# Annotation files folder.
# It must be compatible with the indexed reference genome used for alignment.
# This folder may contain your gtf / gff3 files.
# If tophat2 is used it should contains the known.gff, known.fa, etc…, files
# else they will be generated during the first run.
annotations="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes"
#
# GTF formated annotation file.
# By default looks for genes.gtf located in the ${annotation} folder.
gtf="${annotations}/genes.gtf"
#
# Reference genome (fasta) file.
# It must be compatible with your indexed genomes, and annotation files.
fasta_refgenome="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
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
## Differential Expression options
#
# Select which software to use for differential expression jobs. Leave empty if not needed.
# Corresponding software or modules should be installed in your $PATH
# Valid values: deseq2, cufflinks, both
diffexp="deseq2"
#
# Delimit replicate samples by a COMMA and groups by a SPACE. Sample name is case sensitive.
# example:
# groupedsamples="1d-Control-Rep1,1d-Control-Rep2,1d-Control-Rep3 5d-Control-Rep1,5d-Control-Rep2,5d-Control-Rep3 1d-Experiment-Rep1,1d-Experiment-Rep2,1d-Experiment-Rep3 5d-Experiment-Rep1,5d-Experiment-Rep2,5d-Experiment-Rep3"
groupedsamples=""
#
# Delimit group labels by a COMMA. Groups should be named in the same order than the samples.
# example:
# labels="1d-Control,5d-Control,1d-Experiment,5d-Experiment"
labels=""
#
#
## Local Splicing Variation options
#
# Select which software to use for LSV jobs. Leave empty if not needed.
# Corresponding software or modules should be installed in your $PATH
# or location configured properly below.
# Valid values: majiq, leafcutter, both
splicing="majiq"
#
# MAJIQ options
#
# Length of the RNA-seq reads. Not needed since MAJIQ 2.2, uncomment in the code if using an older MAJIQ version.
# MAJIQ can handle experiments with multiple read lengths, just indicate the longest read length.
# Leave empty for it to be automatically detected, can be overwritten here if needed.
# example:
# readlength="250"
#readlength=""
#
# Preparation strandness.
# Some of the RNASeq preparations are strand specific.
# This preparations can be reverse strand specific [reverse],
# forward strand specific [forward], or non-strand specific [none].
# This parameter is optional, in which case [none] is used as default.
# Valid values: forward, reverse, none
strandness="none"
#
# Transcriptome file with the annotation database to be used with MAJIQ.
# It must be compatible with the indexed reference genome used for alignment.
# Currently, only accepts GFF3 format. By default looks for genes.gff3 located in the ${annotation} folder.
gff3="${annotations}/genes.gff3"
#
# Genome assembly.
# Name of the genome assembly used. Can also be determined automatically from your directory structure (see examples).
# examples:
# genomebuild="hg19"
# genomebuild=$(echo ${gff3} | sed 's/\/Annotation.*//g' | awk -F/ '{print $NF}')
genomebuild=$(echo ${gff3} | sed 's/\/Annotation.*//g' | awk -F/ '{print $NF}')
#
# Perform PSI quantification.
# 0 = No ; 1 = Yes
psi="0"
#
# LeafCutter options
#
# LeafCutter folder, should contain the scripts, clustering and leafviz folders.
leafCutterDir="/Tools/leafcutter"
#
# Optional exon_file location.
# For more details see http://davidaknowles.github.io/leafcutter/articles/Usage.html
exonfile="${annotations}/leafcutter_exons.txt.gz"
#
# Annotation code for Leafcutter Shiny App.
# path to annotation files and <annotation_code> used in step 0 of
# http://davidaknowles.github.io/leafcutter/articles/Visualization.html
# From Step 1 concrete example section this would be:
# anncode="../leafviz/annotation_codes/gencode_hg19/gencode_hg19"
anncode="${annotations}/leafcutter"
#
#
## Notifications options
#
# Event name used in your webhooks or IFTTT recipes.
# The IFTTT maker channel will look for the combination private key + event name to then trigger your recipe.
# You can then create a recipe to send an email, a text message or a push notification.
# In the case of a Slack webhook this will be used in the message title
notifevent="RNAseq"
#
# Trigger a Slack Webhook when script is done.
# You must create an "Incoming WebHooks" associated to your slack workspace on https://slack.com/apps
# Copy your full private webhook URL here. Leave blank to disable this function.
# slackwebhookURL="https://hooks.slack.com/services/T5HF5GTUK85/AHUIK456HJG/GSD27f5gGQ7SD5r2fg" # Not a real webhook URL, you have to use your own private one.
slackwebhookURL=""
#
# Trigger IFTTT when script is done.
# You must register the "Maker channel" on https://ifttt.com/maker and create an event (limit of 3 free event since September 2020...)
# Copy your private key here. Leave blank to disable this function.
# iftttkey="AbCd_15CdhUIvbsFJTHGMcfgjsdHRTgcyjt" # Not a real key, you have to use your own private one.
iftttkey=""
#
#
## Setup done. You should not need to edit below this point ##

# Help!
if [ "${1}" = "--help" ] || [ "${2}" = "--help" ] || [ "${3}" = "--help" ] || [ "${4}" = "--help" ]
then
	echo "Usage: $(basename "${0}") </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini] [OutputDirName]"
	echo ""
	echo "Description"
	echo ""
	echo "This script will process fastq(.gz) files and align them to"
	echo "a reference genome using either STAR (recommended), hishat2"
	echo "or tophat2."
	echo "If STAR is used then RSEM will also be used and differential"
	echo "expression will be analysed using DESeq2."
	echo "Differential expression can also be computed using Cufflinks"
	echo "(Cufflinks is pretty much deprecated, should be avoided"
	echo "unless trying to reproduce old results)"
	echo "If a 4th '[OutputDirName]' argument is provided"
	echo "only the secondary analyses selected in the config file will be queued"
	echo "using the already aligned and processed files from a previous run, "
	echo "and results will be saved in a '_OutputDirName' directory."
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
	echo "$(basename "${0}") version 2.1.2"
	echo "Now easier to add analyses in a rerun and improved default R graphs (2.1.1)"
	echo "Alternative splicing support (2.1)"
	echo "DESeq2 and Single-Read support (2.0)"
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

if [ -z "${labels}" ] || [ -z "${groupedsamples}" ]
then
	echo "Error: both \"groupedsamples\" and \"labels\" variables need to be configured properly"
	exit
fi

if [ -n "${4}" ]
then
	reanalysis="_${4}"
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
echo "-- Alignment options --"
echo ""
if [ "${align}" = "star" ]
then
	echo "Using ${st_refgenome} as reference genome"
elif [ "${align}" = "hisat2" ]
then
	echo "Using ${ha_refgenome} as reference genome"
elif [ "${align}" = "tophat2" ]
then
	echo "Using ${th_refgenome} as reference genome"
fi
if [ "${underdet}" = "0" ]
then
	echo "Undetermined files will not be processed"
fi
if [ "${blank}" = "1" ]
then
	echo "\"BLANK\" files will be processed"
fi
if [ "${merge}" = "1" ]
then
	echo "Individual files will also be merged into a big \"MERGED\" file"
fi
if [ "${trim}" = "1" ]
then
	echo "Sequence will be trimmed"
	echo "Trimmomatic is installed in ${Trimmomatic}"
fi

if [ -n "${diffexp}" ]
then
	echo ""
	echo "-- Differential Expression options --"
	echo ""
	echo "Differential expression will be analysed using $(if [ "${diffexp}" = "deseq2" ] || [ "${diffexp}" = "both" ]; then echo "DESeq2"; fi)$(if [ "${diffexp}" = "both" ]; then echo " and "; fi)$(if [ "${diffexp}" = "cufflinks" ] || [ "${diffexp}" = "both" ]; then echo "Cufflinks"; fi)."
	if [ "${diffexp}" = "cufflinks" ] || [ "${diffexp}" = "both" ]
	then
		echo "WARNING - The use of Cufflinks for RNAseq analysis is pretty much deprecated, should be avoided - WARNING"
	fi
fi

if [ -n "${splicing}" ]
then
	echo ""
	echo "-- Local Splicing Variation options --"
	echo ""
	echo "Local Splicing Variants will be analysed using $(if [ "${splicing}" = "majiq" ] || [ "${splicing}" = "both" ]; then echo "MAJIQ"; fi)$(if [ "${splicing}" = "both" ]; then echo " and "; fi)$(if [ "${splicing}" = "leafcutter" ] || [ "${splicing}" = "both" ]; then echo "LeafCutter"; fi)."
fi

echo ""
echo "-- Hardware --"
echo ""
echo "Up to ${threads} CPUs will be used"
echo "Up to ${mem}GB of memory will be allocated to the programs"

echo ""
echo "---------------"

# Initialize
mkdir -p "${dir2}/"
[ -f "${dir2}/files1" ] && rm "${dir2}/files1"
[ -f "${dir2}/files2" ] && rm "${dir2}/files2"
[ -f "${dir2}/Fastqs" ] && rm "${dir2}/Fastqs"
[ -f "${dir2}/matrices" ] && rm "${dir2}/matrices"
[ -f "${dir2}/assembly_GTF_list.txt" ] && rm "${dir2}/assembly_GTF_list.txt"
[ -d "${dir2}/logs" ] && rm -rf "${dir2}/logs"
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
		find -L "${dir}" -maxdepth 1 -name '*_R1*' -not -name '*ndetermined*' -not -name '*nmatched*' -not -name 'BLANK*' | sed 's#.*/##' | sort -n > "${dir2}/files1"
		find -L "${dir}" -maxdepth 1 -name '*_R2*' -not -name '*ndetermined*' -not -name '*nmatched*' -not -name 'BLANK*' | sed 's#.*/##' | sort -n > "${dir2}/files2"
		paste "${dir2}/files1" "${dir2}/files2" > "${dir2}/Fastqs"
	else
		# Remove all Undetermined_* files
		find -L "${dir}" -maxdepth 1 -name '*_R1*' -not -name '*ndetermined*' -not -name '*nmatched*' | sed 's#.*/##' | sort -n > "${dir2}/files1"
		find -L "${dir}" -maxdepth 1 -name '*_R2*' -not -name '*ndetermined*' -not -name '*nmatched*' | sed 's#.*/##' | sort -n > "${dir2}/files2"
		paste "${dir2}/files1" "${dir2}/files2" > "${dir2}/Fastqs"
	fi
else
	if [ "${blank}" = "0" ]
	then
		# Remove all BLANK* files
		find -L "${dir}" -maxdepth 1 -name '*_R1*' -not -name 'BLANK*' | sort -n | sed 's#.*/##' > "${dir2}/files1"
		find -L "${dir}" -maxdepth 1 -name '*_R2*' -not -name 'BLANK*' | sort -n | sed 's#.*/##' > "${dir2}/files2"
		paste "${dir2}/files1" "${dir2}/files2" > "${dir2}/Fastqs"
	else
		# Process all the files!
		find -L "${dir}" -maxdepth 1 -name '*_R1*' | sed 's#.*/##' | sort -n > "${dir2}/files1"
		find -L "${dir}" -maxdepth 1 -name '*_R2*' | sed 's#.*/##' | sort -n > "${dir2}/files2"
		paste "${dir2}/files1" "${dir2}/files2" > "${dir2}/Fastqs"
	fi
fi



######################
## Processing samples
######################

if [ -z "${reanalysis}" ]
then

	echo ""
	echo "-- Queuing sample jobs --"

	while read -r line;
	do

		# General variables
		read1=$(echo "${line}" | cut -f1)
		read2=$(echo "${line}" | cut -f2)
		samplename=$(echo "${read1}" | awk -F_R1 '{print $1}')
		if [ "${read1}" != "${read2}" ]
		then
			pairedend="1"
		fi

		echo ""
		echo -e "\t-- Processing" "${samplename}"
		echo ""



		######################
		## Prepare and queue trimming job

		job="trim"

		# Variables
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

		if [ -n "${trim}" ] && [ "${trim}" = "1" ]
		then
			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
				echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi) $(if [ -n "${threads}" ] && [ "${threads}" -gt "2" ]; then echo "--cpus-per-task=2"; else echo "--cpus-per-task=${threads}"; fi)"
				echo "#SBATCH --requeue"
				echo "#SBATCH --time=2:00:00"
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
				echo "mkdir -p \"${tmp}\""
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				# Trim fastq files with trimmomatic
				if [ -n "${pairedend}" ] && [ "${pairedend}" = "1" ]
				then
					echo "java -Xmx$(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "8"; else echo "${mem}"; fi)g -Djava.io.tmpdir=\"${tmp}\" -jar ${Trimmomatic}/trimmomatic.jar PE -threads $(if [ -n "${threads}" ] && [ "${threads}" -gt "2" ]; then echo "2"; else echo "${threads}"; fi) -phred33 ${dir}/${read1} ${dir}/${read2} ${trimout1} ${unpaired1} ${trimout2} ${unpaired2} ILLUMINACLIP:${Trimmomatic}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
				else
					echo "java -Xmx$(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "8"; else echo "${mem}"; fi)g -Djava.io.tmpdir=\"${tmp}\" -jar ${Trimmomatic}/trimmomatic.jar SE -threads $(if [ -n "${threads}" ] && [ "${threads}" -gt "2" ]; then echo "2"; else echo "${threads}"; fi) -phred33 ${dir}/${read1} ${trimout1} ILLUMINACLIP:${Trimmomatic}/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
				fi

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/${job}_${samplename}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/${job}_${samplename}.sbatch"

			# Queue job
			SBtrim=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			echo -e "\t Trimming job queued"
		else
			# Need to rename some variables
			trimout1="${dir}/${read1}"
			trimout2="${dir}/${read2}"
		fi



		######################
		## Align reads with hisat2 / tophat2 / star

		job="align"

		# Variables
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
			CN=$(gzip -cd "${dir}/${read1}" | head -n 1 | awk -F: '{print $1}' | sed 's/@//')
			PU=$(gzip -cd "${dir}/${read1}" | head -n 1 | awk -F: '{print $3}')
		else
			CN=$(head -n 1 "${dir}/${read1}" | awk -F: '{print $1}' | sed 's/@//')
			PU=$(head -n 1 "${dir}/${read1}" | awk -F: '{print $3}')
		fi

		{
			# General SLURM parameters
			echo '#!/bin/bash'
			echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
			echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
			echo "#SBATCH --requeue"
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
			if [ -n "${trim}" ] && [ "${trim}" = "1" ]
			then
				echo "#SBATCH --dependency=afterok:${SBtrim##* }"
			fi

			# General commands
			echo "mkdir -p \"${tmp}\""
			if [ -n "${customcmd}" ]
			then
				echo "${customcmd}"
			fi

			# Job specific commands
			if [ ${align} = "star" ]
			then
				# star alignment job using ENCODE standard options
				echo "STAR --runMode alignReads --runThreadN ${threads} --limitBAMsortRAM ${mem}000000000 --outTmpDir \"${tmp}\"/${samplename}_aRtmp --twopassMode Basic --genomeDir ${st_refgenome} --readFilesIn ${trimout1} $(if [ -n "${pairedend}" ] && [ "${pairedend}" = "1" ]; then echo "${trimout2}"; fi) $(if [ -n "${fastqgz}" ] && [ ${trim} -ne "1" ]; then echo "--readFilesCommand zcat"; fi) --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --quantMode TranscriptomeSAM GeneCounts --alignEndsType EndToEnd --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${dir2}/${samplename}. || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
				echo "mv ${dir2}/${samplename}.Aligned.sortedByCoord.out.bam ${dir2}/${bamsortedout}"
			elif [ ${align} = "hisat2" ]
			then
				# hisat2 alignment job
				echo "hisat2 -p ${threads} --phred33 --dta-cufflinks --no-softclip --rg-id ${LB}_${SM} --rg CN:${CN} --rg LB:${LB} --rg PL:${PL} --rg PU:${PU} --rg SM:${SM} -x ${ha_refgenome} -1 ${trimout1} $(if [ -n "${pairedend}" ] && [ "${pairedend}" = "1" ]; then echo "-2 ${trimout2}"; fi) -S ${dir2}/${samout} || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
			elif [ ${align} = "tophat2" ]
			then
				# tophat2 alignment job
				echo "tophat2 -G ${gtf} --transcriptome-index ${annotations}/known -p ${threads} -o ${dir2} ${th_refgenome} ${trimout1} $(if [ -n "${pairedend}" ] && [ "${pairedend}" = "1" ]; then echo "${trimout2}"; fi) || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
			fi
			if [ ${align} = "hisat2" ] || [ ${align} = "tophat2" ]
			then
				# Convert .sam to .bam
				echo "samtools view -bS -o ${dir2}/${bamout} ${dir2}/${samout} || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
				# Sort .bam file
				echo "samtools sort -@ ${threads} -o ${dir2}/${bamsortedout} -O bam -T \"${tmp}\" ${dir2}/${bamout} || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
			fi
			# Index sorted .bam
			echo "samtools index ${dir2}/${bamsortedout} || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

			# Cleaning commands
			# Remove trim.fastq files from destination folder
			if [ ${trim} = "1" ] && [ ${fastqc} = "0" ]
			then
				echo "rm ${trimout1} ${trimout2}"
			fi
			# Remove .sam file
			echo "rm ${dir2}/${samout}"
			# Remove unsorted .bam file
			echo "rm ${dir2}/${bamout}"
			# Clean STAR intermediate files
			echo "rm -rf ${dir2}/${samplename}.Log.progress.out ${dir2}/${samplename}._STARgenome ${dir2}/${samplename}._STARpass1"
			echo "mv ${dir2}/${samplename}.Log.out ${dir2}/${logs}/${samplename}.star.log"
			echo "mv ${dir2}/${samplename}.Log.final.out ${dir2}/${logs}/${samplename}.star.stats"
			# Remove .sbatch
			echo "rm ${dir2}/${job}_${samplename}.sbatch"
			echo "exit 0"

		} > "${dir2}/${job}_${samplename}.sbatch"

		# Queue job
		SBalign=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
		SBalignIDs=${SBalignIDs}:${SBalign##* }
		echo -e "\t Alignment job queued"



		######################
		## Run fastqc to generate quality control files

		job="fqc"

		if [ "${fastqc}" = "1" ]
		then
			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
				echo "#SBATCH $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi)"
				echo "#SBATCH --requeue"
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
				if [ -n "${SBtrim}" ]
				then
					echo "#SBATCH --dependency=afterok:${SBalign##* }"
				fi

				# General commands
				echo "mkdir -p \"${tmp}\""
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				echo "fastqc -o ${dir2}/ --noextract ${trimout1} ${trimout2} || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove trim.fastq files from destination folder
					if [ ${trim} = "1" ]
					then
						echo "rm ${trimout1} ${trimout2}"
					fi
					# Remove .sbatch
					echo "rm ${dir2}/${job}_${samplename}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/${job}_${samplename}.sbatch"

			# Queue job
			SBfqc=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			if [ -n "${SBfqc##* }" ]
			then
				SBfqcIDs=${SBfqcIDs}:${SBfqc##* }
			fi
			echo -e "\t FastQC job queued"
		fi



		######################
		## Count reads

		job="rsem-calcexp"

		if [ "${align}" = "star" ] && { [ "${diffexp}" = "deseq2" ] || [ "${diffexp}" = "both" ]; }
		then
			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
				echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
				echo "#SBATCH --requeue"
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
				echo "#SBATCH --dependency=afterok:${SBalign##* }"

				# General commands
				echo "mkdir -p \"${tmp}\""
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				echo "rsem-calculate-expression -p ${threads} --temporary-folder \"${tmp}\"/${samplename} $(if [ -n "${pairedend}" ] && [ "${pairedend}" = "1" ]; then echo "--paired-end"; fi) --bam --no-bam-output --calc-ci ${dir2}/${samplename}.Aligned.toTranscriptome.out.bam ${rs_refgenome} ${dir2}/${samplename}.rsem || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
				echo "mv ${dir2}/${samplename}.rsem.stat ${dir2}/${samplename}.rsem"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					echo "rm ${dir2}/${samplename}.Aligned.toTranscriptome.out.bam"
					# Remove .sbatch
					echo "rm ${dir2}/${job}_${samplename}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/${job}_${samplename}.sbatch"

			# Queue job
			SBrsem=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			SBrsemIDs=${SBrsemIDs}:${SBrsem##* }
			echo -e "\t RSEM counting job queued"
		fi



		######################
		## Assemble transcriptomes

		job="cufflinks"

		if [ "${diffexp}" = "cufflinks" ] || [ "${diffexp}" = "both" ]
		then
			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
				echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
				echo "#SBATCH --requeue"
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
				echo "#SBATCH --dependency=afterok:${SBalign##* }"

				# General commands
				echo "mkdir -p \"${tmp}\""
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				# Create a sample specific folder, and as cufflinks do not allow to specify an output directory then cd to it.
				echo "mkdir -p ${dir2}/${samplename}.cuff"
				echo "cd ${dir2}/${samplename}.cuff || exit"
				# Create cufflinks job
				echo "cufflinks -u -p ${threads} -g ${gtf} -b ${fasta_refgenome} ${dir2}/${bamsortedout} || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/${job}_${samplename}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/${job}_${samplename}.sbatch"

			# Queue job
			SBcufflinks=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			SBcufflinksIDs=${SBcufflinksIDs}:${SBcufflinks##* }
			echo -e "\t Cufflinks job queued"
		fi



		######################
		## LeafCutter convert bams to juncs

		job="leaf-junc"

		if [ "${splicing}" = "leafcutter" ] || [ "${splicing}" = "both" ]
		then
			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
				echo "#SBATCH $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi)"
				echo "#SBATCH --requeue"
				echo "#SBATCH --time=2:00:00"
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
				echo "mkdir -p \"${tmp}\""
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				# adapted from https://github.com/davidaknowles/leafcutter/blob/master/scripts/bam2junc.sh
				echo "samtools view ${dir2}/${samplename}.sorted.bam | ${leafCutterDir}/scripts/filter_cs.py | ${leafCutterDir}/scripts/sam2bed.pl --use-RNA-strand - ${dir2}/${samplename}.bed || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
				echo "${leafCutterDir}/scripts/bed2junc.pl ${dir2}/${samplename}.bed ${dir2}/${samplename}.junc || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					echo "rm ${dir2}/${samplename}.bed"
					# Remove .sbatch
					echo "rm ${dir2}/${job}_${samplename}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/${job}_${samplename}.sbatch"

			# Queue job
			SBleafjunc=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			SBleafjuncIDs=${SBleafjuncIDs}:${SBleafjunc##* }
			echo -e "\t JUNC conversion job queued"
		fi

	done < "${dir2}/Fastqs"

fi # End of initial run of sample specific jobs



######################
## RSEM

if [ "${align}" = "star" ] && { [ "${diffexp}" = "deseq2" ] || [ "${diffexp}" = "both" ]; }
then

	echo ""
	echo "-- RSEM --"



	######################
	## RSEM Count reads

	# Will only run in reanalysis mode if the .rsem folder was not previously created during the initial run.
	# Variables
	job="rsem-calcexp"
	firstsample="1"

	while read -r line;
	do

		# General variables
		read1=$(echo "${line}" | cut -f1)
		samplename=$(echo "${read1}" | awk -F_R1 '{print $1}')

		if [ -n "${reanalysis}" ] && [ ! -d "${dir2}/${samplename}.rsem" ]
		then
			if [ "${firstsample}" = "1" ]
			then
				echo ""
				echo -e "\t-- Queuing missing count jobs"
				firstsample="0"
			fi

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
				echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
				echo "#SBATCH --requeue"
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
				# No previous jobs expected

				# General commands
				echo "mkdir -p \"${tmp}\""
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				echo "rsem-calculate-expression --seed-length 10 -p ${threads} --temporary-folder \"${tmp}\"/${samplename} $(if [ -n "${pairedend}" ] && [ "${pairedend}" = "1" ]; then echo "--paired-end"; fi) --bam --no-bam-output --calc-ci ${dir2}/${samplename}.Aligned.toTranscriptome.out.bam ${rs_refgenome} ${dir2}/${samplename}.rsem || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
				echo "mv ${dir2}/${samplename}.rsem.stat ${dir2}/${samplename}.rsem"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					echo "rm ${dir2}/${samplename}.Aligned.toTranscriptome.out.bam"
					# Remove .sbatch
					echo "rm ${dir2}/${job}_${samplename}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/${job}_${samplename}.sbatch"

			# Queue job
			SBrsem=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			SBrsemIDs=${SBrsemIDs}:${SBrsem##* }
			echo -e "\t ${samplename} missing RSEM counting job queued"
		fi
	done < "${dir2}/Fastqs"



	######################
	## RSEM Generate data matrix

	# Variables
	job="rsem-generate"
	# We are no longer processing sample files, but all the rsem files to generate the matrix
	samplename="matrix"

	# Variables
	IFS="," read -r -a labelsarray <<< "${labels}"
	IFS=" " read -r -a samplesarray <<< "${groupedsamples}"
	declare -a replicatesarray

	{
		# General SLURM parameters
		echo '#!/bin/bash'
		echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
		echo "#SBATCH $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi)"
		echo "#SBATCH --requeue"
		echo "#SBATCH --time=15:00"
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
		if [ -n "${SBrsemIDs}" ]
		then
			echo "#SBATCH --dependency=afterok${SBrsemIDs}"
		fi

		# General commands
		echo "mkdir -p \"${tmp}\""
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}"
		fi

		# Job specific commands
		mkdir -p "${dir2}"/DESeq2"${reanalysis}"

		# Generate an array with sample replicate names
		repN="0"
		while [ ${repN} -lt ${#labelsarray[@]} ]
		do
			replicatesarray[${repN}]="${samplesarray[${repN}]//,/ }"
			repN=$((repN + 1))
		done

		# Generate pairwise comparison matrices
		condA="0"
		condB="0"
		while [ ${condA} -lt ${#labelsarray[@]} ]
		do
			condB=$((condA + 1))
			while [ ${condB} -lt ${#labelsarray[@]} ]
			do
				echo "rsem-generate-data-matrix ${dir2}/$(echo "${replicatesarray[${condA}]}" "${replicatesarray[${condB}]}" | sed "s# #.rsem.genes.results ${dir2}/#g").rsem.genes.results > ${dir2}/DESeq2${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}.genes.matrix || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
				echo "${dir2}/DESeq2${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}.genes.matrix" >> "${dir2}/matrices"
				{
					echo '#!/usr/bin/env Rscript'
					echo ''
					echo '# Variables'
					echo "table <- \"${dir2}/DESeq2${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}.genes.matrix\""
					echo "cond1 <- \"${labelsarray[${condA}]}\""
					echo "rep1 <- $(wc -w <<< "${replicatesarray[${condA}]}" | sed 's/ //g')			# Number of replicates for cond2, automatically determined"
					echo "cond2 <- \"${labelsarray[${condB}]}\""
					echo "rep2 <- $(wc -w <<< "${replicatesarray[${condB}]}" | sed 's/ //g')			# Number of replicates for cond2, automatically determined"
					echo "output <- \"${dir2}/DESeq2${reanalysis}/\""
					echo "siglfc <- 1			# |Log2| fold change to be used as significant. Set at 1 by default, meaning a 2 fold change (log2(2)=1)"
					echo "sigpadj <- 0.05			# Adjusted p-value used as upper threshold for significance"
				} > "${dir2}/DESeq2${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}.DESeq2.R"
				chmod +x "${dir2}/DESeq2${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}.DESeq2.R"

				condB=$((condB + 1))
			done

			condA=$((condA + 1))
		done

		# Generate comparison matrices with all samples, only if more than 2 conditions
		if [ ${#labelsarray[@]} -gt 2 ]
		then
			condN="0"
			echo "rsem-generate-data-matrix ${dir2}/$(echo "${groupedsamples}" | sed "s# #.rsem.genes.results ${dir2}/#g;s#,#.rsem.genes.results ${dir2}/#g").rsem.genes.results > ${dir2}/DESeq2${reanalysis}/$(basename "${dir2}")_all.genes.matrix || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
			{
				echo '#!/usr/bin/env Rscript'
				echo '# Variables'
				echo "table <- \"${dir2}/DESeq2${reanalysis}/$(basename "${dir2}")_all.genes.matrix\""
				while [ ${condN} -lt ${#labelsarray[@]} ]
				do
					echo "cond$((condN + 1)) <- \"${labelsarray[${condN}]}\""
					conditionsall="$(if [ -n "${conditionsall}" ]; then echo "${conditionsall}, "; fi)rep(cond$((condN + 1)), $(wc -w <<< "${replicatesarray[${condN}]}"))"
					condN=$((condN + 1))
				done
				echo "output <- \"${dir2}/DESeq2${reanalysis}/$(basename "${dir2}")_\""
				echo "siglfc <- 1			# |Log2| fold change to be used as significant. Set at 1 by default, meaning a 2 fold change (log2(2)=1)"
				echo "sigpadj <- 0.05		# Adjusted p-value used as upper threshold for significance"
			} > "${dir2}/DESeq2${reanalysis}/$(basename "${dir2}")_all.DESeq2.R"
			chmod +x "${dir2}/DESeq2${reanalysis}/$(basename "${dir2}")_all.DESeq2.R"
		fi

		# Cleaning commands
		if [ "${debug}" != "1" ]
		then
			# Remove .sbatch
			echo "rm ${dir2}/${job}_${samplename}.sbatch"
		fi
		echo "exit 0"

	} > "${dir2}/${job}_${samplename}.sbatch"

	# Queue job
	SBmatrix=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
	echo ""
	echo -e "\t-- matrix generation job queued"


	######################
	## DESeq2 analysis and basic plots

	job="deseq2"

	echo ""
	echo "-- DESeq2 --"
	echo ""

	# Create R script for all sample matrix, only if more than 2 conditions
	if [ ${#labelsarray[@]} -gt 2 ]
	then
		# Variables
		samplename="$(basename "${dir2}")_all"
		Rscript=${samplename}.DESeq2.R

		{
			# General SLURM parameters
			echo '#!/bin/bash'
			echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
			echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "16" ]; then echo "--mem=16000"; else echo "--mem=${mem}000"; fi)"
			echo "#SBATCH --requeue"
			echo "#SBATCH --time=15:00"
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
			echo "#SBATCH --dependency=afterok:${SBmatrix##* }"

			# General commands
			echo "mkdir -p \"${tmp}\""
			if [ -n "${customcmd}" ]
			then
				echo "${customcmd}"
			fi

			# Job specific commands
			cat >> "${dir2}/DESeq2${reanalysis}/${Rscript}" <<RALLDELIM


# This script uses DESeq2 in R to run differential-expression analysis on a reads count matrix of genes and automatically generate some graphs.


# Load libraries
library(DESeq2)
library(RColorBrewer)
library(calibrate)
library(pheatmap)


# Read in the data and make sample information table
data=read.table(table, header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)											# read the reads count matrix file that was generated by rsem-generate-data-matrix
colnames(data) <- sub('.rsem.genes.results', '', basename(colnames(data)))															# fix sample names in column names
cols=c(1:ncol(data))																												# get the number of columns
data[,cols]=apply(data[,cols], 2, function(x) as.numeric(as.integer(x)))															# format the matrix as integer numbers for DESeq2 to be happy
conditions <- factor(c(${conditionsall}))										# generate a comma-separated list of conditions/treatments for each sample. Replicates should be named the same.
samples=as.data.frame((colnames(data)))																								# make samples list for DESeq2
samples <- cbind(samples, conditions)
colnames(samples)=c("experimental","condition")
groups=factor(samples[,2])


# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromMatrix(countData=as.matrix(data), colData=samples, design =~condition)											# generate a DESeq2 object from the data
colnames(dds) <- colnames(data)
dds <- DESeq(dds)																													# run DESeq2


# Regularized log transformation for clustering/heatmaps, etc
vsd <- vst(dds, blind=FALSE)

# Save differential expression results
res <- results(dds)
table(res\$padj<sigpadj)
## Order by adjusted p-value
res <- res[order(res\$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
## Write results
write.csv(resdata, file=paste0(output, "all-ExpDiff.csv"))


# Variables
conditions_df=data.frame(c(${conditionsall}), row.names=colnames(data))
colnames(conditions_df)=c("Labels")
condcol <- grDevices::rainbow(length(unique(conditions)))												# Distance Matrix conditions colors, use rainbow palette
names(condcol) <- unique(conditions)																	# Proper column names to use in the legend
condcol <- list(Labels = condcol)
dmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)													# Distance Matrix color palette
hmcol <- colorRampPalette(brewer.pal(9, 'RdYlBu'))(100)													# Heatmaps color palette


# Principal Components Analysis
pdf(paste0(output, "all-PCA.pdf"))
plotPCA(vsd, intgroup="condition")
garbage <- dev.off() # Save to file


# Sample Distance Matrix
sampleDists <- dist(t(assay(vsd)))
pdf(paste0(output, "all-DistanceMatrix.pdf"))
pheatmap(as.matrix(sampleDists), col=rev(dmcol), clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames=TRUE, annotation_col=conditions_df, annotation_row=conditions_df, annotation_names_row=FALSE, annotation_names_col=FALSE, annotation_colors=condcol, legend = FALSE, main="all samples Matrix")
garbage <- dev.off() # Save to file


# Heatmap of significant hits ( padj<sigpadj and |log2FoldChange|>=siglfc )
counts <- counts(dds,normalized=TRUE)																								# Get normalized read counts (sorted by padj)
counts <- counts[apply(counts, 1, function(row) all(row !=0 )),]																	# Remove genes with zero reads
## You will probably need to change conditions to perform a meaningful comparison here
respairwise = results(dds, contrast=c("condition",cond1,cond${#labelsarray[@]}))																	# Select comparison to perform (for log2 changes)
sig <- rownames(respairwise[!is.na(respairwise\$padj) & respairwise\$padj<sigpadj & abs(respairwise\$log2FoldChange)>=siglfc,])[1:30]	# Select the first 30 significant hits (sorted by padj)
sig <- sig[!is.na(sig)]
sigcounts <- counts(dds,normalized=TRUE)[sig,]																						# Get normalized read counts (sorted by padj) for significan genes
try({																																# May fail if only 0 or 1 gene is significant
	sigcounts <- sigcounts[apply(sigcounts, 1, function(row) all(row !=0 )),]														# Remove genes with zero reads
})


# By defaults hits are not clustered and thus stay sorted by their padj value
pdf(paste0(output, "all-HeatmapSig.pdf"), onefile=FALSE)
try({																																# May fail if only 0 or 1 gene is significant
	pheatmap(log2(sigcounts), col=rev(hmcol), scale='row', cluster_rows=FALSE, cluster_cols=FALSE, main="all Top Hits")
})
garbage <- dev.off() # Save to file


# But we can cluster them using the following command
pdf(paste0(output, "all-HeatmapSigClust.pdf"), onefile=FALSE)
try({																																# May fail if only 0 or 1 gene is significant
	pheatmap(log2(sigcounts), col=rev(hmcol), scale='row', cluster_rows=TRUE, cluster_cols=TRUE, main="all Clustered Top Hits")
})
garbage <- dev.off() # Save to file


# And even plot every gene (will take time)
pdf(paste0(output, "all-HeatmapAllClust.pdf"), onefile=FALSE)
pheatmap(log2(counts), col=rev(hmcol), scale='row', cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames=FALSE, annotation_col=conditions_df, annotation_names_col=FALSE, annotation_colors=condcol, main="all Clustered Reads Count")
garbage <- dev.off() # Save to file


# Additional optional graphs using DESeqAnalysis https://github.com/acidgenomics/DESeqAnalysis

# Principal Components Analysis with labels
try(library(DESeqAnalysis))
pdf(paste0(output, "all-PCAlabels.pdf"), width=20, height=20)
try(DESeqAnalysis::plotPCA(vsd, label=TRUE, interestingGroups="condition", labels=list(title="all PCA")))
garbage <- dev.off() # Save to file


# Adapted from:
# http://www.bioconductor.org/help/workflows/rnaseqGene/
# https://gist.github.com/stephenturner/f60c1934405c127f09a6
# https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
RALLDELIM
			# Use Rscript
			echo "Rscript ${dir2}/DESeq2${reanalysis}/${Rscript} || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

			# Cleaning commands
			if [ "${debug}" != "1" ]
			then
				# Remove .sbatch
				echo "rm ${dir2}/${job}_${samplename}.sbatch"
			fi
			echo "exit 0"

		} > "${dir2}/${job}_${samplename}.sbatch"

		# Queue job
		SBdeseq2=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
		if [ -n "${SBdeseq2##* }" ]
		then
			SBdeseq2IDs=${SBdeseq2IDs}:${SBdeseq2##* }
		fi
		echo -e "\t-- ${samplename} analysis job queued"
		echo  ""
	fi

	echo -e "\t-- Queuing pairwise comparison analysis jobs"
	echo ""

	# Create R scripts for pairwise matrices
	while read -r line;
	do

		# Variables
		samplename=$(basename "${line}" | sed 's/.genes.matrix//g')
		Rscript=${samplename}.DESeq2.R

		{
			# General SLURM parameters
			echo '#!/bin/bash'
			echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
			echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "16" ]; then echo "--mem=16000"; else echo "--mem=${mem}000"; fi)"
			echo "#SBATCH --requeue"
			echo "#SBATCH --time=10:00"
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
			echo "#SBATCH --dependency=afterok:${SBmatrix##* }"

			# General commands
			echo "mkdir -p \"${tmp}\""
			if [ -n "${customcmd}" ]
			then
				echo "${customcmd}"
			fi

			# Job specific commands
			cat >> "${dir2}/DESeq2${reanalysis}/${Rscript}" <<RDELIM


# This script uses DESeq2 in R to run differential-expression analysis on a reads count matrix of genes and automatically generate some graphs.


# Load libraries
library(DESeq2)
library(RColorBrewer)
library(calibrate)
library(pheatmap)


# Read in the data and make sample information table
data=read.table(table, header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)				# read the reads count matrix file that was generated by rsem-generate-data-matrix
colnames(data) <- sub('.rsem.genes.results', '', basename(colnames(data)))								# fix sample names in column names
cols=c(1:ncol(data))																					# get the number of columns
data[,cols]=apply(data[,cols], 2, function(x) as.numeric(as.integer(x)))								# format the matrix as integer numbers for DESeq2 to be happy
conditions <- factor(c(rep(cond1, rep1), rep(cond2, rep2)))												# generate a comma-separated list of conditions/treatments for each sample. Replicates should be named the same.
samples=as.data.frame((colnames(data)))																	# make samples list for DESeq2
samples <- cbind(samples, conditions)
colnames(samples)=c("experimental","condition")
groups=factor(samples[,2])


# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromMatrix(countData=as.matrix(data), colData=samples, design =~condition)				# generate a DESeq2 object from the data
colnames(dds) <- colnames(data)
dds <- DESeq(dds)																						# run DESeq2


# Regularized log transformation for clustering/heatmaps, etc
vsd <- vst(dds, blind=FALSE)

# Save differential expression results
res <- results(dds)
table(res\$padj<sigpadj)
## Order by adjusted p-value
res <- res[order(res\$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
## Write results
write.csv(resdata, file=paste0(output, cond1, "_vs_", cond2, "-ExpDiff.csv"))


# Variables
conditions_df=as.data.frame(colData(dds)[c("condition")])
colnames(conditions_df)=c("Labels")
condcol <- brewer.pal(8, "Set1")[1:length(unique(conditions))]											# Distance Matrix conditions colors
names(condcol) <- unique(conditions)																	# Proper column names to use in the legend
condcol <- list(Labels = condcol)
dmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)													# Distance Matrix color palette
hmcol <- colorRampPalette(brewer.pal(9, 'RdYlBu'))(100)													# Heatmaps color palette


# MA plot (without and with "significant" genes labeled)
maplot <- function (res, sigthresh=0.05, labelsig=TRUE, textcx=1, ...) {								# sigthresh=0.05 is the default value if ommitted when creating the graph. No need to edit here.
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<sigthresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
resNorm <- lfcShrink(dds, coef=2, type="normal")
try({																									# May fail if apeglm is not installed ( BiocManager::install("apeglm") )
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")										# Will take a while
})

# DESeq2 MA plot function and custom MA plot with significant gene names. Here we use the sigpadj variable value to generate the graph
pdf(paste0(output, cond1, "_vs_", cond2, "-MAplot.pdf"))
plotMA(res, ylim=c(-2,2), main=paste0(cond1, " vs ", cond2, " MA Plot"))
plotMA(resNorm, ylim=c(-2,2), main=paste0(cond1, " vs ", cond2, " normal shrunken log2 MA Plot"))
try({																									# May fail if resLFC was not computed before
plotMA(resLFC, ylim=c(-2,2), main=paste0(cond1, " vs ", cond2, " apeglm shrunken log2 MA Plot"))
})
maplot(resdata, sigpadj, xlab="mean of normalized counts", ylab=expression(log[2]~fold~change), main=paste0(cond1, " vs ", cond2, " MA Plot"))
garbage <- dev.off() # Save to file


# Sample Distance Matrix
sampleDists <- dist(t(assay(vsd)))
pdf(paste0(output, cond1, "_vs_", cond2, "-DistanceMatrix.pdf"))
pheatmap(as.matrix(sampleDists), col=rev(dmcol), clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames=TRUE, annotation_col=conditions_df, annotation_row=conditions_df, annotation_names_row=FALSE, annotation_names_col=FALSE, annotation_colors=condcol, legend = FALSE, main=paste0(cond1, " vs ", cond2, " Matrix"))
garbage <- dev.off() # Save to file


# Plot dispersions
pdf(paste0(output, cond1, "_vs_", cond2, "-DispersionPlot.pdf"))
plotDispEsts(dds, main=paste0(cond1, " vs ", cond2, " Dispersion Plot"))
garbage <- dev.off() # Save to file


# Principal Components Analysis
pdf(paste0(output, cond1, "_vs_", cond2, "-PCA.pdf"))
plotPCA(vsd, intgroup="condition")
garbage <- dev.off() # Save to file


# Volcano plot with "significant" genes labeled
# sigthresh=0.05 and lfcthresh=1 are the default values if ommitted when creating the graph. No need to edit here.
volcanoplot <- function (res, lfcthresh=1, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "Both"), pch=20, col=c("red","orange","green"))
}


# Custom Volcano Plot with significant gene names. Here we use the sigpadj and siglfc variable values to generate the graph
pdf(paste0(output, cond1, "_vs_", cond2, "-VolcanoPlot.pdf"))
volcanoplot(resdata, siglfc, sigpadj, textcx=.8, xlim=c(-2, 2), main=paste0(cond1, " vs ", cond2, " Volcano Plot"))
garbage <- dev.off() # Save to file


# Heatmap of significant hits ( padj<sigpadj and |log2FoldChange|>=siglfc )
sig <- rownames(res[!is.na(res\$padj) & res\$padj<sigpadj & abs(res\$log2FoldChange)>=siglfc,])[1:30]	# Select the first 30 significant hits (sorted by padj)
sig <- sig[!is.na(sig)]																					# If we have less than 30 hits, clean "NA" fields
counts <- counts(dds,normalized=TRUE)																	# Get normalized read counts (sorted by padj)
counts <- counts[apply(counts, 1, function(row) all(row !=0 )),]										# Remove genes with zero reads
sigcounts <- counts(dds,normalized=TRUE)[sig,]															# Get normalized read counts (sorted by padj) for significant genes
try({																									# May fail if only 0 or 1 gene is significant
sigcounts <- sigcounts[apply(sigcounts, 1, function(row) all(row !=0 )),]								# Remove genes with zero reads
})


# By defaults hits are not clustered and thus stay sorted by their padj value
pdf(paste0(output, cond1, "_vs_", cond2, "-HeatmapSig.pdf"), onefile=FALSE)
try({																									# May fail if only 0 or 1 gene is significant
pheatmap(log2(sigcounts), col=rev(hmcol), scale='row', cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=conditions_df, annotation_names_col=FALSE, annotation_colors=condcol, main=paste0(cond1, " vs ", cond2, " Top Hits"))
})
garbage <- dev.off() # Save to file


# But we can cluster them using the following command
pdf(paste0(output, cond1, "_vs_", cond2, "-HeatmapSigClust.pdf"), onefile=FALSE)
try({																									# May fail if only 0 or 1 gene is significant
pheatmap(log2(sigcounts), col=rev(hmcol), scale='row', cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=conditions_df, annotation_names_col=FALSE, annotation_colors=condcol, main=paste0(cond1, " vs ", cond2, " Clustered Top Hits"))
})
garbage <- dev.off() # Save to file


# And even plot every gene
pdf(paste0(output, cond1, "_vs_", cond2, "-HeatmapAllClust.pdf"), onefile=FALSE)
pheatmap(log2(counts), col=rev(hmcol), scale='row', cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames=FALSE, annotation_col=conditions_df, annotation_names_col=FALSE, annotation_colors=condcol, main=paste0(cond1, " vs ", cond2, " Clustered Reads Count"))
garbage <- dev.off() # Save to file


# Additional optional graphs using DESeqAnalysis https://github.com/acidgenomics/DESeqAnalysis

# Principal Components Analysis with labels
try(library(DESeqAnalysis))
pdf(paste0(output, cond1, "_vs_", cond2, "-PCAlabels.pdf"), width=20, height=20)
try(DESeqAnalysis::plotPCA(vsd, label=TRUE, interestingGroups="condition", labels=list(title=paste0(cond1, " vs ", cond2, " PCA"))))
garbage <- dev.off() # Save to file


# Adapted from:
# http://www.bioconductor.org/help/workflows/rnaseqGene/
# https://gist.github.com/stephenturner/f60c1934405c127f09a6
# https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
RDELIM
			# Use Rscript
			echo "Rscript ${dir2}/DESeq2${reanalysis}/${Rscript} || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

			# Cleaning commands
			if [ "${debug}" != "1" ]
			then
				# Remove .sbatch
				echo "rm ${dir2}/${job}_${samplename}.sbatch"
			fi
			echo "exit 0"

		} > "${dir2}/${job}_${samplename}.sbatch"

		# Queue job
		SBdeseq2=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
		if [ -n "${SBdeseq2##* }" ]
		then
			SBdeseq2IDs=${SBdeseq2IDs}:${SBdeseq2##* }
		fi
		echo -e "\t ${samplename} DESeq2 job queued"

	done < "${dir2}"/matrices

fi



if [ "${diffexp}" = "cufflinks" ] || [ "${diffexp}" = "both" ]
then

	echo ""
	echo "-- Cufflinks --"



	######################
	## Assemble transcriptomes

	# Will only run in reanalysis mode if the .cuff folder was not previously created during the initial run.
	# Variables
	job="cufflinks"
	firstsample="1"

	while read -r line;
	do

		# General variables
		read1=$(echo "${line}" | cut -f1)
		samplename=$(echo "${read1}" | awk -F_R1 '{print $1}')

		if [ -n "${reanalysis}" ] && [ ! -d "${dir2}/${samplename}.cuff" ]
		then

			if [ "${firstsample}" = "1" ]
			then
				echo ""
				echo -e "\t-- Queuing missing cufflinks jobs"
				firstsample="0"
			fi

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
				echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
				echo "#SBATCH --requeue"
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
				# No previous jobs expected

				# General commands
				echo "mkdir -p \"${tmp}\""
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				# Create a sample specific folder, and as cufflinks do not allow to specify an output directory then cd to it.
				echo "mkdir -p ${dir2}/${samplename}.cuff"
				echo "cd ${dir2}/${samplename}.cuff || exit"
				# Create cufflinks job
				echo "cufflinks -u -p ${threads} -g ${gtf} -b ${fasta_refgenome} ${dir2}/${bamsortedout} || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/${job}_${samplename}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/${job}_${samplename}.sbatch"

			# Queue job
			SBcufflinks=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			SBcufflinksIDs=${SBcufflinksIDs}:${SBcufflinks##* }
			echo -e "\t ${samplename} missing Cufflinks job queued"
		fi
	done < "${dir2}/Fastqs"



	if [ -z "${reanalysis}" ] || [ "${firstsample}" = "0" ]
	then

		######################
		## Merge transcripts

		# Variables
		job="cuffmerge"
		# We are no longer processing sample files, but GTF files
		samplename="GTFs"

		{
			# General SLURM parameters
			echo '#!/bin/bash'
			echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
			echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
			echo "#SBATCH --requeue"
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
			if [ -n "${SBcufflinksIDs}" ]
			then
				echo "#SBATCH --dependency=afterok${SBcufflinksIDs}"
			fi

			# General commands
			echo "mkdir -p \"${tmp}\""
			if [ -n "${customcmd}" ]
			then
				echo "${customcmd}"
			fi

			# Job specific commands
			# Generate a list of sample specific transcripts.gtf files to be merged
			echo "find ${dir2}/*.cuff/ -maxdepth 1 -type f -name 'transcripts.gtf' > ${dir2}/assembly_GTF_list.txt"
			echo "cuffmerge -p ${threads} -g ${gtf} -s ${fasta_refgenome} -o ${dir2} ${dir2}/assembly_GTF_list.txt || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

			# Cleaning commands
			if [ "${debug}" != "1" ]
			then
				echo "[ -f ${dir2}/assembly_GTF_list.txt ] && rm ${dir2}/assembly_GTF_list.txt"
				# Remove .sbatch
				echo "rm ${dir2}/${job}_${samplename}.sbatch"
			fi
			echo "exit 0"

		} > "${dir2}/${job}_${samplename}.sbatch"

		# Queue job
		SBcuffmerge=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
		echo ""
		echo -e "\t-- Cuffmerge job queued"



		######################
		## Quantify transcripts

		job="cuffquant"

		echo ""
		echo -e "\t-- Queuing quantification jobs"
		echo ""

		while read -r line;
		do
			# Variables
			read1=$(echo "${line}" | cut -f1)
			read2=$(echo "${line}" | cut -f2)
			samplename=$(echo "${read1}" | awk -F_R1 '{print $1}')
			bamsortedout=$(basename "${read1}" | sed "s/_R1${fileext}/.sorted.bam/g")

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
				echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
				echo "#SBATCH --requeue"
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
				echo "#SBATCH --dependency=afterok:${SBcuffmerge##* }"

				# General commands
				echo "mkdir -p \"${tmp}\""
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				echo "cuffquant -p ${threads} -o ${dir2}/${samplename}.cuff ${dir2}/merged.gtf ${dir2}/${bamsortedout} || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/${job}_${samplename}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/${job}_${samplename}.sbatch"

			# Queue job
			SBcuffquant=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			SBcuffquantIDs=${SBcuffquantIDs}:${SBcuffquant##* }
			echo -e "\t ${samplename} cuffquant job queued"

		done < "${dir2}/Fastqs"

	fi



	######################
	## Compare expression levels

	job="cuffdiff"

	# Variables
	# We are no longer processing sample files, but CXBs files
	samplename="CXBs"

	{
		# General SLURM parameters
		echo '#!/bin/bash'
		echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
		echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
		echo "#SBATCH --requeue"
		echo "#SBATCH --time=2:00:00"
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
		if [ -z "${reanalysis}" ]
		then
			echo "#SBATCH --dependency=afterok${SBcuffquantIDs}"
		fi

		# General commands
		echo "mkdir -p \"${tmp}\""
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}"
		fi

		# Job specific commands
		echo "cuffdiff -u -p ${threads} -b ${fasta_refgenome} ${dir2}/merged.gtf ${dir2}/$(echo "${groupedsamples}" | sed "s#,#.cuff/abundances.cxb,${dir2}/#g;s# #.cuff/abundances.cxb ${dir2}/#g").cuff/abundances.cxb -L ${labels} -o ${dir2}/Cufflinks${reanalysis} || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

		# Cleaning commands
		if [ "${debug}" != "1" ]
		then
			# Remove .sbatch
			echo "rm ${dir2}/${job}_${samplename}.sbatch"
		fi
		echo "exit 0"

	} > "${dir2}/${job}_${samplename}.sbatch"

	# Queue job
	SBcuffdiff=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
	echo ""
	echo -e "\t-- Cuffdiff job queued"



	######################
	## Normalize expression levels

	job="cuffnorm"

	{
		# General SLURM parameters
		echo '#!/bin/bash'
		echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
		echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
		echo "#SBATCH --requeue"
		echo "#SBATCH --time=1:00:00"
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
		if [ -z "${reanalysis}" ]
		then
			echo "#SBATCH --dependency=afterok${SBcuffquantIDs}"
		fi

		# General commands
		echo "mkdir -p \"${tmp}\""
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}"
		fi

		# Job specific commands
		echo "cuffnorm -p ${threads} ${dir2}/merged.gtf ${dir2}/$(echo "${groupedsamples}" | sed "s#,#.cuff/abundances.cxb,${dir2}/#g;s# #.cuff/abundances.cxb ${dir2}/#g").cuff/abundances.cxb -L ${labels} -o ${dir2}/Cufflinks${reanalysis}/Normalized || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

		# Cleaning commands
		if [ "${debug}" != "1" ]
		then
			# Remove .sbatch
			echo "rm ${dir2}/${job}_${samplename}.sbatch"
		fi
		echo "exit 0"

	} > "${dir2}/${job}_${samplename}.sbatch"

	# Queue job
	SBcuffnorm=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
	echo ""
	echo -e "\t-- Cuffnorm job queued"

fi



######################
## MAJIQ Splicing variant analysis

if [ "${splicing}" = "majiq" ] || [ "${splicing}" = "both" ]
then

	echo ""
	echo "-- MAJIQ --"

	mkdir -p "${dir2}/MAJIQ${reanalysis}"



	######################
	## Generate config file

	# Variables
	job="majiq"
	IFS="," read -r -a labelsarray <<< "${labels}"
	IFS=" " read -r -a samplesarray <<< "${groupedsamples}"
	declare -a replicatesarray

	cat > "${dir2}/MAJIQ${reanalysis}/config" <<MBDELIM
[info]
bamdirs=${dir2}
genome=${genomebuild}
strandness=${strandness}
[experiments]
MBDELIM

	# Generate an array with sample replicate names
	repN="0"
	while [ ${repN} -lt ${#labelsarray[@]} ]
	do
		replicatesarray[${repN}]="${samplesarray[${repN}]//,/.sorted,}"
		repN=$((repN + 1))
	done

	# Generate the majiq "experiments" lines from the config file
	condA="0"
	while [ ${condA} -lt ${#labelsarray[@]} ]
	do
		echo "${labelsarray[${condA}]}=${replicatesarray[${condA}]}.sorted" >> "${dir2}/MAJIQ${reanalysis}/config"
		condA=$((condA + 1))
	done



	######################
	## LSV detection using majiq build

	# Variables
	# We are no longer processing sample files, but doing majiq
	samplename="build"

	{
		# General SLURM parameters
		echo '#!/bin/bash'
		echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
		echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
		echo "#SBATCH --requeue"
		echo "#SBATCH --time=2:00:00"
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
		if [ -z "${reanalysis}" ]
		then
			echo "#SBATCH --dependency=afterok${SBalignIDs}"
		fi

		# General commands
		echo "mkdir -p \"${tmp}\""
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}"
		fi

		# Job specific commands
		# MAJIQ build
		echo "majiq build ${gff3} -c ${dir2}/MAJIQ${reanalysis}/config -j ${threads} -o ${dir2}/MAJIQ${reanalysis}/ --min-intronic-cov 1 --simplify || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

		# Cleaning commands
		if [ "${debug}" != "1" ]
		then
			echo "rm ${dir2}/MAJIQ${reanalysis}/majiq.log"
			# Keeping .sj files by default since MAJIQ 2.2 for incremental build compatibility.
			#echo "rm ${dir2}/MAJIQ${reanalysis}/*.sorted.sj"
			# Remove .sbatch
			echo "rm ${dir2}/${job}_${samplename}.sbatch"
		fi
		echo "exit 0"

	} > "${dir2}/${job}_${samplename}.sbatch"

	# Queue job
	SBmajiqbuild=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
	echo ""
	echo -e "\t-- Build job queued"



	######################
	## PSI quantification

	if [ "${psi}" = "1" ]
	then

		echo ""
		echo -e "\t-- Queuing PSI jobs"
		echo ""

		# Variables
		job="majiq-psi"
		IFS=" " read -r -a samplesarray <<< "${groupedsamples//,/ }"

		# Perform PSI quantification on each individual sample
		condA="0"
		while [ ${condA} -lt ${#samplesarray[@]} ]
		do
			# Variables
			samplename="${samplesarray[${condA}]}"

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
				echo "#SBATCH $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi)"
				echo "#SBATCH --requeue"
				echo "#SBATCH --time=30:00"
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
				echo "#SBATCH --dependency=afterok:${SBmajiqbuild##* }"

				# General commands
				echo "mkdir -p \"${tmp}\""
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				# MAJIQ PSI analysis
				echo "majiq psi ${dir2}/MAJIQ${reanalysis}/${samplesarray[${condA}]}.sorted.majiq -j $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "1"; else echo "${threads}"; fi) -o ${dir2}/MAJIQ${reanalysis}/ -n ${samplesarray[${condA}]} || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/${job}_${samplename}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/${job}_${samplename}.sbatch"

			# Queue job
			SBmajiqpsi=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			SBmajiqpsiIDs=${SBmajiqpsiIDs}:${SBmajiqpsi##* }
			echo -e "\t ${samplesarray[${condA}]} PSI job queued"

			condA=$((condA + 1))
		done
	fi



	######################
	## Delta PSI quantification

	echo ""
	echo -e "\t-- Queuing Delta PSI jobs"
	echo ""

	# Variables
	job="majiq-Dpsi"
	IFS="," read -r -a labelsarray <<< "${labels}"
	IFS=" " read -r -a samplesarray <<< "${groupedsamples}"
	declare -a replicatesarray

	# Generate an array with sample replicate names
	repN="0"
	while [ ${repN} -lt ${#labelsarray[@]} ]
	do
		replicatesarray[${repN}]="${samplesarray[${repN}]//,/.sorted.majiq ${dir2}/MAJIQ${reanalysis}/}"
		repN=$((repN + 1))
	done

	# Generate pair-wise comparison matrices
	condA="0"
	condB="0"
	while [ ${condA} -lt ${#labelsarray[@]} ]
	do
		condB=$((condA + 1))
		while [ ${condB} -lt ${#labelsarray[@]} ]
		do
			# Variables
			samplename="${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}"

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
				echo "#SBATCH $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi)"
				echo "#SBATCH --requeue"
				echo "#SBATCH --time=30:00"
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
				echo "#SBATCH --dependency=afterok:${SBmajiqbuild##* }"

				# General commands
				echo "mkdir -p \"${tmp}\""
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				# By default we only use 1 cpu, easier to be backfilled, change if needed.
				# MAJIQ Delta PSI analysis
				# Job likely to fail for some of the odd pair-wise comparisons. This way will fail silently and remove the .sbatch and the empty folder.
				echo "majiq deltapsi -grp1 ${dir2}/MAJIQ${reanalysis}/${replicatesarray[${condA}]}.sorted.majiq -grp2 ${dir2}/MAJIQ${reanalysis}/${replicatesarray[${condB}]}.sorted.majiq -j $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "1"; else echo "${threads}"; fi) -o ${dir2}/MAJIQ${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}/ -n ${labelsarray[${condA}]} ${labelsarray[${condB}]} || { rm -rf ${dir2}/MAJIQ${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}/ ${dir2}/${job}_${samplename}.sbatch; exit 0; }"
				# Generate corresponding VOILA tsv files with default filtering options
				echo "voila tsv -j $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "1"; else echo "${threads}"; fi) ${dir2}/MAJIQ${reanalysis}/splicegraph.sql ${dir2}/MAJIQ${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}/${labelsarray[${condA}]}_${labelsarray[${condB}]}.deltapsi.voila -f ${dir2}/MAJIQ${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}.filtered.tsv || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					echo "rm ${dir2}/MAJIQ${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}/deltapsi_majiq.log"
					# Remove .sbatch
					echo "rm ${dir2}/${job}_${samplename}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/${job}_${samplename}.sbatch"

			# Queue job
			SBmajiqdeltapsi=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			SBmajiqdeltapsiIDs=${SBmajiqdeltapsiIDs}:${SBmajiqdeltapsi##* }
			echo -e "\t ${labelsarray[${condA}]}_vs_${labelsarray[${condB}]} Delta PSI job queued"

			condB=$((condB + 1))
		done

		condA=$((condA + 1))
	done

fi



######################
## LeafCutter Splicing variant analysis

if [ "${splicing}" = "leafcutter" ] || [ "${splicing}" = "both" ]
then

	echo ""
	echo "-- LeafCutter --"

	mkdir -p "${dir2}/LeafCutter${reanalysis}"



	######################
	## LeafCutter convert bams to juncs

	# Will only run in reanalysis mode if the .junc file was not previously created during the initial run.
	# Variables
	job="leaf-junc"
	firstsample="1"

	while read -r line;
	do

		# General variables
		read1=$(echo "${line}" | cut -f1)
		samplename=$(echo "${read1}" | awk -F_R1 '{print $1}')

		if [ -n "${reanalysis}" ] && [ ! -f "${dir2}/${samplename}.junc" ]
		then

			if [ "${firstsample}" = "1" ]
			then
				echo ""
				echo -e "\t-- Queuing missing leafcutter jobs"
				firstsample="0"
			fi

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
				echo "#SBATCH $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi)"
				echo "#SBATCH --requeue"
				echo "#SBATCH --time=2:00:00"
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
				# No previous jobs expected

				# General commands
				echo "mkdir -p \"${tmp}\""
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				# adapted from https://github.com/davidaknowles/leafcutter/blob/master/scripts/bam2junc.sh
				echo "samtools view ${dir2}/${samplename}.sorted.bam | ${leafCutterDir}/scripts/filter_cs.py | ${leafCutterDir}/scripts/sam2bed.pl --use-RNA-strand - ${dir2}/${samplename}.bed || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
				echo "${leafCutterDir}/scripts/bed2junc.pl ${dir2}/${samplename}.bed ${dir2}/${samplename}.junc || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					echo "rm ${dir2}/${samplename}.bed"
					# Remove .sbatch
					echo "rm ${dir2}/${job}_${samplename}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/${job}_${samplename}.sbatch"

			# Queue job
			SBleafjunc=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			SBleafjuncIDs=${SBleafjuncIDs}:${SBleafjunc##* }
			echo -e "\t ${samplename} missing JUNC conversion job queued"
		fi
	done < "${dir2}/Fastqs"



	######################
	## Intron clustering

	if [ -z "${reanalysis}" ] || [ "${firstsample}" = "0" ]
	then

		# Variables
		job="leaf-clust"
		# We are no longer processing sample files, but doing intron clustering
		samplename="introns"

		{
			# General SLURM parameters
			echo '#!/bin/bash'
			echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
			echo "#SBATCH $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi)"
			echo "#SBATCH --requeue"
			echo "#SBATCH --time=30:00"
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
			if [ -n "${SBleafjuncIDs}" ]
			then
				echo "#SBATCH --dependency=afterok${SBleafjuncIDs}"
			fi


			# General commands
			echo "mkdir -p \"${tmp}\""
			if [ -n "${customcmd}" ]
			then
				echo "${customcmd}"
			fi

			# Job specific commands
			echo "find -L ${dir2} -maxdepth 1 -name '*.junc' > ${dir2}/juncfiles.txt"
			echo "python ${leafCutterDir}/clustering/leafcutter_cluster.py -j ${dir2}/juncfiles.txt -r ${dir2}/ || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

			# Cleaning commands
			if [ "${debug}" != "1" ]
			then
				echo "rm ${dir2}/juncfiles.txt"
				echo "rm ${dir2}/*.junc.leafcutter.sorted.gz"
				echo "rm ${dir2}/leafcutter_perind.counts.gz"
				echo "rm ${dir2}/leafcutter_pooled"
				echo "rm ${dir2}/leafcutter_refined"
				echo "rm ${dir2}/leafcutter_sortedlibs"
				# Remove .sbatch
				echo "rm ${dir2}/${job}_${samplename}.sbatch"
			fi
			echo "exit 0"

		} > "${dir2}/${job}_${samplename}.sbatch"

		# Queue job
		SBleafintclust=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
		echo ""
		echo -e "\t-- Intron clustering job queued"
	fi



	######################
	## Grouped differential intron excision analysis, plot and leafviz prep jobs for faster processing

	echo ""
	echo -e "\t-- Queuing differential analysis, plotting and leafviz prep jobs"
	echo ""

	# Variables
	job="leaf-diffex"
	IFS="," read -r -a labelsarray <<< "${labels}"
	IFS=" " read -r -a samplesarray <<< "${groupedsamples}"
	condA="0"
	condB="0"

	while [ ${condA} -lt ${#labelsarray[@]} ]
	do
		condB=$((condA + 1))
		while [ ${condB} -lt ${#labelsarray[@]} ]
		do
			# Variables
			samplename="${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}"

			# Generate pair-wise comparison groups_file files
			echo -e "${samplesarray[${condA}]//,/\t${labelsarray[${condA}]}\n}\t${labelsarray[${condA}]}" > "${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_group.txt"
			echo -e "${samplesarray[${condB}]//,/\t${labelsarray[${condB}]}\n}\t${labelsarray[${condB}]}" >> "${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_group.txt"

			# Check number of replicates in each groups to fine tune -i and -g on leafcutter_ds.R is less than 4 samples
			IFS="," read -r -a replicatesA <<< "${samplesarray[${condA}]}"
			IFS="," read -r -a replicatesB <<< "${samplesarray[${condB}]}"
			if [ ${#replicatesA[@]} -lt ${#replicatesB[@]} ]
			then
				numsamp="${#replicatesA[@]}"
			else
				numsamp="${#replicatesB[@]}"
			fi

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
				echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
				echo "#SBATCH --requeue"
				echo "#SBATCH --time=1:00:00"
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
				if [ -n "${SBleafintclust##* }" ]
				then
					echo "#SBATCH --dependency=afterok:${SBleafintclust##* }"
				fi

				# General commands
				echo "mkdir -p \"${tmp}\""
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				# Differential intron excision analysis
				echo "${leafCutterDir}/scripts/leafcutter_ds.R -p ${threads} -o ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]} $(if [ -n "${exonfile}" ]; then echo "-e ${exonfile}"; fi) $(if [ "${numsamp}" -le 4 ]; then echo "-i ${numsamp} -g ${numsamp}"; fi) ${dir2}/leafcutter_perind_numers.counts.gz ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_group.txt || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"
				# Prepare results for visualisation
				# Job will fail if there are no significant clusters, instead silently exit. Change the FDR cutoff and run again using the ShinyLeafviz script.
				echo "${leafCutterDir}/leafviz/prepare_results.R -f 0.05 -m ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_group.txt -c ${labelsarray[${condA}]}_vs_${labelsarray[${condB}]} -o ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}.RData ${dir2}/leafcutter_perind_numers.counts.gz ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_cluster_significance.txt ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_effect_sizes.txt ${anncode} || { rm ${dir2}/${job}_${samplename}.sbatch; exit 0; }"
				# Plotting splice junctions
				echo "${leafCutterDir}/scripts/ds_plots.R -f 0.05 -o ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}.pdf $(if [ -n "${exonfile}" ]; then echo "-e ${exonfile}"; fi) ${dir2}/leafcutter_perind_numers.counts.gz ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_group.txt ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_cluster_significance.txt || if [ -f ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err ]; then echo \"\${SLURM_JOB_NODELIST}\" >> ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && exit 1; else echo \"\${SLURM_JOB_NODELIST}\" > ${dir2}/\"\${SLURM_JOBID}\"-${job}_${samplename}.err && scontrol requeue \"\${SLURM_JOBID}\" && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/${job}_${samplename}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/${job}_${samplename}.sbatch"

			# Queue job
			SBleafdiffex=$(until sbatch "${dir2}/${job}_${samplename}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			SBleafdiffexIDs=${SBleafdiffexIDs}:${SBleafdiffex##* }
			echo -e "\t ${labelsarray[${condA}]}_vs_${labelsarray[${condB}]} analysis job queued"

			# Create an easy to launch script for the Shiny app, still need to deal with port forfarding on your own
			{
				echo '#!/bin/bash'
				echo "${customcmd}"
				echo ""
				echo 'leafFDR="0.05"'
				echo ""
				echo '# Uncomment to plot splice junctions with a different adjusted p value threshold'
				echo "#${leafCutterDir}/scripts/ds_plots.R -f \${leafFDR} -o ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}.pdf $(if [ -n "${exonfile}" ]; then echo "-e ${exonfile}"; fi) ${dir2}/leafcutter_perind_numers.counts.gz ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_group.txt ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_cluster_significance.txt"
				echo ""
				echo '# Uncomment to prepare results for visualisation with a different adjusted p value threshold'
				echo "#${leafCutterDir}/leafviz/prepare_results.R -f \${leafFDR} -m ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_group.txt -c ${labelsarray[${condA}]}_vs_${labelsarray[${condB}]} -o ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}.RData ${dir2}/leafcutter_perind_numers.counts.gz ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_cluster_significance.txt ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_effect_sizes.txt ${anncode} || exit 0"
				echo ""
				echo '# Actual command to prepare results for visualisation'
				echo "cd ${leafCutterDir}/leafviz/"
				echo "echo 'Launching leafviz: you may need to have to deal with port forwarding on your own to access the website'"
				echo "./run_leafviz.R ${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}.RData"
			} > "${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_ShinyLeafviz.sh"
			chmod +x "${dir2}/LeafCutter${reanalysis}/${labelsarray[${condA}]}_vs_${labelsarray[${condB}]}_ShinyLeafviz.sh"

			condB=$((condB + 1))
		done

		condA=$((condA + 1))
	done

fi



######################
## Final processing of result files

job="notif"

# Variables
# We are no longer processing sample files, using project name
samplename="$(basename "${dir2}")-RNAseq"

{
	# General SLURM parameters
	echo '#!/bin/bash'
	echo "#SBATCH --job-name=${job}_${samplename} --output=${dir2}/${logs}/${job}_${samplename}.out --error=${dir2}/${logs}/${job}_${samplename}.err --open-mode=append"
	echo "#SBATCH $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi)"
	echo "#SBATCH --requeue"
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
	echo "#SBATCH --dependency=afterok$(if [ -n "${SBcuffdiff##* }" ]; then echo ":${SBcuffdiff##* }"; fi)$(if [ -n "${SBcuffnorm##* }" ]; then echo ":${SBcuffnorm##* }"; fi)$(if [ -n "${SBfqcIDs}" ]; then echo "${SBfqcIDs}"; fi)$(if [ -n "${SBdeseq2IDs}" ]; then echo "${SBdeseq2IDs}"; fi)$(if [ -n "${SBmajiqpsiIDs}" ]; then echo "${SBmajiqpsiIDs}"; fi)$(if [ -n "${SBmajiqdeltapsiIDs}" ]; then echo "${SBmajiqdeltapsiIDs}"; fi)$(if [ -n "${SBleafdiffexIDs}" ]; then echo "${SBleafdiffexIDs}"; fi)"

	# General commands
	echo "mkdir -p \"${tmp}\""
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}"
	fi

	# Job specific commands
	# Creating folders for archive
	echo "[ -f ${dir2}/Results_$(basename "${dir2}").tar.gz ] && rm ${dir2}/Results_$(basename "${dir2}").tar.gz"
	echo "[ -d ${dir2}/Results_$(basename "${dir2}") ] && rm -rf ${dir2}/Results_$(basename "${dir2}")"
	echo "mkdir -p ${dir2}/Results_$(basename "${dir2}")/"
	# List all csv, html, pdf, results and zip files
	echo "csvarchive=\$(ls ${dir2}/*.csv)"			# Various
	echo "htmlarchive=\$(ls ${dir2}/*.html)"		# FastQC
	echo "pdfarchive=\$(ls ${dir2}/*.pdf)"			# Various
	echo "resultsarchive=\$(ls ${dir2}/*.results)"	# RSEM files
	echo "ziparchive=\$(ls ${dir2}/*.zip)"			# FastQC
	# Create an archive of all csv, html, pdf and zip files for easy download
	echo "for i in \${csvarchive}; do cp -rf \"\${i}\" ${dir2}/Results_$(basename "${dir2}")/\"\$(basename \"\${i}\")\"; done"
	echo "for i in \${htmlarchive}; do cp -rf \"\${i}\" ${dir2}/Results_$(basename "${dir2}")/\"\$(basename \"\${i}\")\"; done"
	echo "for i in \${pdfarchive}; do cp -rf \"\${i}\" ${dir2}/Results_$(basename "${dir2}")/\"\$(basename \"\${i}\")\"; done"
	echo "for i in \${resultsarchive}; do cp -rf \"\${i}\" ${dir2}/Results_$(basename "${dir2}")/\"\$(basename \"\${i}\")\"; done"
	echo "for i in \${ziparchive}; do cp -rf \"\${i}\" ${dir2}/Results_$(basename "${dir2}")/\"\$(basename \"\${i}\")\"; done"
	echo "cp -rf ${dir2}/Cufflinks${reanalysis} ${dir2}/Results_$(basename "${dir2}")/"
	echo "cp -rf ${dir2}/merged.gtf ${dir2}/Results_$(basename "${dir2}")/"
	echo "cp -rf ${dir2}/DESeq2${reanalysis} ${dir2}/Results_$(basename "${dir2}")/"
	echo "cp -rf ${dir2}/MAJIQ${reanalysis} ${dir2}/Results_$(basename "${dir2}")/"
	echo "rm \$(find ${dir2}/Results_$(basename "${dir2}")/MAJIQ${reanalysis} -maxdepth 1 -type f -not -name '*.sql')"
	echo "cp -rf ${dir2}/LeafCutter${reanalysis} ${dir2}/Results_$(basename "${dir2}")/"
	echo "tar --remove-files -C ${dir2} -pczf ${dir2}/Results_$(basename "${dir2}").tar.gz Results_$(basename "${dir2}")"
	if [ -n "${iftttkey}" ] && [ -n "${notifevent}" ]
	then
		# Trigger IFTTT maker channel event when it's ready, nice isn't it?
		echo "curl -X POST -H \"Content-Type: application/json\" -d '{ \"value1\" : \"$(basename "${dir2}")\" , \"value2\" : \"$(whoami)\"}' https://maker.ifttt.com/trigger/${notifevent}/with/key/${iftttkey}"
	fi
	if [ -n "${slackwebhookURL}" ]
	then
		# Trigger a Slack Webhook when it's ready, nice isn't it?
		echo "curl -X POST -d \"payload={'blocks': [ { 'type': 'section','text': {'type': 'mrkdwn','text': '*${notifevent}\$(date +' on %A %B %d %Y at %H:%M')*' } }, {'type': 'section', 'text': {'type': 'mrkdwn','text': '@$(whoami) - $(basename "${dir2}") is done'} } ] }\" ${slackwebhookURL}"
	fi

	# Cleaning commands
	if [ "${debug}" != "1" ]
	then
		# Remove error files upon successful completion. Comment to disable.
		echo "rm ${dir2}/*.err"
		# Remove logs folder upon successfull completion. Comment to disable.
		echo "rm -rf ${dir2}/logs"
		# Remove Temporary directory
		echo "rm -rf \"${tmp}\""
		# Remove .sbatch
		echo "rm ${dir2}/${job}_${samplename}.sbatch"
	fi
	echo "exit 0"

} > "${dir2}/${job}_${samplename}.sbatch"

# Queue job
until sbatch "${dir2}/${job}_${samplename}.sbatch" >/dev/null 2>&1; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done
echo ""
echo "-- Results ready notification job queued --"
echo ""



# Clean temporary files
if [ "${debug}" != "1" ]
then
	[ -f "${dir2}/files1" ] && rm "${dir2}/files1"
	[ -f "${dir2}/files2" ] && rm "${dir2}/files2"
	[ -f "${dir2}/Fastqs" ] && rm "${dir2}/Fastqs"
	[ -f "${dir2}/matrices" ] && rm "${dir2}/matrices"
fi

# That's all folks!
exit 0
