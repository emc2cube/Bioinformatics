#!/bin/bash
# shellcheck disable=SC2028,SC2030,SC2031
#
# Usage: sh_CRISPR.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]
#
##############################################################
##                       Description                        ##
##############################################################
#
# This script will process the fastq(.gz) files generated in
# a typical CRISPR screen using either casTLE or MAGeCK.
#
# If using casTLE, a reference file of all the Indices will be
# automatically created using bowtie (NOT bowtie2).
# It will then use casTLE scripts to analyze the screen and
# generate basic graphs
# Download casTLE from https://bitbucket.org/dmorgens/castle/
#
# If using MAGeCK counts, tests, mle and pathway analysis
# will be performed. It will also run the R package MAGeCKFlute
# and in all cases generate basic graphs.
# Download MAGeCK from https://sourceforge.net/projects/mageck/
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
## Screen options
#
# day0-label.
# Specify the label for control sample (usually day 0 or plasmid).
# If using MAGeCK it will also turn on the negative selection QC for every other sample label.
# The negative selection QC will compare each other sample with day0 sample, and thus estimate the degree of negative selections in essential genes.
# If using casTLE, only the first TWO replicates will be used.
# day0="Plasmid,Control_t0_rep1,Control_t0_rep2,Treated_t0_rep1,Treated_t0_rep2"
day0=""
#
# Tests groups.
# Enter your sample names (not including ".fastq" or ".fastq.gz") for comparisons.
# Delimit replicate samples by a COMMA and groups by a SPACE.
# If using casTLE, only the first TWO replicates will be used.
# testgroups="Control_rep1,Control_rep2 Treated_condA_rep1,Treated_condA_rep Treated_condB_rep1,Treated_condB_rep2"
testgroups=""
#
#
## casTLE options
# download the last version from https://bitbucket.org/dmorgens/castle/
#
# Use casTLE.
# 0 = No ; 1 = Yes
usecastle="0"
#
# Python version.
# Are you using python2.7 (original) or python3 (included in this repo) casTLE scripts?
# If you also use MAGeCK python3 is REQUIRED
# use "python" or "python3"
python="python3"
#
# casTLE folder location.
castlepath="/Tools/dmorgens-castle"
#
# Number of permutations to generate p-values.
# For a first pass, use 5x the permutations (so for 20,000 genes that is 100,000 permutations).
# For publication, use 50x (Default, takes time).
permutations="1000000"
#
# Perform permutations on the individual result files.
# Permutations are always calculated for the combo file.
# 0 = No ; 1 = Yes
permres="0"
#
# Output format for graphs.
# If left empty default value is png.
# (png, pdf, eps)
graphformat="pdf"
#
# Screen performed in mouse cells.
# Default is 0, for human cells.
# 0 = No ; 1 = Yes
mouse="0"
#
# Path to bowtie, bowtie-build needs to be in the same folder (probably is already the case).
# bowtie="/usr/local/bin/bowtie"
bowtie="/usr/local/bin/bowtie"
#
# Type of screen. Will be used to create Indices for the guides.
screentype="Cas9-10"
#
# Name of the guides index file. Will be saved in the Indices folder.
# It will overwrite any files with this name prefix.
outputbowtieindex=""
#
# Oligo file location.
# Leave empty if it was previously used and the corresponding Index are already generated for this type of screen.
oligofile=""
#
#
## MAGeCK options.
# download the last version from https://sourceforge.net/projects/mageck/
#
# Use MAGeCK.
# 0 = No ; 1 = Yes
usemageck="1"
#
# MAGeCK list of sgRNA names (see https://sourceforge.net/p/mageck/wiki/input/#sgrna-library-file ) location.
magecksgRNAlibrary=""
#
# Use the reverse complement of the MAGeCK list of sgRNA names.
# 0 = No ; 1 = Yes
mageckrevcomplib="0"
#
# MAGeCK list of control sgRNA names (see https://sourceforge.net/p/mageck/wiki/input/#negative-control-sgrna-list ) location.
mageckcontrolsgrna=""
#
# GMT file for MAGeCK pathway analysis (see https://sourceforge.net/p/mageck/wiki/input/#pathway-file-gmt ) location.
gmtfile=""
#
# Matrix file for mle analysis ( see https://sourceforge.net/p/mageck/wiki/input/#design-matrix-file ) location.
# By default will look for a "matrix.txt" file stored with the FastQ files.
matrixfile=$([ -f "${1}"/matrix.txt ] && echo "${1}/matrix.txt")
#
#
## Notifications options
#
# Event name used in your webhooks or IFTTT recipes.
# The IFTTT maker channel will look for the combination private key + event name to then trigger your recipe.
# You can then create a recipe to send an email, a text message or a push notification.
# In the case of a Slack webhook this will be used in the message title
notifevent="CRISPR"
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
if [ "${1}" == "--help" ] || [ "${2}" == "--help" ] || [ "${3}" == "--help" ] || [ "${4}" == "--help" ]
then
	echo "Usage: $(basename "${0}") </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]"
	echo ""
	echo "Description"
	echo ""
	echo "This script will process the fastq(.gz) files generated in"
	echo "a typical CRISPR screen using either casTLE or MAGeCK."
	echo ""
	echo "If using casTLE a reference file of all the Indices will be"
	echo "automatically created using bowtie (NOT bowtie2)."
	echo "It will then use casTLE scripts to analyze the screen and"
	echo "generate basic graphs"
	echo "Download casTLE from https://bitbucket.org/dmorgens/castle/"
	echo ""
	echo "If using MAGeCK counts, tests, mle and pathway analysis"
	echo "will be performed. It will also run the R package MAGeCKFlute"
	echo "and in all cases generate basic graphs."
	echo "Download MAGeCK from https://sourceforge.net/projects/mageck/"
	echo ""
	echo "Options:"
	echo "$(basename "${0}") --help : Display this help message."
	echo "$(basename "${0}") --version : Display version number."
	echo ""
	exit
fi

# Version
if [ "${1}" == "--version" ] || [ "${2}" == "--version" ] || [ "${3}" == "--version" ] || [ "${4}" == "--version" ]
then
	echo "$(basename "${0}") version 2.0.2"
	echo "Major code cleaning (2.0.1)"
	echo "casTLE and MAGeCK support (2.0)"
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

if [ -z "${usecastle}" ] || [ -z "${usemageck}" ]
then
	if [ -z "${usecastle}" ]
	then
		usecastle="0"
	else
		usemageck="0"
	fi
fi

if [ ${usecastle} -eq "0" ] && [ ${usemageck} -eq "0" ]
then
	echo "Neither casTLE or MAGeCK are set up to be used right now, edit your config file."
	echo ""
	exit
fi

if [ ${usecastle} -eq "1" ]  ## casTLE checks
then
	if [ -n "${oligofile}" ]
	then
		if [ "${oligofile: -4}" != ".csv"  ]
		then
			echo "Invalid population oligo file file detected. Is ${oligofile} a .csv file?"
			echo ""
			exit
		fi
	fi
	if [ -n "${testgroups}" ]
	then
		IFS=" " read -r -a testsarray <<< "${day0} ${testgroups}"
		for replicates in "${testsarray[@]}"
		do
			IFS="," read -r -a  replicatesarray <<< "${replicates}"
			if [ "${#replicatesarray[@]}" -gt "2" ]
			then
				echo ""
				echo "WARNING"
				echo "This casTLE pipeline is only designed to work with a maximum of TWO replicates, only the first 2 replicates will be used"
				echo ""
			fi
		done
	else
		echo "You need to provide groups so casTLE can perform group by group comparisons."
		echo ""
		exit
	fi
fi

if [ ${usemageck} -eq "1" ]  ## MAGeCK checks
then
	if [ -z "${magecksgRNAlibrary}" ]
	then
		echo "Providing MAGeCK with a list of sgRNA names and sequences is mandatory, edit your config file."
		echo ""
		exit
	fi
	if [ -z "${testgroups}" ]
	then
		echo "You need to provide groups so MAGeCK can perform group by group comparisons."
		echo ""
		exit
	fi
	if [ "${python}" != "python3" ]
	then
		echo "MAGeCK require python 3, please install and rerun or disable MAGeCK in your config file."
		echo ""
		exit
	fi
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
echo "-- Hardware"
echo ""
echo "Up to ${threads} CPUs will be used"
echo "Up to ${mem}GB of memory will be allocated to the programs"

echo ""
echo "---------------"

# Initialize
mkdir -p "${dir2}/"
[ -f "${dir2}/Fastqs" ] && rm "${dir2}/Fastqs"
[ -d "${dir2}/casTLE/${logs}" ] && rm -rf "${dir2}/casTLE/${logs}"
[ -d "${dir2}/MAGeCK/${logs}" ] && rm -rf "${dir2}/MAGeCK/${logs}"

if [ -z "${tmp}" ]
then
	tmp="${dir2}/tmp"
fi



######################
## Processing samples
######################

find -L "${dir}" -maxdepth 1 -name "*${fileext}" | sed 's#.*/##' | sort -n > "${dir2}/Fastqs"

if [ ${usecastle} -eq "1" ]  ## Start of casTLE
then



###################
##    casTLE     ##
###################

	echo ""
	echo "-- casTLE --"

	mkdir -p "${dir2}/casTLE/${logs}"



	if [ -n "${oligofile}" ]
	then

		######################
		## Make bowtie index

		job="bowtie_casTLE"

		# Variables
		samplename="index"

		{
			# General SLURM parameters
			echo '#!/bin/bash'
			echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/casTLE/${logs}/${samplename}_${job}.out --error=${dir2}/casTLE/${logs}/${samplename}_${job}.err --open-mode=append"
			echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi) $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi)"
			echo "#SBATCH --requeue"
			echo "#SBATCH --time=10:00"
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
			echo "cd ${castlepath}"
			echo "${python} ./Scripts/makeIndices.py $(if [ -n "${bowtie}" ]; then echo "-b ${bowtie}-build"; fi) -t -o ${oligofile} ${screentype} ${outputbowtieindex} || if [ -f ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

			# Cleaning commands
			if [ "${debug}" != "1" ]
			then
				# Remove .sbatch
				echo "rm ${dir2}/casTLE/${samplename}_${job}.sbatch"
			fi
			echo "exit 0"

		} > "${dir2}/casTLE/${samplename}_${job}.sbatch"

		# Queue job
		SBbowdex=$(until sbatch "${dir2}/casTLE/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
		echo ""
		echo -e "\t-- makeIndices job queued --"

	fi



######################
## Generate count files

	job="counts_casTLE"

	echo ""
	echo -e "\t-- Queuing makeCounts jobs"
	echo ""

	while read -r line;
	do
		# Variables
		samplename=$(echo "${line}" | awk -F${fileext} '{print $1}')

		{
			# General SLURM parameters
			echo '#!/bin/bash'
			echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/casTLE/${logs}/${samplename}_${job}.out --error=${dir2}/casTLE/${logs}/${samplename}_${job}.err --open-mode=append"
			echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
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
			if [ -n "${SBbowdex}" ]
			then
				echo "#SBATCH --dependency=afterok:${SBbowdex##* }"
			fi

			# General commands
			echo "mkdir -p ${tmp}"
			if [ -n "${customcmd}" ]
			then
				echo "${customcmd}"
			fi

			# Job specific commands
			echo "cd ${castlepath}"
			echo "${python} ./Scripts/makeCounts.py $(if [ -n "${bowtie}" ]; then echo "-b ${bowtie}"; fi) -p ${threads} ${dir}/${line} ${dir2}/casTLE/${samplename} ${screentype} || if [ -f ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

			# Cleaning commands
			if [ "${debug}" != "1" ]
			then
				# Remove map files, not used after that step.
				echo "rm ${dir2}/casTLE/${samplename}.map.gz ${dir2}/casTLE/${samplename}.unmapped.gz"
				# Remove .sbatch
				echo "rm ${dir2}/casTLE/${samplename}_${job}.sbatch"
			fi
			echo "exit 0"

		} > "${dir2}/casTLE/${samplename}_${job}.sbatch"

		# Queue job
		SBcounts=$(until sbatch "${dir2}/casTLE/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
		if [ -n "${SBcounts##* }" ]
		then
			SBcountsIDs=${SBcountsIDs}:${SBcounts##* }
		fi
		echo -e "\t ${samplename} makeCounts job queued"

	done < "${dir2}/Fastqs"



######################
## Plot distribution of elements from count file

	job="plot_casTLE"

	# Variables
	samplename="dist"
	countfiles=$(find "${dir2}"/casTLE -name '*_counts_casTLE.sbatch' | sed ':a;N;$!ba;s/\n/ /g' | sed "s/_counts_casTLE.sbatch/_counts.csv/g")
	distlegends=$(sed ':a;N;$!ba;s/\n/ /g' "${dir2}"/Fastqs | sed "s/${fileext}//g")

	{
		# General SLURM parameters
		echo '#!/bin/bash'
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/casTLE/${logs}/${samplename}_${job}.out --error=${dir2}/casTLE/${logs}/${samplename}_${job}.err --open-mode=append"
		echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi)"
		echo "#SBATCH --requeue"
		echo "#SBATCH --time=15:00"
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
		if [ -n "${SBcountsIDs}" ]
		then
			echo "#SBATCH --dependency=afterok${SBcountsIDs}"
		fi

		# General commands
		echo "mkdir -p ${tmp}"
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}"
		fi

		# Job specific commands
		echo "cd ${castlepath}"
		echo "${python} ./Scripts/plotDist.py -of ${dir2}/casTLE/DiversityPlot ${countfiles} -l ${distlegends} $(if [ -n "${graphformat}" ]; then echo "-f ${graphformat}"; fi) || if [ -f ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

		# Cleaning commands
		if [ "${debug}" != "1" ]
		then
			# Remove .sbatch
			echo "rm ${dir2}/casTLE/${samplename}_${job}.sbatch"
		fi
		echo "exit 0"

	} > "${dir2}/casTLE/${samplename}_${job}.sbatch"

	# Queue job
	SBplotdist=$(until sbatch "${dir2}/casTLE/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
	echo ""
	echo -e "\t-- plotDist job queued --"



	######################
	## Compare count files with casTLE

	echo ""
	echo -e "\t-- Queuing analysis$(if [ -n "${permres}" ] && [ ${permres} -eq "1" ]; then echo " and permutations"; fi) jobs"
	echo ""

	# Variables
	IFS=" " read -r -a testsarray <<< "${day0} ${testgroups}"
	condA="0"
	condB="0"

	while [ ${condA} -lt ${#testsarray[@]} ]
	do
		condB=$((condA + 1))
		IFS="," read -r -a condAarray <<< "${testsarray[${condA}]}"

		while [ ${condB} -lt ${#testsarray[@]} ]
		do
			IFS="," read -r -a condBarray <<< "${testsarray[${condB}]}"

			for replicates in 0 1
			do
				# Variables
				job="analyze_casTLE"
				samplename=$(if [ -z "${condAarray[${replicates}]}" ]; then echo "${condAarray[0]}"; else echo "${condAarray[${replicates}]}";fi)_vs_$(if [ -z "${condBarray[${replicates}]}" ]; then echo "${condBarray[0]}"; else echo "${condBarray[${replicates}]}";fi)

				{
					# General SLURM parameters
					echo '#!/bin/bash'
					echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/casTLE/${logs}/${samplename}_${job}.out --error=${dir2}/casTLE/${logs}/${samplename}_${job}.err --open-mode=append"
					echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
					echo "#SBATCH --requeue"
					echo "#SBATCH --time=30:00"
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
					if [ -n "${SBcountsIDs}" ]
					then
						echo "#SBATCH --dependency=afterok${SBcountsIDs}"
					fi

					# General commands
					echo "mkdir -p ${tmp}"
					if [ -n "${customcmd}" ]
					then
						echo "${customcmd}"
					fi

					# Job specific commands
					echo "cd ${castlepath}"
					echo "${python} ./Scripts/analyzeCounts.py $(if [ -n "${mouse}" ] && [ "${mouse}" -eq "1" ] ; then echo "-m"; fi) -p ${threads} ${dir2}/casTLE/$(if [ -z "${condAarray[${replicates}]}" ]; then echo "${condAarray[0]}"; else echo "${condAarray[${replicates}]}";fi)_counts.csv ${dir2}/casTLE/$(if [ -z "${condBarray[${replicates}]}" ]; then echo "${condBarray[0]}"; else echo "${condBarray[${replicates}]}";fi)_counts.csv ${dir2}/casTLE/${samplename}_results || if [ -f ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

					# Cleaning commands
					if [ "${debug}" != "1" ]
					then
						# Remove .sbatch
						echo "rm ${dir2}/casTLE/${samplename}_${job}.sbatch"
					fi
					echo "exit 0"

				} > "${dir2}/casTLE/${samplename}_${job}.sbatch"

				# Queue job
				SBanalyze=$(until sbatch "${dir2}/casTLE/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
				if [ -n "${SBanalyze##* }" ]
				then
					SBanalyzeIDs=${SBanalyzeIDs}:${SBanalyze##* }
				fi
				echo -e "\t ${samplename} analyzeCounts job queued"

				if [ -n "${permres}" ] && [ ${permres} -eq "1" ]
				then



					######################
					## Calculate p-values for casTLE result file

					job="permut_casTLE"

					{
						# General SLURM parameters
						echo '#!/bin/bash'
						echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/casTLE/${logs}/${samplename}_${job}.out --error=${dir2}/casTLE/${logs}/${samplename}_${job}.err --open-mode=append"
						echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
						echo "#SBATCH --requeue"
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

						# Require previous job successful completion
						if [ -n "${SBanalyze}" ]
						then
							echo "#SBATCH --dependency=afterok:${SBanalyze##* }"
						fi

						# General commands
						echo "mkdir -p ${tmp}"
						if [ -n "${customcmd}" ]
						then
							echo "${customcmd}"
						fi

						# Job specific commands
						echo "cd ${castlepath}"
						echo "${python} ./Scripts/addPermutations.py $(if [ -n "${mouse}" ] && [ "${mouse}" -eq "1" ] ; then echo "-m"; fi) -p ${threads} ${dir2}/casTLE/${samplename}_results.csv ${permutations} || if [ -f ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

						# Cleaning commands
						if [ "${debug}" != "1" ]
						then
							# Remove .sbatch
							echo "rm ${dir2}/casTLE/${samplename}_${job}.sbatch"
						fi
						echo "exit 0"

					} > "${dir2}/casTLE/${samplename}_${job}.sbatch"

					# Queue job
					SBpermres=$(until sbatch "${dir2}/casTLE/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
					if [ -n "${SBpermres##* }" ]
					then
						SBpermresIDs=${SBpermresIDs}:${SBpermres##* }
					fi
					echo -e "\t ${samplename} addPermutations job queued"
					echo ""
				fi

			done

			# Move to the next comparison
			condB=$((condB + 1))
		done

		# Move to the next comparison
		condA=$((condA + 1))
	done



######################
## Combine multiple casTLE result files

	echo ""
	echo -e "\t-- Queuing combo analysis, permutation and graph jobs"
	echo ""

	# Variables
	IFS=" " read -r -a testsarray <<< "${day0} ${testgroups}"
	condA="0"
	condB="0"

	while [ ${condA} -lt ${#testsarray[@]} ]
	do
		condB=$((condA + 1))
		IFS="," read -r -a condAarray <<< "${testsarray[${condA}]}"

		while [ ${condB} -lt ${#testsarray[@]} ]
		do
			# Variables
			job="analyze_casTLE"
			IFS="," read -r -a condBarray <<< "${testsarray[${condB}]}"
			samplename="combo_${condAarray[0]}_vs_${condBarray[0]}_VS_$(if [ -z "${condAarray[1]}" ]; then echo "${condAarray[0]}"; else echo "${condAarray[1]}";fi)_vs_$(if [ -z "${condBarray[1]}" ]; then echo "${condBarray[0]}"; else echo "${condBarray[1]}";fi)"

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/casTLE/${logs}/${samplename}_${job}.out --error=${dir2}/casTLE/${logs}/${samplename}_${job}.err --open-mode=append"
				echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
				echo "#SBATCH --requeue"
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

				# Require previous job successful completion
				if [ -n "${SBanalyzeIDs}" ]
				then
					echo "#SBATCH --dependency=afterok${SBanalyzeIDs}"
				fi

				# General commands
				echo "mkdir -p ${tmp}"
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				echo "cd ${castlepath}"
				echo "${python} ./Scripts/analyzeCombo.py $(if [ -n "${mouse}" ] && [ "${mouse}" -eq "1" ] ; then echo "-m"; fi) -p ${threads} ${dir2}/casTLE/${condAarray[0]}_vs_${condBarray[0]}_results.csv ${dir2}/casTLE/$(if [ -z "${condAarray[1]}" ]; then echo "${condAarray[0]}"; else echo "${condAarray[1]}";fi)_vs_$(if [ -z "${condBarray[1]}" ]; then echo "${condBarray[0]}"; else echo "${condBarray[1]}";fi)_results.csv ${dir2}/casTLE/${samplename}  || if [ -f ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/casTLE/${samplename}_${job}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/casTLE/${samplename}_${job}.sbatch"

			# Queue job
			SBcombo=$(until sbatch "${dir2}/casTLE/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			if [ -n "${SBcombo##* }" ]
			then
				SBcomboIDs=${SBcomboIDs}:${SBcombo##* }
			fi
			echo -e "\t ${samplename} analyzeCombo job queued"



			######################
			## Calculate p-values for combination of multiple casTLE results

			# Variables
			job="permut_casTLE"

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/casTLE/${logs}/${samplename}_${job}.out --error=${dir2}/casTLE/${logs}/${samplename}_${job}.err --open-mode=append"
				echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
				echo "#SBATCH --requeue"
				echo "#SBATCH --time=6:00:00"
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
				if [ -n "${SBcombo}" ]
				then
					echo "#SBATCH --dependency=afterok:${SBcombo##* }"
				fi

				# General commands
				echo "mkdir -p ${tmp}"
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				echo "cd ${castlepath}"
				echo "${python} ./Scripts/addCombo.py $(if [ -n "${mouse}" ] && [ "${mouse}" -eq "1" ] ; then echo "-m"; fi) -p ${threads} ${dir2}/casTLE/${samplename}.csv ${permutations} || if [ -f ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/casTLE/${samplename}_${job}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/casTLE/${samplename}_${job}.sbatch"

			# Queue job
			SBpermcombo=$(until sbatch "${dir2}/casTLE/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			if [ -n "${SBpermcombo##* }" ]
			then
				SBpermcomboIDs=${SBpermcomboIDs}:${SBpermcombo##* }
			fi
			echo -e "\t ${samplename} addCombo job queued"



			######################
			## Plot casTLE result file

			# Variables
			job="volcano_plot_casTLE"

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/casTLE/${logs}/${samplename}_${job}.out --error=${dir2}/casTLE/${logs}/${samplename}_${job}.err --open-mode=append"
				echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi)"
				echo "#SBATCH --requeue"
				echo "#SBATCH --time=15:00"
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
				if [ -n "${SBcombo}" ]
				then
					echo "#SBATCH --dependency=afterok:${SBcombo##* }"
				fi

				# General commands
				echo "mkdir -p ${tmp}"
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				echo "cd ${castlepath}"
				echo "${python} ./Scripts/plotVolcano.py $(if [ -n "${mouse}" ] && [ "${mouse}" -eq "1" ] ; then echo "-m"; fi) ${dir2}/casTLE/${samplename}.csv -n \$(csvcut -c Symbol,'Combo casTLE Score','# elements 1','# elements 2' ${dir2}/casTLE/${samplename}.csv | csvgrep -c '# elements 1','# elements 2' -r ^'[6-9]|10' | csvsort -c 'Combo casTLE Score' | csvcut -c Symbol | tail -n 10 | sed ':a;N;$!ba;s/\n/ /g') $(if [ -n "${graphformat}" ]; then echo "-f ${graphformat}"; fi) || if [ -f ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/casTLE/${samplename}_${job}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/casTLE/${samplename}_${job}.sbatch"

			# Queue job
			SBplotvolcano=$(until sbatch "${dir2}/casTLE/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			if [ -n "${SBplotvolcano##* }" ]
			then
				SBplotvolcanoIDs=${SBplotvolcanoIDs}:${SBplotvolcano##* }
			fi
			echo -e "\t ${samplename} plotVolcano job queued"



			######################
			## Plot individual gene results from casTLE result file

			# Variables
			job="genes_plot_casTLE"

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/casTLE/${logs}/${samplename}_${job}.out --error=${dir2}/casTLE/${logs}/${samplename}_${job}.err --open-mode=append"
				echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi)"
				echo "#SBATCH --requeue"
				echo "#SBATCH --time=15:00"
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
				if [ -n "${SBcombo}" ]
				then
					echo "#SBATCH --dependency=afterok:${SBcombo##* }"
				fi

				# General commands
				echo "mkdir -p ${tmp}"
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				echo "cd ${castlepath}"
				echo "${python} ./Scripts/plotGenes.py $(if [ -n "${mouse}" ] && [ "${mouse}" -eq "1" ] ; then echo "-m"; fi) ${dir2}/casTLE/${condAarray[0]}_vs_${condBarray[0]}_results.csv \$(csvcut -c Symbol,'Combo casTLE Score','# elements 1','# elements 2' ${dir2}/casTLE/${samplename}.csv | csvgrep -c '# elements 1','# elements 2' -r ^'[6-9]|10' | csvsort -c 'Combo casTLE Score' | csvcut -c Symbol | tail -n 10 | sed ':a;N;$!ba;s/\n/ /g') $(if [ -n "${graphformat}" ]; then echo "-f ${graphformat}"; fi) || if [ -f ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"
				echo "${python} ./Scripts/plotGenes.py $(if [ -n "${mouse}" ] && [ "${mouse}" -eq "1" ] ; then echo "-m"; fi) ${dir2}/casTLE/$(if [ -z "${condAarray[1]}" ]; then echo "${condAarray[0]}"; else echo "${condAarray[1]}";fi)_vs_$(if [ -z "${condBarray[1]}" ]; then echo "${condBarray[0]}"; else echo "${condBarray[1]}";fi)_results.csv \$(csvcut -c Symbol,'Combo casTLE Score','# elements 1','# elements 2' ${dir2}/casTLE/${samplename}.csv | csvgrep -c '# elements 1','# elements 2' -r ^'[6-9]|10' | csvsort -c 'Combo casTLE Score' | csvcut -c Symbol | tail -n 10 | sed ':a;N;$!ba;s/\n/ /g') $(if [ -n "${graphformat}" ]; then echo "-f ${graphformat}"; fi) || if [ -f ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/casTLE/${samplename}_${job}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/casTLE/${samplename}_${job}.sbatch"

			# Queue job
			SBplotgenes=$(until sbatch "${dir2}/casTLE/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			if [ -n "${SBplotgenes##* }" ]
			then
				SBplotgenesIDs=${SBplotgenesIDs}:${SBplotgenes##* }
			fi
			echo -e "\t ${samplename} plotGenes job queued"



			######################
			## Compare enrichment of individual elements between multiple result files

			# Variables
			job="elements_plot_casTLE"

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/casTLE/${logs}/${samplename}_${job}.out --error=${dir2}/casTLE/${logs}/${samplename}_${job}.err --open-mode=append"
				echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi) $(if [ -n "${threads}" ] && [ "${threads}" -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi)"
				echo "#SBATCH --requeue"
				echo "#SBATCH --time=15:00"
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
				if [ -n "${SBcombo}" ]
				then
					echo "#SBATCH --dependency=afterok:${SBcombo##* }"
				fi

				# General commands
				echo "mkdir -p ${tmp}"
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				echo "cd ${castlepath}"
				echo "${python} ./Scripts/plotElements.py ${dir2}/casTLE/${condAarray[0]}_vs_${condBarray[0]}_results.csv ${dir2}/casTLE/$(if [ -z "${condAarray[1]}" ]; then echo "${condAarray[0]}"; else echo "${condAarray[1]}";fi)_vs_$(if [ -z "${condBarray[1]}" ]; then echo "${condBarray[0]}"; else echo "${condBarray[1]}";fi)_results.csv ${dir2}/casTLE/${samplename}_Elements -x ${condAarray[0]}_vs_${condBarray[0]} -y $(if [ -z "${condAarray[1]}" ]; then echo "${condAarray[0]}"; else echo "${condAarray[1]}";fi)_vs_$(if [ -z "${condBarray[1]}" ]; then echo "${condBarray[0]}"; else echo "${condBarray[1]}";fi) $(if [ -n "${graphformat}" ]; then echo "-f ${graphformat}"; fi) || if [ -f ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/casTLE/${samplename}_${job}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/casTLE/${samplename}_${job}.sbatch"

			# Queue job
			SBplotelem=$(until sbatch "${dir2}/casTLE/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			if [ -n "${SBplotelem##* }" ]
			then
				SBplotelemIDs=${SBplotelemIDs}:${SBplotelem##* }
			fi
			echo -e "\t ${samplename} plotElements job queued"



			######################
			## Compare effect size and confidence between multiple result files

			# Variables
			job="rep_plot_casTLE"

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/casTLE/${logs}/${samplename}_${job}.out --error=${dir2}/casTLE/${logs}/${samplename}_${job}.err --open-mode=append"
				echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi)"
				echo "#SBATCH --requeue"
				echo "#SBATCH --time=15:00"
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
				if [ -n "${SBcombo}" ]
				then
					echo "#SBATCH --dependency=afterok:${SBcombo##* }"
				fi

				# General commands
				echo "mkdir -p ${tmp}"
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				echo "cd ${castlepath}"
				echo "${python} ./Scripts/plotRep.py $(if [ -n "${mouse}" ] && [ "${mouse}" -eq "1" ] ; then echo "-m"; fi) ${dir2}/casTLE/${condAarray[0]}_vs_${condBarray[0]}_results.csv ${dir2}/casTLE/$(if [ -z "${condAarray[1]}" ]; then echo "${condAarray[0]}"; else echo "${condAarray[1]}";fi)_vs_$(if [ -z "${condBarray[1]}" ]; then echo "${condBarray[0]}"; else echo "${condBarray[1]}";fi)_results.csv ${dir2}/casTLE/${samplename}_Replicates -n \$(csvcut -c Symbol,'Combo casTLE Score','# elements 1','# elements 2' ${dir2}/casTLE/${samplename}.csv | csvgrep -c '# elements 1','# elements 2' -r ^'[6-9]|10' | csvsort -c 'Combo casTLE Score' | csvcut -c Symbol | tail -n 10 | sed ':a;N;$!ba;s/\n/ /g') -x ${condAarray[0]}_vs_${condBarray[0]} -y $(if [ -z "${condAarray[1]}" ]; then echo "${condAarray[0]}"; else echo "${condAarray[1]}";fi)_vs_$(if [ -z "${condBarray[1]}" ]; then echo "${condBarray[0]}"; else echo "${condBarray[1]}";fi) $(if [ -n "${graphformat}" ]; then echo "-f ${graphformat}"; fi) || if [ -f ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/casTLE/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/casTLE/${samplename}_${job}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/casTLE/${samplename}_${job}.sbatch"

			# Queue job
			SBplotrep=$(until sbatch "${dir2}/casTLE/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			if [ -n "${SBplotrep##* }" ]
			then
				SBplotrepIDs=${SBplotrepIDs}:${SBplotrep##* }
			fi
			echo -e "\t ${samplename} plotRep job queued"
			echo ""

			# Move to the next comparison
			condB=$((condB + 1))
		done

		# Move to the next comparison
		condA=$((condA + 1))
	done

fi ## End of casTLE



if [ ${usemageck} -eq "1" ]  ## Start of MAGeCK
then

###################
##    MAGeCK     ##
###################

	echo ""
	echo "-- MAGeCK --"

	mkdir -p "${dir2}/MAGeCK/${logs}" "${dir2}/MAGeCK/results/count" "${dir2}/MAGeCK/results/test"



	######################
	## Generate count files

	# Variables
	job="MAGeCK"
	samplename="count"

	{
		# General SLURM parameters
		echo '#!/bin/bash'
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/MAGeCK/${logs}/${samplename}_${job}.out --error=${dir2}/MAGeCK/${logs}/${samplename}_${job}.err --open-mode=append"
		echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi)"
		echo "#SBATCH --requeue"
		echo "#SBATCH --time=8:00:00"
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
		echo "mageck count --list-seq ${magecksgRNAlibrary} $(if [ -n "${mageckcontrolsgrna}" ]; then echo " --control-sgrna ${mageckcontrolsgrna}"; fi)$(if [ -n "${mageckrevcomplib}" ] && [ ${mageckrevcomplib} -eq "1" ]; then echo " --reverse-complement"; fi) --output-prefix ${dir2}/MAGeCK/results/count/all --sample-label $(sed ':a;N;$!ba;s/\n/,/g' "${dir2}"/Fastqs | sed "s/${fileext}//g")$(if [ -n "${day0}" ]; then echo " --day0-label ${day0}"; fi) --fastq ${dir}/$(sed ':a;N;$!ba;s/\n/ /g' "${dir2}"/Fastqs | sed "s| | ${dir}/|g") --pdf-report || if [ -f ${dir2}/MAGeCK/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/MAGeCK/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/MAGeCK/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

		# Cleaning commands
		if [ "${debug}" != "1" ]
		then
			# Remove .sbatch
			echo "rm ${dir2}/MAGeCK/${samplename}_${job}.sbatch"
		fi
		echo "exit 0"

	} > "${dir2}/MAGeCK/${samplename}_${job}.sbatch"

	# Queue job
	SBmageckcount=$(until sbatch "${dir2}/MAGeCK/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
	echo ""
	echo -e "\t-- MAGeCK count job queued --"



	######################
	## Perform tests and pathway analysis

	# Variables
	IFS=" " read -r -a testsarray <<< "${testgroups}"

	echo ""
	echo -e "\t-- Queuing MAGeCK test$(if [ -n "${gmtfile}" ]; then echo " and pathway analysis"; fi) jobs --"
	echo ""

	# Generate pair-wise comparison matrix
	condA="0"
	condB="0"
	while [ ${condA} -lt ${#testsarray[@]} ]
	do

		condB=$((condA + 1))



		######################
		## Perform tests

		while [ ${condB} -lt ${#testsarray[@]} ]
		do

			# Variables
			job="MAGeCK_test"
			samplename="${testsarray[${condA}]}_vs_${testsarray[${condB}]}"

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/MAGeCK/${logs}/${samplename}_${job}.out --error=${dir2}/MAGeCK/${logs}/${samplename}_${job}.err --open-mode=append"
				echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi)"
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

				# Require previous job successful completion
				if [ -n "${SBmageckcount}" ]
				then
					echo "#SBATCH --dependency=afterok:${SBmageckcount##* }"
				fi

				# General commands
				echo "mkdir -p ${tmp}"
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				echo "mageck test -c ${testsarray[${condA}]} -t ${testsarray[${condB}]} --output-prefix ${dir2}/MAGeCK/results/test/${testsarray[${condA}]}_vs_${testsarray[${condB}]} --count-table ${dir2}/MAGeCK/results/count/all.count.txt$(if [ -n "${mageckcontrolsgrna}" ]; then echo " --control-sgrna ${mageckcontrolsgrna}"; fi) --pdf-report || if [ -f ${dir2}/MAGeCK/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/MAGeCK/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/MAGeCK/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/MAGeCK/${samplename}_${job}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/MAGeCK/${samplename}_${job}.sbatch"

			# Queue job
			SBtestgroups=$(until sbatch "${dir2}/MAGeCK/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			if [ -n "${SBtestgroups##* }" ]
			then
				SBtestgroupsIDs=${SBtestgroupsIDs}:${SBtestgroups##* }
			fi
			echo -e "\t ${samplename} test job queued"

			if [ -n "${gmtfile}" ]
			then



				######################
				## Perform pathway analysis

				mkdir -p "${dir2}/MAGeCK/results/pathway"

				# Variables
				job="MAGeCK_pathway"

				{
					# General SLURM parameters
					echo '#!/bin/bash'
					echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/MAGeCK/${logs}/${samplename}_${job}.out --error=${dir2}/MAGeCK/${logs}/${samplename}_${job}.err --open-mode=append"
					echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi)"
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
					if [ -n "${SBtestgroups}" ]
					then
						echo "#SBATCH --dependency=afterok:${SBtestgroups##* }"
					fi

					# General commands
					echo "mkdir -p ${tmp}"
					if [ -n "${customcmd}" ]
					then
						echo "${customcmd}"
					fi

					# Job specific commands
					echo "mageck pathway --gene-ranking ${dir2}/MAGeCK/results/test/${testsarray[${condA}]}_vs_${testsarray[${condB}]}.gene_summary.txt --gmt-file ${gmtfile} --output-prefix ${dir2}/MAGeCK/results/pathway/${testsarray[${condA}]}_vs_${testsarray[${condB}]} || if [ -f ${dir2}/MAGeCK/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/MAGeCK/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/MAGeCK/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi"

					# Cleaning commands
					if [ "${debug}" != "1" ]
					then
						# Remove .sbatch
						echo "rm ${dir2}/MAGeCK/${samplename}_${job}.sbatch"
					fi
					echo "exit 0"

				} > "${dir2}/MAGeCK/${samplename}_${job}.sbatch"

				# Queue job
				SBmageckpathway=$(until sbatch "${dir2}/MAGeCK/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
				if [ -n "${SBmageckpathway##* }" ]
				then
					SBmageckpathwayIDs=${SBmageckpathwayIDs}:${SBmageckpathway##* }
				fi
				echo -e "\t ${samplename} pathway job queued"

			fi

			condB=$((condB + 1))
		done

		condA=$((condA + 1))
	done



	######################
	## Perform mle analysis

	# Variables
	job="MAGeCK"
	samplename="mle"

	{
		# General SLURM parameters
		echo '#!/bin/bash'
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/MAGeCK/${logs}/${samplename}_${job}.out --error=${dir2}/MAGeCK/${logs}/${samplename}_${job}.err --open-mode=append"
		echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}"
		echo "#SBATCH --requeue"
		echo "#SBATCH --time=6:00:00"
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
		if [ -n "${SBmageckcount}" ]
		then
			echo "#SBATCH --dependency=afterok:${SBmageckcount##* }"
		fi

		# General commands
		echo "mkdir -p ${tmp}"
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}"
		fi

		# Job specific commands
		# mageck mle tends to freeze randomly a lot. So here we launch it as a background process and then monitor the log file every minute.
		# Job will be automatically rescheduled if no changes have been made to the log file for the last hour, as it probably just hang up for no good reason again...
		# Update - It seems leaving --genes-varmodeling undefined and using the default value (which is either 0, 1000, 2000 or all depending where you look in the code) is the cause of the crash.
		# --genes-varmodeling is now manually set to 0, as done by Snakemake in mageck_vispr. Still leaving the log check in place, tuned up to 1h for now just in case, may need to be adjusted on other clusters.
		echo "mageck mle --genes-varmodeling 0 --count-table ${dir2}/MAGeCK/results/count/all.count.txt $(if [ -n "${matrixfile}" ]; then echo "-d ${matrixfile}"; else echo "--day0-label ${day0}"; fi) --output-prefix ${dir2}/MAGeCK/results/test/mle$(if [ -n "${mageckcontrolsgrna}" ]; then echo " --control-sgrna ${mageckcontrolsgrna}"; fi) --threads ${threads} || if [ -f ${dir2}/MAGeCK/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/MAGeCK/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/MAGeCK/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi &"
		echo 'while ps -p ${!} >/dev/null'
		echo 'do'
		echo "	checkmlealive=\`find ${dir2}/MAGeCK/${logs}/${samplename}_${job}.err -cmin -60\`"
		echo '	sleep 5s'
		echo '	if [ -z ${checkmlealive} ]'
		echo '	then'
		echo '		scontrol requeue ${SLURM_JOBID}'
		echo '		sleep 42m'
		echo '	fi'
		echo '	sleep 1m'
		echo 'done'
		echo 'sleep 30s'

		# Cleaning commands
		if [ "${debug}" != "1" ]
		then
			# Remove .sbatch
			echo "rm ${dir2}/MAGeCK/${samplename}_${job}.sbatch"
		fi
		echo "exit 0"

	} > "${dir2}/MAGeCK/${samplename}_${job}.sbatch"

	# Queue job
	SBmageckmle=$(until sbatch "${dir2}/MAGeCK/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
	echo ""
	echo -e "\t-- MAGeCK mle job queued --"



	######################
	## MAGeCKFlute

	# Variables
	job="MAGeCKFlute"
	if [ -n "${matrixfile}" ]
	then
		IFS=" " read -r -a testsarray <<< $(head -1 "${matrixfile}")
	else
		IFS=" " read -r -a testsarray <<< "${testgroups}"
	fi

	echo ""
	echo -e "\t-- Queuing MAGeCKFlute jobs"
	echo ""

	# Generate pair-wise comparison matrix
	if [ -n "${matrixfile}" ]
	then
		condA="2"	# We exclude the first 2 columns, first being the Samples and second being the Baseline. Name of the first 2 columns do not matter.
		condB="2"
	else
		condA="0"
		condB="0"
	fi

	while [ ${condA} -lt ${#testsarray[@]} ]
	do

		condB=$((condA + 1))

		while [ ${condB} -lt ${#testsarray[@]} ]
		do
			# Variables
			samplename="${testsarray[${condA}]}_vs_${testsarray[${condB}]}"

			{
				# General SLURM parameters
				echo '#!/bin/bash'
				echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/MAGeCK/${logs}/${samplename}_${job}.out --error=${dir2}/MAGeCK/${logs}/${samplename}_${job}.err --open-mode=append"
				echo "#SBATCH $(if [ -n "${mem}" ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi)"
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
				if [ -n "${SBmageckmle}" ]
				then
					echo "#SBATCH --dependency=afterok:${SBmageckmle##* }"
				fi

				# General commands
				echo "mkdir -p ${tmp}"
				if [ -n "${customcmd}" ]
				then
					echo "${customcmd}"
				fi

				# Job specific commands
				echo "cd ${tmp}"

				{
					echo '#!/usr/bin/env Rscript'
					echo ''
					echo 'library(MAGeCKFlute)	# Load MAGeCKFlute, see https://bioconductor.org/packages/release/bioc/html/MAGeCKFlute.html for install instructions'
					echo 'library(ggplot2)	# Load ggplot2, no longer automatically loaded by MAGeCKFlute'
					echo "MLE_Data <- read.table(\"${dir2}/MAGeCK/results/test/mle.gene_summary.txt\", header=TRUE)	# Load the mle analysis"
					echo "countsummary <- read.table(\"${dir2}/MAGeCK/results/count/all.countsummary.txt\", header=TRUE)	# Load the count summary table"
					echo "ctrlname=c(\"$(echo "${testsarray[${condA}]}" | sed 's/,/\",\"/g' | sed 's/-/./g')\")	# Input control sample names, MAGeCKFlute change filenames containing - to a ., so fixing the filename accordingly"
					echo "treatname=c(\"$(echo "${testsarray[${condB}]}" | sed 's/,/\",\"/g' | sed 's/-/./g')\")	# Same for treated sample names"
					echo "try(FluteMLE(MLE_Data, ctrlname, treatname, prefix=\"${samplename}\", organism=\"hsa\"))	# Pied Piper time (MAGeCK v1.4.X)"
					echo "try(FluteMLE(MLE_Data, ctrlname, treatname, proj=\"${samplename}\", organism=\"hsa\"))	# Pied Piper time (MAGeCK v1.6.X)"
				}  > "${dir2}"/MAGeCK/${job}-"${samplename}".R
				chmod +x "${dir2}"/MAGeCK/${job}-"${samplename}".R

				echo "Rscript ${dir2}/MAGeCK/${job}-${samplename}.R"
				echo "cp -R ${samplename}_Flute_Results ${dir2}/MAGeCK/results/"
				echo "cp -R MAGeCKFlute_${samplename} ${dir2}/MAGeCK/results/"
				# Move .R script inside the Flute_Results folder
				echo "mv ${dir2}/MAGeCK/${job}-${samplename}.R ${dir2}/MAGeCK/results/${samplename}_Flute_Results/"
				echo "mv ${dir2}/MAGeCK/${job}-${samplename}.R ${dir2}/MAGeCK/results/MAGeCKFlute_${samplename}/"

				# Cleaning commands
				if [ "${debug}" != "1" ]
				then
					# Remove .sbatch
					echo "rm ${dir2}/MAGeCK/${samplename}_${job}.sbatch"
				fi
				echo "exit 0"

			} > "${dir2}/MAGeCK/${samplename}_${job}.sbatch"

			# Queue job
			SBmageckflute=$(until sbatch "${dir2}/MAGeCK/${samplename}_${job}.sbatch"; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done)
			if [ -n "${SBmageckflute##* }" ]
			then
				SBmageckfluteIDs=${SBmageckfluteIDs}:${SBmageckflute##* }
			fi
			echo -e "\t ${samplename} MAGeCKFlute job queued"

			condB=$((condB + 1))
		done

		condA=$((condA + 1))
	done

fi ## End of MAGeCK



######################
## Final processing of result files

# Variables
job="notif"
samplename=$(basename "${dir2}")

{
	# General SLURM parameters
	echo '#!/bin/bash'
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${samplename}_${job}.out"
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
	echo "#SBATCH --dependency=afterok$(if [ -n "${SBanalyzeIDs}" ]; then echo "${SBanalyzeIDs}"; fi)$(if [ -n "${SBplotdist##* }" ]; then echo ":${SBplotdist##* }"; fi)$(if [ -n "${SBplotvolcanoIDs}" ]; then echo "${SBplotvolcanoIDs}"; fi)$(if [ -n "${SBplotgenesIDs}" ]; then echo "${SBplotgenesIDs}"; fi)$(if [ -n "${SBplotelemIDs}" ]; then echo "${SBplotelemIDs}"; fi)$(if [ -n "${SBplotrepIDs}" ]; then echo "${SBplotrepIDs}"; fi)$(if [ -n "${SBpermresIDs}" ]; then echo "${SBpermresIDs}"; fi)$(if [ -n "${SBcomboIDs}" ]; then echo "${SBcomboIDs}"; fi)$(if [ -n "${SBmageckcount##* }" ]; then echo ":${SBmageckcount##* }"; fi)$(if [ -n "${SBtestgroupsIDs}" ]; then echo "${SBtestgroupsIDs}"; fi)$(if [ -n "${SBmageckpathwayIDs}" ]; then echo "${SBmageckpathwayIDs}"; fi)$(if [ -n "${SBmageckmle##* }" ]; then echo ":${SBmageckmle##* }"; fi)$(if [ -n "${SBmageckfluteIDs}" ]; then echo "${SBmageckfluteIDs}"; fi)$(if [ -n "${SBpermcomboIDs}" ]; then echo ",afterany${SBpermcomboIDs}"; fi)"

	# General commands
	echo "mkdir -p ${tmp}"
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}"
	fi

	# Job specific commands
	echo "[ -f ${dir2}/Results.tar.gz ] && rm ${dir2}/Results.tar.gz"
	echo "[ -d ${dir2}/Results_${samplename}/ ] && rm -rf ${dir2}/Results_${samplename}/"
	if [ ${usecastle} -eq "1" ]  ## casTLE Results
	then
		# Creating folders for archive
		echo "mkdir -p ${dir2}/Results_${samplename}/casTLE/"
		# Fetch the few useful reports for easy download
		echo "cp -rf ${dir2}/casTLE/combo_*.csv ${dir2}/Results_${samplename}/casTLE/"
		echo "cp -rf ${dir2}/casTLE/*_results.csv ${dir2}/Results_${samplename}/casTLE/"
		echo "cp -rf ${dir2}/casTLE/*.${graphformat} ${dir2}/Results_${samplename}/casTLE/"
	fi
	if [ ${usemageck} -eq "1" ]  ## MAGeCK Results
	then
		# Creating folders for archive
		echo "mkdir -p ${dir2}/Results_${samplename}/MAGeCK/"
		# Fetch the few useful reports for easy download
		echo "cp -rf ${dir2}/MAGeCK/results/count/*_countsummary.pdf ${dir2}/Results_${samplename}/MAGeCK/"
		echo "cp -rf ${dir2}/MAGeCK/results/test/*_summary.pdf ${dir2}/Results_${samplename}/MAGeCK/"
		echo "cp -rf ${dir2}/MAGeCK/results/pathway/*.pathway_summary.txt ${dir2}/Results_${samplename}/MAGeCK/"
		echo "cp -rf ${dir2}/MAGeCK/results/*_Flute_Results/*_summary.pdf ${dir2}/Results_${samplename}/MAGeCK/"
	fi
	echo "tar --remove-files -C ${dir2} -pczf ${dir2}/Results.tar.gz Results_${samplename}"
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
		echo "rm ${dir2}/*/*.err"
		# Remove logs folder upon successfull completion. Comment to disable.
		echo "rm -rf ${dir2}/*/logs"
		echo "rm -rf ${dir2}/MAGeCK/results/*/*.log"  ## MAGeCK created log files, can be good to keep
		echo "rm -rf ${dir2}/${samplename}_${job}.out"
		# Remove Temporary directory
		echo "rm -rf ${tmp}"
		# Remove .sbatch
		echo "rm ${dir2}/${samplename}_${job}.sbatch"
	fi
	echo "exit 0"

} > "${dir2}/${samplename}_${job}.sbatch"

# Queue job
until sbatch "${dir2}/${samplename}_${job}.sbatch" >/dev/null 2>&1; do echo "Job submission failed (exit code ${?}) - Trying again in 5s"; sleep 5; done
echo ""
echo "-- Results ready notification job queued --"
echo ""

# Clean temporary files
if [ "${debug}" != "1" ]
then
	[ -f "${dir2}/Fastqs" ] && rm "${dir2}/Fastqs"
fi

# That's all folks!
exit 0
