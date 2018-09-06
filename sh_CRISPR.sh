#!/bin/bash
#
# Usage: sh_CRISPR.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]
#
##############################################################
##                       Description                        ##
##############################################################
#
# This script will process the fastq(.gz) files generated in
# a typical CRISPR screen. It will, if needed, create a
# reference file of all the Indices using bowtie (NOT bowtie2).
# It will then use casTLE scripts to analyze the screen and
# generate basic graphs https://bitbucket.org/dmorgens/castle/
#
##############################################################
##                  Configurable variables                  ##
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
mem="64"
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
#
## bowtie options
# path to bowtie, bowtie-build needs to be in the same folder (probably is the case)
# bowtie="/usr/local/bin/bowtie""
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
# casTLE options
# casTLE folder location
# download the last version from https://bitbucket.org/dmorgens/castle/
castlepath="/home/user/scripts/dmorgens-castle/"
#
# Number of permutations to generate p-values.
# For a first pass, use 5x the permutations (so for 20,000 genes that is 100,000 permutations).
# For publication, use 50x (Default, takes time)
permutations="1000000"
#
# Perform permutations on the individual result files.
# Permutations are always calculated for the combo file.
# 0 = No ; 1 = Yes
permres="0"
#
# Output format for graphs. If left empty default value is png
# (png, pdf, eps)
graphformat="pdf"
#
# Screen performed in mouse cells?
# Default is 0, for human cells
# 0 = No ; 1 = Yes
mouse="0"
#
# Enter your sample names (not including .fastq or .fastq.gz) for comparisons.
# You should only have 4 samples, organized in 2 pairs
# analyzecounts="Untreated1,Treated1 Untreated2,Treated2"
analyzecounts=""
#
#
# IFTTT options
#
# Trigger IFTTT when script is done.
# You must register the "Maker channel" on https://ifttt.com/maker
# Copy your private key here. Leave blank to disable this function.
iftttkey="AbCd_15CdhUIvbsFJTHGMcfgjsdHRTgcyjt" # Not a real key, you have to use your own private key
#
# Event name used in your IFTTT recipes.
# The maker channel will look for the combination private key + event name to then trigger your recipe
# You can create a recipe to send an email, a text message or a push notification.
iftttevent="CRISPR"
#
#
## Setup done. You should not need to edit below this point ##

echo ""

# Get fastq directory
dir="$1"

# Get destination directory
dir2="$2"

# Get config file location
config="$3"

# Check paths and trailing / in directories
if [ -z "${dir}" -o -z "${dir2}" ]
then
	echo "Usage: sh_CRISPR.sh </path/to/fastq(.gz)/folder> </path/to/casTLE/results/folder> [/path/to/config/file.ini]"
	echo ""
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
		echo "Usage: sh_CRISPR.sh </path/to/fastq(.gz)/folder> </path/to/casTLE/results/folder> [/path/to/config/file.ini]"
		echo ""
		exit
	fi
fi

if [ -n "${oligofile}" ]
then
	if [ ${oligofile: -4} != ".csv"  ]
	then
		echo "Invalid population oligo file file detected. Is ${oligofile} a .csv file?"
		echo ""
		exit
	fi
fi

if [ -n "${analyzecounts}" ]
then
	declare -a samplearray=(`echo ${analyzecounts} | sed 's/,/ /g'`)
	if [ ${#samplearray[@]} != "4" ]
	then
		echo "This pipeline is only designed to work with 2 replicates, 4 samples total"
		echo "Current sample names:" ${samplearray[@]}
		echo ""
		exit
	fi
else
	echo "Sample names and comparisons should be entered in analyzecounts"
	echo ""
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
echo "-- Hardware"
echo ""
echo "Up to ${threads} CPUs will be used"
echo "Up to ${mem}GB of memory will be allocated to the programs"

# Initialize
[ -f ${dir2}/Fastqs ] && rm ${dir2}/Fastqs
[ -d ${dir2}/logs ] && rm -rf ${dir2}/logs
mkdir -p ${dir2}/
mkdir -p ${dir2}/${logs}

if [ -z ${tmp} ]
then
	tmp="${dir2}/tmp"
fi



######################
## Processing samples
######################



if [ -n "${oligofile}" ]
then

	###################### 
	## Make bowtie index
	job="bowtie"

	# Variables
	samplename="index"

	# General SLURM parameters
	echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH `if [ ! -z ${mem} ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi` `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --time=10:00" >> ${dir2}/${samplename}_${job}.sbatch
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
	echo "cd ${castlepath}" >> ${dir2}/${samplename}_${job}.sbatch
	echo "python Scripts/makeIndices.py `if [ -n \"${bowtie}\" ]; then echo \"-b ${bowtie}-build\"; fi` -t -o ${oligofile} ${screentype} ${outputbowtieindex} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch

	# Cleaning commands
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
	echo "exit 0" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBbowdex=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	echo ""
	echo "-- makeIndices job queued"

fi



echo ""
echo "-- Queuing makeCounts jobs --"
echo ""

###################### 
## Align sequence files
job="counts"

ls ${dir}/ | grep ${fileext} > ${dir2}/Fastqs

while read line;
do

	# Variables
	samplename=`echo ${line} | awk -F${fileext} '{print $1}'`

	# General SLURM parameters
	echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}" >> ${dir2}/${samplename}_${job}.sbatch
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
	if [ -n "${SBbowdex}" ]
	then
		echo "#SBATCH --dependency=afterok:${SBbowdex##* }" >> ${dir2}/${samplename}_${job}.sbatch
	fi

	# General commands
	echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
	fi

	# Job specific commands
	echo "cd ${castlepath}" >> ${dir2}/${samplename}_${job}.sbatch
	echo "python Scripts/makeCounts.py `if [ -n \"${bowtie}\" ]; then echo \"-b ${bowtie}\"; fi` -p ${threads} ${dir}/${line} ${dir2}/${samplename} ${screentype} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch

	# Cleaning commands
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
	echo "exit 0" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBcounts=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	if [ -n "${SBcounts##* }" ]
	then
		SBcountsIDs=${SBcountsIDs}:${SBcounts##* }
	fi
	echo -e "\t ${samplename} makeCounts job queued"

done < ${dir2}/Fastqs



###################### 
## Plot distribution of elements from count file
job="dist"

# Variables
samplename="plot"
countfiles=`find ${dir2} -name '*_counts.sbatch' | sed ':a;N;$!ba;s/\n/ /g' | sed "s/.sbatch/.csv/g"`
distlegends=`cat ${dir2}/Fastqs | sed ':a;N;$!ba;s/\n/ /g' | sed "s/${fileext}//g"`

# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH `if [ ! -z ${mem} ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi` `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=15:00" >> ${dir2}/${samplename}_${job}.sbatch
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
if [ -n "${SBcountsIDs}" ]
then
	echo "#SBATCH --dependency=afterok${SBcountsIDs}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Job specific commands
echo "cd ${castlepath}" >> ${dir2}/${samplename}_${job}.sbatch
echo "python Scripts/plotDist.py -of ${dir2}/DiversityPlot ${countfiles} -l ${distlegends} `if [ -n \"${graphformat}\" ]; then echo \"-f ${graphformat}\"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch

# Cleaning commands
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
echo "exit 0" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBplotdist=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo ""
echo "-- plotDist job queued"



echo ""
echo "-- Queuing analysis jobs --"
echo ""

# Variables
declare -a resultarray=(`echo ${analyzecounts}`)

# Generate comparison matrices
res="0"
while [ ${res} -lt ${#resultarray[@]} ]
do
	###################### 
	## Compare two count files with casTLE
	job="analyze"

	# Variables
	declare -a comparray=(`echo ${resultarray[${res}]} | sed 's/,/ /g'`)
	samplename=`echo ${comparray[0]}_vs_${comparray[1]}`

	# General SLURM parameters
	echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --time=30:00" >> ${dir2}/${samplename}_${job}.sbatch
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
	if [ -n "${SBcountsIDs}" ]
	then
		echo "#SBATCH --dependency=afterok${SBcountsIDs}" >> ${dir2}/${samplename}_${job}.sbatch
	fi

	# General commands
	echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
	fi

	# Job specific commands
	echo "cd ${castlepath}" >> ${dir2}/${samplename}_${job}.sbatch
	echo "python Scripts/analyzeCounts.py `if [ -n \"${mouse}\" ] && [ \"${mouse}\" -eq "1" ] ; then echo \"-m\"; fi` -p ${threads} ${dir2}/${comparray[0]}_counts.csv ${dir2}/${comparray[1]}_counts.csv ${dir2}/${comparray[0]}_vs_${comparray[1]}_results || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch

	# Cleaning commands
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
	echo "exit 0" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBanalyze=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	if [ -n "${SBanalyze##* }" ]
	then
		SBanalyzeIDs=${SBanalyzeIDs}:${SBanalyze##* }
	fi
	echo -e "\t ${comparray[0]}_vs_${comparray[1]} analyzeCounts job queued"



	if [ ! -z ${permres} ] && [ ${permres} -eq "1" ]
	then
		###################### 
		## Calculate p-values for casTLE result file
		job="permut"

		# General SLURM parameters
		echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH --cpus-per-task=${threads} --mem-per-cpu=6000MB" >> ${dir2}/${samplename}_${job}.sbatch
		echo "#SBATCH --time=3:00:00" >> ${dir2}/${samplename}_${job}.sbatch
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
		if [ -n "${SBanalyze}" ]
		then
			echo "#SBATCH --dependency=afterok:${SBanalyze##* }" >> ${dir2}/${samplename}_${job}.sbatch
		fi

		# General commands
		echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
		fi

		# Job specific commands
		echo "cd ${castlepath}" >> ${dir2}/${samplename}_${job}.sbatch
		echo "python Scripts/addPermutations.py `if [ -n \"${mouse}\" ] && [ \"${mouse}\" -eq "1" ] ; then echo \"-m\"; fi` -p ${threads} ${dir2}/${comparray[0]}_vs_${comparray[1]}_results.csv ${permutations} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch

		# Cleaning commands
		# remove .sbatch
		echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
		echo "exit 0" >> ${dir2}/${samplename}_${job}.sbatch

		# Queue job
		SBpermres=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
		if [ -n "${SBpermres##* }" ]
		then
			SBpermresIDs=${SBpermresIDs}:${SBpermres##* }
		fi
		echo -e "\t ${comparray[0]}_vs_${comparray[1]} addPermutations job queued"
	fi

	# Move to the next comparison
	res=$((res + 1))
done



###################### 
## Combine multiple casTLE result files
job="analyze"

# Variables
samplename="combo"
declare -a comboarray=(`echo ${analyzecounts} | sed 's/,/ /g'`)

# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=3:00:00" >> ${dir2}/${samplename}_${job}.sbatch
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
if [ -n "${SBanalyzeIDs}" ]
then
	echo "#SBATCH --dependency=afterok${SBanalyzeIDs}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Job specific commands
echo "cd ${castlepath}" >> ${dir2}/${samplename}_${job}.sbatch
echo "python Scripts/analyzeCombo.py `if [ -n \"${mouse}\" ] && [ \"${mouse}\" -eq "1" ] ; then echo \"-m\"; fi` -p ${threads} ${dir2}/${comboarray[0]}_vs_${comboarray[1]}_results.csv ${dir2}/${comboarray[2]}_vs_${comboarray[3]}_results.csv ${dir2}/${samplename}  || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch

# Cleaning commands
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
echo "exit 0" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBcombo=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo ""
echo "-- analyzeCombo job queued"



###################### 
## Calculate p-values for combination of multiple casTLE results
job="permut"

# Variables
samplename="combo"
	
# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --cpus-per-task=${threads} --mem-per-cpu=6000MB" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=3:00:00" >> ${dir2}/${samplename}_${job}.sbatch
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
if [ -n "${SBcombo}" ]
then
	echo "#SBATCH --dependency=afterok:${SBcombo##* }" >> ${dir2}/${samplename}_${job}.sbatch
fi

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi
		
# Job specific commands
echo "cd ${castlepath}" >> ${dir2}/${samplename}_${job}.sbatch
echo "python Scripts/addCombo.py `if [ -n \"${mouse}\" ] && [ \"${mouse}\" -eq "1" ] ; then echo \"-m\"; fi` -p ${threads} ${dir2}/combo.csv ${permutations} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch

# Cleaning commands
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
echo "exit 0" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBpermcombo=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo ""
echo "-- addCombo job queued"



###################### 
## Plot casTLE result file
job="volcano"

# Variables
samplename="plot"
	
# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH `if [ ! -z ${mem} ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi` `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=15:00" >> ${dir2}/${samplename}_${job}.sbatch
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
if [ -n "${SBcombo}" ]
then
	echo "#SBATCH --dependency=afterok:${SBcombo##* }" >> ${dir2}/${samplename}_${job}.sbatch
fi

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Job specific commands
echo "cd ${castlepath}" >> ${dir2}/${samplename}_${job}.sbatch
echo "python Scripts/plotVolcano.py `if [ -n \"${mouse}\" ] && [ \"${mouse}\" -eq "1" ] ; then echo \"-m\"; fi` ${dir2}/combo.csv -n \`csvcut -c Symbol,'Combo casTLE Score' ${dir2}/combo.csv | csvsort -c 'Combo casTLE Score' | csvcut -c Symbol | tail -n 10 | sed ':a;N;$!ba;s/\n/ /g'\` `if [ -n \"${graphformat}\" ]; then echo \"-f ${graphformat}\"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch

# Cleaning commands
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
echo "exit 0" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBplotvolcano=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo ""
echo "-- plotVolcano job queued"



###################### 
## Plot individual gene results from casTLE result file
job="genes"

# Variables
samplename="plot"

# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH `if [ ! -z ${mem} ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi` `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=15:00" >> ${dir2}/${samplename}_${job}.sbatch
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
if [ -n "${SBcombo}" ]
then
	echo "#SBATCH --dependency=afterok:${SBcombo##* }" >> ${dir2}/${samplename}_${job}.sbatch
fi

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi
		
# Job specific commands
echo "cd ${castlepath}" >> ${dir2}/${samplename}_${job}.sbatch
echo "python Scripts/plotGenes.py `if [ -n \"${mouse}\" ] && [ \"${mouse}\" -eq "1" ] ; then echo \"-m\"; fi` ${dir2}/${comboarray[0]}_vs_${comboarray[1]}_results.csv \`csvcut -c Symbol,'Combo casTLE Score' ${dir2}/combo.csv | csvsort -c 'Combo casTLE Score' | csvcut -c Symbol | tail -n 10 | sed ':a;N;$!ba;s/\n/ /g'\` `if [ -n \"${graphformat}\" ]; then echo \"-f ${graphformat}\"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
echo "python Scripts/plotGenes.py `if [ -n \"${mouse}\" ] && [ \"${mouse}\" -eq "1" ] ; then echo \"-m\"; fi` ${dir2}/${comboarray[2]}_vs_${comboarray[3]}_results.csv \`csvcut -c Symbol,'Combo casTLE Score' ${dir2}/combo.csv | csvsort -c 'Combo casTLE Score' | csvcut -c Symbol | tail -n 10 | sed ':a;N;$!ba;s/\n/ /g'\` `if [ -n \"${graphformat}\" ]; then echo \"-f ${graphformat}\"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch

# Cleaning commands
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
echo "exit 0" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBplotgenes=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo ""
echo "-- plotGenes job queued"



###################### 
## Compare enrichment of individual elements between multiple result files
job="elements"

# Variables
samplename="plot"

# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH `if [ ! -z ${mem} ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi` `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=15:00" >> ${dir2}/${samplename}_${job}.sbatch
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
if [ -n "${SBcombo}" ]
then
	echo "#SBATCH --dependency=afterok:${SBcombo##* }" >> ${dir2}/${samplename}_${job}.sbatch
fi

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Job specific commands
echo "cd ${castlepath}" >> ${dir2}/${samplename}_${job}.sbatch
echo "python Scripts/plotElements.py ${dir2}/${comboarray[0]}_vs_${comboarray[1]}_results.csv ${dir2}/${comboarray[2]}_vs_${comboarray[3]}_results.csv ${dir2}/Elements -x ${comboarray[0]}_vs_${comboarray[1]} -y ${comboarray[2]}_vs_${comboarray[3]} `if [ -n \"${graphformat}\" ]; then echo \"-f ${graphformat}\"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch

# Cleaning commands
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
echo "exit 0" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBplotelem=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo ""
echo "-- plotElementsjob queued"



###################### 
## Compare effect size and confidence between multiple result files
job="rep"

# Variables
samplename="plot"

# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH `if [ ! -z ${mem} ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi` `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=15:00" >> ${dir2}/${samplename}_${job}.sbatch
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
if [ -n "${SBcombo}" ]
then
	echo "#SBATCH --dependency=afterok:${SBcombo##* }" >> ${dir2}/${samplename}_${job}.sbatch
fi

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi
		
# Job specific commands
echo "cd ${castlepath}" >> ${dir2}/${samplename}_${job}.sbatch
echo "python Scripts/plotRep.py `if [ -n \"${mouse}\" ] && [ \"${mouse}\" -eq "1" ] ; then echo \"-m\"; fi` ${dir2}/${comboarray[0]}_vs_${comboarray[1]}_results.csv ${dir2}/${comboarray[2]}_vs_${comboarray[3]}_results.csv ${dir2}/Replicates -n \`csvcut -c Symbol,'Combo casTLE Score' ${dir2}/combo.csv | csvsort -c 'Combo casTLE Score' | csvcut -c Symbol | tail -n 10 | sed ':a;N;$!ba;s/\n/ /g'\` -x ${comboarray[0]}_vs_${comboarray[1]} -y ${comboarray[2]}_vs_${comboarray[3]} `if [ -n \"${graphformat}\" ]; then echo \"-f ${graphformat}\"; fi` || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch

# Cleaning commands
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
echo "exit 0" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBplotrep=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo ""
echo "-- plotRep job queued"



######################
## Final processing of result files

# Variables
job="notif"
samplename="`basename ${dir2}`-CRISPR"

# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH `if [ ! -z ${mem} ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi` `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=10:00" >> ${dir2}/${samplename}_${job}.sbatch
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
echo "#SBATCH --dependency=afterok`if [ -n \"${SBanalyzeIDs}\" ]; then echo \"${SBanalyzeIDs}\"; fi``if [ -n \"${SBplotvolcano##* }\" ]; then echo \":${SBplotvolcano##* }\"; fi``if [ -n \"${SBplotgenes##* }\" ]; then echo \":${SBplotgenes##* }\"; fi``if [ -n \"${SBplotelem##* }\" ]; then echo \":${SBplotelem##* }\"; fi``if [ -n \"${SBplotrep##* }\" ]; then echo \":${SBplotrep##* }\"; fi``if [ -n \"${SBpermresIDs}\" ]; then echo \"${SBpermresIDs}\"; fi``if [ -n \"${SBcombo##* }\" ]; then echo \":${SBcombo##* }\"; fi``if [ -n \"${SBpermcombo##* }\" ]; then echo \":${SBpermcombo##* }\"; fi`" >> ${dir2}/${samplename}_${job}.sbatch

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Job specific commands
# Creating folders for archive
echo "[ -f ${dir2}/Results.tar.gz ] && rm ${dir2}/Results.tar.gz" >> ${dir2}/${samplename}_${job}.sbatch
echo "[ -d ${dir2}/Results_`basename ${dir2}` ] && rm -rf ${dir2}/Results" >> ${dir2}/${samplename}_${job}.sbatch
echo "mkdir -p ${dir2}/Results_`basename ${dir2}`/" >> ${dir2}/${samplename}_${job}.sbatch
# List all csv, png, pdf and zip files
echo "csvarchive=\`ls ${dir2}/*.csv\`" >> ${dir2}/${samplename}_${job}.sbatch
echo "epsarchive=\`ls ${dir2}/*.eps\`" >> ${dir2}/${samplename}_${job}.sbatch
echo "pdfarchive=\`ls ${dir2}/*.pdf\`" >> ${dir2}/${samplename}_${job}.sbatch
echo "pngarchive=\`ls ${dir2}/*.png\`" >> ${dir2}/${samplename}_${job}.sbatch
echo "ziparchive=\`ls ${dir2}/*.zip\`" >> ${dir2}/${samplename}_${job}.sbatch
# Create an archive of all csv, html, pdf and zip files for easy download
echo "for i in \${csvarchive}; do cp -rf \$i ${dir2}/Results_`basename ${dir2}`/\`basename \$i\`; done" >> ${dir2}/${samplename}_${job}.sbatch
echo "for i in \${epsarchive}; do cp -rf \$i ${dir2}/Results_`basename ${dir2}`/\`basename \$i\`; done" >> ${dir2}/${samplename}_${job}.sbatch
echo "for i in \${pdfarchive}; do cp -rf \$i ${dir2}/Results_`basename ${dir2}`/\`basename \$i\`; done" >> ${dir2}/${samplename}_${job}.sbatch
echo "for i in \${pngarchive}; do cp -rf \$i ${dir2}/Results_`basename ${dir2}`/\`basename \$i\`; done" >> ${dir2}/${samplename}_${job}.sbatch
echo "for i in \${ziparchive}; do cp -rf \$i ${dir2}/Results_`basename ${dir2}`/\`basename \$i\`; done" >> ${dir2}/${samplename}_${job}.sbatch
echo "tar --remove-files -C ${dir2} -pczf ${dir2}/Results.tar.gz Results_`basename ${dir2}`" >> ${dir2}/${samplename}_${job}.sbatch
# Trigger IFTTT maker channel event when it's ready, nice isn't it?
if [ -n "${iftttkey}" ] && [ -n "${iftttevent}" ]
then
		echo "curl -X POST -H \"Content-Type: application/json\" -d '{ \"value1\" : \"`basename ${dir2}`\" , \"value2\" : \"`whoami`\"}' https://maker.ifttt.com/trigger/${iftttevent}/with/key/${iftttkey}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Cleaning commands
# Remove error files upon successful completion. Comment to disable.
echo "rm ${dir2}/*.err" >> ${dir2}/${samplename}_${job}.sbatch
# Remove logs folder upon successfull completion. Comment to disable.
#echo "rm -rf ${dir2}/logs" >> ${dir2}/${samplename}_${job}.sbatch
# Remove Temporary directory
echo "rm -rf ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch
echo "exit 0" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBnotif=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo ""
echo "-- Results ready notification job queued --"
echo ""

# Clean temporary files
[ -f ${dir2}/Fastqs ] && rm ${dir2}/Fastqs

# That's all folks!
exit 0
