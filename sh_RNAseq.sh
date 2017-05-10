#!/bin/bash
#
# Usage: sh_RNAseq.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]
#
##############################################################
##                       Description                        ##
##############################################################
#
# This script will process fastq(.gz) files and align them to
# a reference genome using either STAR (recommended), hishat2
# or tophat2.
# Differential expression will then be computed using
# cufflinks. If STAR is used then RSEM will also be used to
# generate additional files for DESeq2.
#
##############################################################
##                  Configurable variables                  ##
##############################################################
#
## Hardware options
#
# Maximum number of threads (or CPUs) to request and allocate to programs.
# In some case less than this value may automatically be allowed.
threads=$(nproc --all --ignore=1)
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
# By default only send FAIL notifications for all jobs and FAIL,END for the last one.
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
# Process undetermined/unmatched files (tag not properly recognized)
# 0 = No ; 1 = Yes
underdet="0"
#
# Process BLANK files (Sample name should be "BLANK")
# 0 = No ; 1 = Yes
blank="0"
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
## Alignment options
#
# tophat2 (bowtie2) indexed reference genome location
th_refgenome="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
#
# hisat2 indexed reference genome location
ha_refgenome="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Sequence/Hisat2Index/genome"
#
# STAR indexed reference genome location
st_refgenome="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Sequence/StarIndex/"
#
# RSEM indexed reference genome location
rs_refgenome="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Sequence/RsemIndex/genome"
#
# Annotation files folder.
# It must be compatible with the indexed bowtie2 file used for alignment
# or used while creating the hishat2 / STAR indexes.
# This folder should contain at least a gene.gtf file.
# If bowtie2 is used it should also contains the known.gff, known.fa, etc, files
# else they will be generated during the first run.
gtf="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes"
#
# Reference genome (fasta) file.
# It must be the same one that was indexed by bowtie2 for alignment
fasta_refgenome="/Tools/RefGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
#
# Read group parameters
# Library: if empty LB will be set to the destination folder name.
LB=""
#
# Platform used to produce the reads. Valid values: CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT and PACBIO.
PL="ILLUMINA"
#
#
## Cuffdiff and differential expression options
#
# Separe replicate samples by a comma and groups by a space. Sample name is case sensitive.
# example:
# groupedsamples="1d-Control-Rep1,1d-Control-Rep2,1d-Control-Rep3 5d-Control-Rep1,5d-Control-Rep2,5d-Control-Rep3 1d-Experiment-Rep1,1d-Experiment-Rep2,1d-Experiment-Rep3 5d-Experiment-Rep1,5d-Experiment-Rep2,5d-Experiment-Rep3"
groupedsamples=""
#
# Separe group labels by a comma. Groups should be named in the same order than the samples
# example:
# labels="1d-Control,5d-Control,1d-Experiment,5d-Experiment"
labels=""
#
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
iftttevent="RNAseq"
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
	echo "Usage: sh_RNAseq.sh </path/to/fastq(.gz)/folder> </path/to/Aligned(.bam)/destination/folder> [/path/to/config/file.ini]"
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
		echo "Usage: sh_RNAseq.sh </path/to/fastq(.gz)/folder> </path/to/Aligned(.bam)/destination/folder> [/path/to/config/file.ini]"
		exit
	fi
fi

if [ -z "${labels}" ] || [ -z "${groupedsamples}" ]
then
	echo "Error: both \"groupedsamples\" and \"labels\" variables need to be configured properly"
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

echo ""
echo "-- Hardware"
echo ""
echo "Up to ${threads} CPUs will be used"
echo "Up to ${mem}GB of memory will be allocated to the programs"

# Initialize
mkdir -p ${dir2}/
[ -f ${dir2}/files1 ] && rm ${dir2}/files1
[ -f ${dir2}/files2 ] && rm ${dir2}/files2
[ -f ${dir2}/Fastqs ] && rm ${dir2}/Fastqs
[ -f ${dir2}/matrices ] && rm ${dir2}/matrices
[ -f ${dir2}/assembly_GTF_list.txt ] && rm ${dir2}/assembly_GTF_list.txt
[ -d ${dir2}/logs ] && rm -rf ${dir2}/logs
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
		ls ${dir}/ --hide=*ndetermined* --hide=*nmatched* --hide=BLANK* | grep _R1 > ${dir2}/files1
		ls ${dir}/ --hide=*ndetermined* --hide=*nmatched* --hide=BLANK* | grep _R2 > ${dir2}/files2
		paste ${dir2}/files1 ${dir2}/files2 > ${dir2}/Fastqs
	else
		# Remove all Undetermined_* files
		ls ${dir}/ --hide=*ndetermined* --hide=*nmatched* | grep _R1 > ${dir2}/files1
		ls ${dir}/ --hide=*ndetermined* --hide=*nmatched* | grep _R2 > ${dir2}/files2
		paste ${dir2}/files1 ${dir2}/files2 > ${dir2}/Fastqs
	fi
else
	if [ ${blank} -eq "0" ]
	then
		# Remove all BLANK* files
		ls ${dir}/ --hide=BLANK* | grep _R1 > ${dir2}/files1
		ls ${dir}/ --hide=BLANK* | grep _R2 > ${dir2}/files2
		paste ${dir2}/files1 ${dir2}/files2 > ${dir2}/Fastqs
	else
		# Process all the files!
		ls ${dir}/ | grep _R1 > ${dir2}/files1
		ls ${dir}/ | grep _R2 > ${dir2}/files2
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
		echo "#SBATCH `if [ ! -z ${mem} ] && [ ${mem} -gt "8" ]; then echo "--mem=8000"; else echo "--mem=${mem}000"; fi` `if [ ! -z ${threads} ] && [ ${threads} -gt "2" ]; then echo "--cpus-per-task=2"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
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
		echo "java -Xmx`if [ ! -z ${mem} ] && [ ${mem} -gt "8" ]; then echo "8"; else echo "${mem}"; fi`g -Djava.io.tmpdir=${tmp} -jar ${Trimmomatic}/trimmomatic.jar PE -threads `if [ ! -z ${threads} ] && [ ${threads} -gt "2" ]; then echo "2"; else echo "${threads}"; fi` -phred33 ${dir}/${read1} ${dir}/${read2} ${dir2}/${trimout1} ${unpaired1} ${dir2}/${trimout2} ${unpaired2} ILLUMINACLIP:${Trimmomatic}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch

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
	## Align reads with hisat2 / tophat2 / star
	job="align"

	# variables
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
	echo "#SBATCH --dependency=afterok:${SBtrim##* }" >> ${dir2}/${samplename}_${job}.sbatch
	
	# General commands
	echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# Job specific commands
	if [ ! -z ${align} ] && [ ${align} = "star" ]
	then
		# star alignment job using ENCODE standard options
		echo "STAR --runMode alignReads --runThreadN ${threads} --limitBAMsortRAM ${mem}000000000 --outTmpDir ${tmp}/${samplename}_aRtmp --twopassMode Basic --genomeDir ${st_refgenome} --readFilesIn ${dir2}/${trimout1} ${dir2}/${trimout2} `if [ -n "${fastqgz}" ] && [ ${trim} -ne "1" ]; then echo "--readFilesCommand zcat"; fi` --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --quantMode TranscriptomeSAM GeneCounts --alignEndsType EndToEnd --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${dir2}/${samplename}. || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
		echo "mv ${dir2}/${samplename}.Aligned.sortedByCoord.out.bam ${dir2}/${bamsortedout}" >> ${dir2}/${samplename}_${job}.sbatch
#		echo "STAR --runMode inputAlignmentsFromBAM --runThreadN ${threads} --outTmpDir ${tmp}/${samplename}_iAFBtmp --inputBAMfile ${dir2}/${bamsortedout} --outWigType bedGraph --outWigStrand Stranded --outWigReferencesPrefix chr --outFileNamePrefix ${dir2}/${samplename}. || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
#		echo "STAR --runMode inputAlignmentsFromBAM --runThreadN ${threads} --outTmpDir ${tmp}/${samplename}_iAFBtmp --inputBAMfile ${dir2}/${bamsortedout} --outWigType wiggle --outWigStrand Stranded --outWigReferencesPrefix chr --outFileNamePrefix ${dir2}/${samplename}. || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
	elif [ ! -z ${align} ] && [ ${align} = "hisat2" ]
	then
		# hisat2 alignment job
		echo "hisat2 -p ${threads} --phred33 --dta-cufflinks --no-softclip --rg-id ${LB}_${SM} --rg CN:${CN} --rg LB:${LB} --rg PL:${PL} --rg PU:${PU} --rg SM:${SM} -x ${ha_refgenome} -1 ${dir2}/${trimout1} -2 ${dir2}/${trimout2} -S ${dir2}/${samout} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
	elif [ ! -z ${align} ] && [ ${align} -= "tophat2" ]
	then
		# tophat2 alignment job
		echo "tophat2 -G ${gtf}/genes.gtf --transcriptome-index ${gtf}/known -p ${threads} -o ${dir2} ${th_refgenome} ${dir2}/${trimout1} ${dir2}/${trimout2} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	if [ ! -z ${align} ] && ([ ${align} = "hisat2" ] || [ ${align} = "tophat2" ])
	then
		# Convert .sam to .bam
		echo "samtools view -bS -o ${dir2}/${bamout} ${dir2}/${samout} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
		# Sort .bam file
		echo "samtools sort -@ ${threads} -o ${dir2}/${bamsortedout} -O bam -T ${tmp} ${dir2}/${bamout} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	# Index sorted .bam
	echo "samtools index ${dir2}/${bamsortedout} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch

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
	# clean STAR intermediate files
	echo "rm -rf ${dir2}/${samplename}.Log.progress.out ${dir2}/${samplename}._STARgenome ${dir2}/${samplename}._STARpass1" >> ${dir2}/${samplename}_${job}.sbatch
	echo "mv ${dir2}/${samplename}.Log.out ${dir2}/${logs}/${samplename}.star.log" >> ${dir2}/${samplename}_${job}.sbatch
	echo "mv ${dir2}/${samplename}.Log.final.out ${dir2}/${logs}/${samplename}.star.stats" >> ${dir2}/${samplename}_${job}.sbatch
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBalign=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	echo -e "\t Alignment job queued"



	###################### 
	## Run fastqc to generate quality control files
	job="fqc"
			
	if [ ! -z ${fastqc} ] && [ ${fastqc} -eq "1" ]
	then
	
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
		echo "fastqc -o ${dir2}/ --noextract ${dir2}/${trimout1} ${dir2}/${trimout2} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
		
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
		SBfqcIDs=${SBfqcIDs}:${SBfqc##* }
		echo -e "\t FastQC job queued"
	fi



	###################### 
	## Count reads
	job="rsem"
	
	if [ ! -z ${align} ] && [ ${align} = "star" ]
	then
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
		echo "#SBATCH --dependency=afterok:${SBalign##* }" >> ${dir2}/${samplename}_${job}.sbatch
	
		# General commands
		echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
		if [ -n "${customcmd}" ]
		then
			echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
		fi
	
		# Job specific commands
		echo "rsem-calculate-expression -p ${threads} --temporary-folder ${tmp}/${samplename} --paired-end --bam --no-bam-output --calc-ci ${dir2}/${samplename}.Aligned.toTranscriptome.out.bam ${rs_refgenome} ${dir2}/${samplename}.rsem || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
		echo "mv ${dir2}/${samplename}.rsem.stat ${dir2}/${samplename}.rsem" >> ${dir2}/${samplename}_${job}.sbatch

		# Cleaning commands
		# remove .sbatch
		echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

		# Queue job
		SBrsem=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
		SBrsemIDs=${SBrsemIDs}:${SBrsem##* }
		echo -e "\t RSEM counting job queued"
	fi



	###################### 
	## Assemble transcriptomes
	job="cufflinks"
	
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
	echo "#SBATCH --dependency=afterok:${SBalign##* }" >> ${dir2}/${samplename}_${job}.sbatch
	
	# General commands
	echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# Job specific commands
	# Create a sample specific folder, and as cufflinks do not allow to specify an output directory then cd to it.
	echo "mkdir -p ${dir2}/${samplename}.cuff" >> ${dir2}/${samplename}_${job}.sbatch
	echo "cd ${dir2}/${samplename}.cuff" >> ${dir2}/${samplename}_${job}.sbatch
	# Create cufflinks job
	echo "cufflinks -u -p ${threads} -g ${gtf}/genes.gtf -b ${fasta_refgenome} ${dir2}/${bamsortedout} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
	# Generate a list of sample specific transcripts.gtf files to be merged
	echo "${dir2}/${samplename}.cuff/transcripts.gtf" >> ${dir2}/assembly_GTF_list.txt

	# Cleaning commands
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBcufflinks=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	SBcufflinksIDs=${SBcufflinksIDs}:${SBcufflinks##* }
	echo -e "\t Cufflinks job queued"

done < ${dir2}/Fastqs



echo ""
echo "-- All sample jobs queued! --"
echo ""



######################
## Generate data matrix
job="matrix"
# We are no longer processing sample files, but rsem files
samplename="rsem"

# variables
declare -a samplesarray=(`echo ${groupedsamples} | sed 's/,/ /g'`)
declare -a labelsarray=(`echo ${labels} | sed 's/,/ /g'`)
replicates=$(( ${#samplesarray[@]} / ${#labelsarray[@]} ))

# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=5:00" >> ${dir2}/${samplename}_${job}.sbatch
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
echo "#SBATCH --dependency=afterok${SBrsemIDs}" >> ${dir2}/${samplename}_${job}.sbatch

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi
	
# Job specific commands
mkdir -p ${dir2}/DESeq2

# Generate comparison matrices with all samples, will probably need to customize R script afterwards, only if more than 2 conditions
if [ ${#labelsarray[@]} -lt 2 ]
then
	condN="0"
	echo "rsem-generate-data-matrix ${dir2}/`echo ${groupedsamples} | sed \"s# #.rsem.genes.results ${dir2}/#g;s#,#.rsem.genes.results ${dir2}/#g\"`.rsem.genes.results > ${dir2}/DESeq2/`basename ${dir2}`_all.genes.matrix || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
	echo '#!/usr/bin/env Rscript' > ${dir2}/DESeq2/`basename ${dir2}`_all.DESeq2.R
	chmod +x ${dir2}/DESeq2/`basename ${dir2}`_all.DESeq2.R
	echo '# variables' >> ${dir2}/DESeq2/`basename ${dir2}`_all.DESeq2.R
	echo "table <- \"${dir2}/DESeq2/`basename ${dir2}`_all.genes.matrix\"" >> ${dir2}/DESeq2/`basename ${dir2}`_all.DESeq2.R
	while [ ${condN} -lt ${#labelsarray[@]} ]
	do
		echo "cond$(( $condN + 1 )) <- \"${labelsarray[${condN}]}\"" >> ${dir2}/DESeq2/`basename ${dir2}`_all.DESeq2.R
		conditions="`if [ -n "${conditions}" ]; then echo "${conditions}, "; fi`rep(cond$(( $condN + 1 )), ncol(data)/${#labelsarray[@]})"
		condN=$(( $condN + 1 ))
	done
	echo "output <- \"${dir2}/DESeq2/`basename ${dir2}`_\"" >> ${dir2}/DESeq2/`basename ${dir2}`_all.DESeq2.R
fi

# Generate pair-wise comparison matrices
condA="0"
condB="0"
while [ ${condA} -lt ${#samplesarray[@]} ]
do
	condB=$((condA + replicates ))
	while [ ${condB} -lt ${#samplesarray[@]} ]
	do
		echo "rsem-generate-data-matrix ${dir2}/`echo ${samplesarray[@]:${condA}:${replicates}} ${samplesarray[@]:${condB}:${replicates}} | sed \"s# #.rsem.genes.results ${dir2}/#g\"`.rsem.genes.results > ${dir2}/DESeq2/${labelsarray[$((condA / replicates))]}"_vs_"${labelsarray[$((condB / replicates))]}.genes.matrix || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
		echo ${dir2}/DESeq2/${labelsarray[$((condA / replicates))]}_vs_${labelsarray[$((condB / replicates))]}.genes.matrix >> ${dir2}/matrices
		echo '#!/usr/bin/env Rscript' > ${dir2}/DESeq2/${labelsarray[$((condA / replicates))]}_vs_${labelsarray[$((condB / replicates))]}.DESeq2.R
		chmod +x ${dir2}/DESeq2/${labelsarray[$((condA / replicates))]}_vs_${labelsarray[$((condB / replicates))]}.DESeq2.R
		echo '' >> ${dir2}/DESeq2/${labelsarray[$((condA / replicates))]}_vs_${labelsarray[$((condB / replicates))]}.DESeq2.R
		echo '# variables' >> ${dir2}/DESeq2/${labelsarray[$((condA / replicates))]}_vs_${labelsarray[$((condB / replicates))]}.DESeq2.R
		echo "table <- \"${dir2}/DESeq2/${labelsarray[$((condA / replicates))]}_vs_${labelsarray[$((condB / replicates))]}.genes.matrix\"" >> ${dir2}/DESeq2/${labelsarray[$((condA / replicates))]}_vs_${labelsarray[$((condB / replicates))]}.DESeq2.R
		echo "cond1 <- \"${labelsarray[$((condA / replicates))]}\"" >> ${dir2}/DESeq2/${labelsarray[$((condA / replicates))]}_vs_${labelsarray[$((condB / replicates))]}.DESeq2.R
		echo "cond2 <- \"${labelsarray[$((condB / replicates))]}\"" >> ${dir2}/DESeq2/${labelsarray[$((condA / replicates))]}_vs_${labelsarray[$((condB / replicates))]}.DESeq2.R
		echo "output <- \"${dir2}/DESeq2/\"" >> ${dir2}/DESeq2/${labelsarray[$((condA / replicates))]}_vs_${labelsarray[$((condB / replicates))]}.DESeq2.R
		condB=$((condB + replicates))
	done
	condA=$((condA + replicates))
done

# Cleaning commands
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBmatrix=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo -e "-- Matrix generation job queued --"



######################
## DESeq2 analysis and basic plots
job="deseq2"

echo ""
echo "-- Queuing DESeq2 jobs --"
echo ""

# Create R script for all sample matrix, only if more than 2 conditions
if [ ${#labelsarray[@]} -lt 2 ]
then
	# variables
	samplename="`basename ${dir2}`_all"
	Rscript=${samplename}.DESeq2.R

	# General SLURM parameters
	echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --time=5:00" >> ${dir2}/${samplename}_${job}.sbatch
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
	echo "#SBATCH --dependency=afterok:${SBmatrix##* }" >> ${dir2}/${samplename}_${job}.sbatch

	# General commands
	echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# Job specific commands
	cat >> ${dir2}/DESeq2/${Rscript} <<RALLDELIM


# This script uses DESeq2 in R to run differential-expression analysis on a reads count matrix of genes and automatically generate some graphs.


# Load libraries
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(calibrate)
library(pheatmap)

# Read in the data and make sample information table
data=read.table(table, header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)		# read the reads count matrix file that was generated by rsem-generate-data-matrix
colnames(data) <- sub('.rsem.genes.results', '', basename(colnames(data)))						# fix sample names in column names
cols=c(1:ncol(data))																			# get the number of columns
data[,cols]=apply(data[,cols], 2, function(x) as.numeric(as.integer(x)))						# format the matrix as integer numbers for DESeq2 to be happy
conditions <- factor(c(${conditions}))	# generate a comma-separated list of conditions/treatments for each sample. Replicates should be named the same.
samples=as.data.frame((colnames(data)))															# make samples list for DESeq2
samples <- cbind(samples, conditions)
colnames(samples)=c("sample","condition") 
groups=factor(samples[,2])

# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromMatrix(countData=as.matrix(data), colData=samples, design =~condition)		# generate a DESeq2 object from the data
colnames(dds) <- colnames(data)
dds <- DESeq(dds)																				# run DESeq2

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)


# Save differential expression results
res <- results(dds)
table(res\$padj<0.05)
## Order by adjusted p-value
res <- res[order(res\$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
## Write results
write.csv(resdata, file=paste0(output, "all-ExpDiff.csv"))


# Principal Components Analysis
pdf(paste0(output, "all-PCA.pdf"))
plotPCA(rld, intgroup="condition")
garbage <- dev.off() # Save to file


# Sample Distance Matrix
condcol <- brewer.pal(8, "Set1")[1:length(unique(conditions))]
dmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
sampleDists <- as.matrix(dist(t(assay(rld))))

pdf(paste0(output, "all-DistanceMatrix.pdf"))
heatmap.2(as.matrix(sampleDists), key=FALSE, symm=TRUE, trace="none", col=rev(dmcol), ColSideColors=condcol[conditions], RowSideColors=condcol[conditions], margin=c(15, 15), main="all samples Matrix")
garbage <- dev.off() # Save to file


# Heatmap of significant hits ( padj<0.05 and |log2FoldChange|>=1 )
hmcol <- colorRampPalette(brewer.pal(9, 'RdYlBu'))(100)
counts <- counts(dds,normalized=TRUE)															# Get normalized read counts (sorted by padj)
counts <- counts[apply(counts, 1, function(row) all(row !=0 )),]								# remove genes with zero reads
## You will probably need to change cond1 and cond2 to perform a meaningful comparison here
respairwise = results(dds, contrast=c("condition",cond2,cond4))									# Select comparison to perform (for log2 changes)
sig <- rownames(respairwise[!is.na(respairwise\$padj) & respairwise\$padj<0.05 & abs(respairwise\$log2FoldChange)>=1,])[1:30]		# Select the first 30 significant hits (sorted by padj)
sig <- sig[!is.na(sig)]					
sigcounts <- counts(dds,normalized=TRUE)[sig,]													# Get normalized read counts (sorted by padj) for significan genes
try({																							# May fail if only 0 or 1 gene is significant
	sigcounts <- sigcounts[apply(sigcounts, 1, function(row) all(row !=0 )),]					# remove genes with zero reads
})

# By defaults hits are not clustered and thus stay sorted by their padj value
pdf(paste0(output, "all-HeatmapSig.pdf"), onefile=FALSE)
try({																							# May fail if only 0 or 1 gene is significant
	pheatmap(log2(sigcounts), col=rev(hmcol), scale='row', cluster_rows=FALSE, cluster_cols=FALSE, main="all Top Hits")
})
garbage <- dev.off() # Save to file

# But we can cluster them using the following command
pdf(paste0(output, "all-HeatmapSigClust.pdf"), onefile=FALSE)
try({																							# May fail if only 0 or 1 gene is significant
	pheatmap(log2(sigcounts), col=rev(hmcol), scale='row', cluster_rows=TRUE, cluster_cols=TRUE, main="all Clustered Top Hits")
})
garbage <- dev.off() # Save to file

# And even plot everything, this time using the old school Red/Green color palette
pdf(paste0(output, "all-HeatmapAllClust.pdf"), onefile=FALSE)
pheatmap(log2(counts), col=redgreen(75), scale='row', cluster_rows=TRUE, cluster_cols=TRUE, main="all Clustered Reads Count")
garbage <- dev.off() # Save to file


# Adapted from:
# http://www.bioconductor.org/help/workflows/rnaseqGene/
# https://gist.github.com/stephenturner/f60c1934405c127f09a6
# https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/

RALLDELIM
	# Use Rscript
	echo "Rscript ${dir2}/DESeq2/${Rscript}" >> ${dir2}/${samplename}_${job}.sbatch

	# Cleaning commands
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBdeseq2=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	SBdeseq2IDs=${SBdeseq2IDs}:${SBdeseq2##* }
	echo -e "\t ${samplename} DESeq2 job queued"
fi

# Create R scripts for pair wise matrices
while read line;
do

	# variables
	samplename=`basename ${line} | sed 's/.genes.matrix//g'`
	Rscript=${samplename}.DESeq2.R

	# General SLURM parameters
	echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
	echo "#SBATCH --time=5:00" >> ${dir2}/${samplename}_${job}.sbatch
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
	echo "#SBATCH --dependency=afterok:${SBmatrix##* }" >> ${dir2}/${samplename}_${job}.sbatch

	# General commands
	echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# Job specific commands
	cat >> ${dir2}/DESeq2/${Rscript} <<RDELIM


# This script uses DESeq2 in R to run differential-expression analysis on a reads count matrix of genes and automatically generate some graphs.


# Load libraries
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(calibrate)
library(pheatmap)

# Read in the data and make sample information table
data=read.table(table, header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)		# read the reads count matrix file that was generated by rsem-generate-data-matrix
colnames(data) <- sub('.rsem.genes.results', '', basename(colnames(data)))						# fix sample names in column names
cols=c(1:ncol(data))																			# get the number of columns
data[,cols]=apply(data[,cols], 2, function(x) as.numeric(as.integer(x)))						# format the matrix as integer numbers for DESeq2 to be happy
conditions <- factor(c(rep(cond1, ncol(data)/2), rep(cond2, ncol(data)/2)))						# generate a comma-separated list of conditions/treatments for each sample. Replicates should be named the same.
samples=as.data.frame((colnames(data)))															# make samples list for DESeq2
samples <- cbind(samples, conditions)
colnames(samples)=c("sample","condition") 
groups=factor(samples[,2])

# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromMatrix(countData=as.matrix(data), colData=samples, design =~condition)		# generate a DESeq2 object from the data
colnames(dds) <- colnames(data)
dds <- DESeq(dds)																				# run DESeq2

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)


# Save differential expression results
res <- results(dds)
table(res\$padj<0.05)
## Order by adjusted p-value
res <- res[order(res\$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
## Write results
write.csv(resdata, file=paste0(output, cond1, "_vs_", cond2, "-ExpDiff.csv"))


# Plot dispersions
pdf(paste0(output, cond1, "_vs_", cond2, "-DispersionPlot.pdf"))
plotDispEsts(dds, main=paste0(cond1, " vs ", cond2, " Dispersion Plot"))
garbage <- dev.off() # Save to file


# Principal Components Analysis
pdf(paste0(output, cond1, "_vs_", cond2, "-PCA.pdf"))
plotPCA(rld, intgroup="condition")
garbage <- dev.off() # Save to file


# MA plot (without and with "significant" genes labeled)
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}

pdf(paste0(output, cond1, "_vs_", cond2, "-MAplot.pdf"))
plotMA(res, ylim=c(-2,2), main=paste0(cond1, " vs ", cond2, " MA Plot"))																# DESeq2 MA plot function
maplot(resdata, xlab="mean of normalized counts", ylab=expression(log[2]~fold~change), main=paste0(cond1, " vs ", cond2, " MA Plot"))	# Custom MA with gene names
garbage <- dev.off() # Save to file


# Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
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

pdf(paste0(output, cond1, "_vs_", cond2, "-VolcanoPlot.pdf"))
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2, 2), main=paste0(cond1, " vs ", cond2, " Volcano Plot"))
garbage <- dev.off() # Save to file


# Sample Distance Matrix
condcol <- brewer.pal(8, "Set1")[1:length(unique(conditions))]
dmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
sampleDists <- as.matrix(dist(t(assay(rld))))

pdf(paste0(output, cond1, "_vs_", cond2, "-DistanceMatrix.pdf"))
heatmap.2(as.matrix(sampleDists), key=FALSE, symm=TRUE, trace="none", col=rev(dmcol), ColSideColors=condcol[conditions], RowSideColors=condcol[conditions], margin=c(15, 15), main=paste0(cond1, " vs ", cond2, " Matrix"))
garbage <- dev.off() # Save to file

# Heatmap of significant hits ( padj<0.05 and |log2FoldChange|>=1 )
hmcol <- colorRampPalette(brewer.pal(9, 'RdYlBu'))(100)
sig <- rownames(res[!is.na(res\$padj) & res\$padj<0.05 & abs(res\$log2FoldChange)>=1,])[1:30]	# Select the first 30 significant hits (sorted by padj)
sig <- sig[!is.na(sig)]																			# If we have less than 30 hits, clean "NA" fields

counts <- counts(dds,normalized=TRUE)															# Get normalized read counts (sorted by padj)
counts <- counts[apply(counts, 1, function(row) all(row !=0 )),]								# remove genes with zero reads
sigcounts <- counts(dds,normalized=TRUE)[sig,]													# Get normalized read counts (sorted by padj) for significant genes
try({																							# May fail if only 0 or 1 gene is significant
sigcounts <- sigcounts[apply(sigcounts, 1, function(row) all(row !=0 )),]						# remove genes with zero reads
})

# By defaults hits are not clustered and thus stay sorted by their padj value
pdf(paste0(output, cond1, "_vs_", cond2, "-HeatmapSig.pdf"), onefile=FALSE)
try({																							# May fail if only 0 or 1 gene is significant
pheatmap(log2(sigcounts), col=rev(hmcol), scale='row', cluster_rows=FALSE, cluster_cols=FALSE, main=paste0(cond1, " vs ", cond2, " Top Hits"))
})
garbage <- dev.off() # Save to file

# But we can cluster them using the following command
pdf(paste0(output, cond1, "_vs_", cond2, "-HeatmapSigClust.pdf"), onefile=FALSE)
try({																							# May fail if only 0 or 1 gene is significant
pheatmap(log2(sigcounts), col=rev(hmcol), scale='row', cluster_rows=TRUE, cluster_cols=TRUE, main=paste0(cond1, " vs ", cond2, " Clustered Top Hits"))
})
garbage <- dev.off() # Save to file

# And even plot everything, this time using the old school Red/Green color palette
pdf(paste0(output, cond1, "_vs_", cond2, "-HeatmapAllClust.pdf"), onefile=FALSE)
pheatmap(log2(counts), col=redgreen(75), scale='row', cluster_rows=TRUE, cluster_cols=TRUE, main=paste0(cond1, " vs ", cond2, " Clustered Reads Count"))
garbage <- dev.off() # Save to file


# Adapted from:
# http://www.bioconductor.org/help/workflows/rnaseqGene/
# https://gist.github.com/stephenturner/f60c1934405c127f09a6
# https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/

RDELIM
	# Use Rscript
	echo "Rscript ${dir2}/DESeq2/${Rscript}" >> ${dir2}/${samplename}_${job}.sbatch

	# Cleaning commands
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBdeseq2=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	SBdeseq2IDs=${SBdeseq2IDs}:${SBdeseq2##* }
	echo -e "\t ${samplename} DESeq2 job queued"

done < ${dir2}/matrices



echo ""
echo "-- All DESeq2 jobs queued! --"
echo ""



######################
## Merge transcripts
job="cuffmerge"
# We are no longer processing sample files, but GTF files
samplename="GTFs"
	
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
echo "#SBATCH --dependency=afterok${SBcufflinksIDs}" >> ${dir2}/${samplename}_${job}.sbatch

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi
	
# Job specific commands
#echo "cuffmerge -p ${threads} -g ${gtf}/genes.gtf -s ${fasta_refgenome} -o ${dir2} ${dir2}/assembly_GTF_list.txt" >> ${dir2}/${samplename}_${job}.sbatch
echo "cuffmerge -p ${threads} -g ${gtf}/genes.gtf -s ${fasta_refgenome} -o ${tmp} ${dir2}/assembly_GTF_list.txt || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
echo "cp ${tmp}/merged.gtf ${dir2}/merged.gtf" >> ${dir2}/${samplename}_${job}.sbatch

# Cleaning commands
echo "[ -f ${dir2}/assembly_GTF_list.txt ] && rm ${dir2}/assembly_GTF_list.txt" >> ${dir2}/${samplename}_${job}.sbatch
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBcuffmerge=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo -e "-- Cuffmerge job queued --"



######################
## Quantify transcripts
job="cuffquant"

echo ""
echo "-- Queuing quantification jobs --"
echo ""

while read line;
do

	# General variables
	read1=`echo ${line} | cut -d" " -f1`
	read2=`echo ${line} | cut -d" " -f2`
	samplename=`echo ${read1} | awk -F_R1 '{print $1}'`
	bamsortedout=`basename ${read1} | sed "s/_R1${fileext}/.sorted.bam/g"`

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
	echo "#SBATCH --dependency=afterok:${SBcuffmerge##* }" >> ${dir2}/${samplename}_${job}.sbatch
	
	# General commands
	echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
	if [ -n "${customcmd}" ]
	then
		echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
	fi
	
	# Job specific commands
	echo "cuffquant -p ${threads} -o ${tmp}/${samplename} ${dir2}/merged.gtf ${dir2}/${bamsortedout} || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
	echo "cp -rf ${tmp}/${samplename}/abundances.cxb ${dir2}/${samplename}.cuff/abundances.cxb" >> ${dir2}/${samplename}_${job}.sbatch

	# Cleaning commands
	# remove .sbatch
	echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

	# Queue job
	SBcuffquant=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
	SBcuffquantIDs=${SBcuffquantIDs}:${SBcuffquant##* }
	echo -e "\t ${samplename} cuffquant job queued"

done < ${dir2}/Fastqs



echo ""
echo "-- All quantification jobs queued! --"
echo ""



######################
## Compare expression levels
job="cuffdiff"
# We are no longer processing sample files, but CXBs files
samplename="CXBs"

# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=2:00:00" >> ${dir2}/${samplename}_${job}.sbatch
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
echo "#SBATCH --dependency=afterok${SBcuffquantIDs}" >> ${dir2}/${samplename}_${job}.sbatch
	
# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi
	
# Job specific commands
echo "[ -d ${dir2}/ExpDiff ] && rm -rf ${dir2}/ExpDiff" >> ${dir2}/${samplename}_${job}.sbatch
#echo "cuffdiff -u -p ${threads} -b ${fasta_refgenome} ${dir2}/merged.gtf ${dir2}/`echo ${groupedsamples} | sed "s#,#.cuff/abundances.cxb,${dir2}/#g;s# #.cuff/abundances.cxb ${dir2}/#g"`.cuff/abundances.cxb -L ${labels} -o ${dir2}/ExpDiff || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
echo "cuffdiff -u -p ${threads} -b ${fasta_refgenome} ${dir2}/merged.gtf ${dir2}/`echo ${groupedsamples} | sed "s#,#.cuff/abundances.cxb,${dir2}/#g;s# #.cuff/abundances.cxb ${dir2}/#g"`.cuff/abundances.cxb -L ${labels} -o ${tmp}/ExpDiff || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
echo "cp -rf ${tmp}/ExpDiff ${dir2}/ExpDiff" >> ${dir2}/${samplename}_${job}.sbatch

# Cleaning commands
echo "[ -f ${dir2}/assembly_GTF_list.txt ] && rm ${dir2}/assembly_GTF_list.txt" >> ${dir2}/${samplename}_${job}.sbatch
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBcuffdiff=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo -e "-- Cuffdiff job queued --"



######################
## Compare expression levels
job="cuffnorm"

# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --mem=${mem}000 --cpus-per-task=${threads}" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=1:00:00" >> ${dir2}/${samplename}_${job}.sbatch
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
echo "#SBATCH --dependency=afterok${SBcuffquantIDs}" >> ${dir2}/${samplename}_${job}.sbatch
	
# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi
	
# Job specific commands
echo "[ -d ${dir2}/Normalized ] && rm -rf ${dir2}/Normalized" >> ${dir2}/${samplename}_${job}.sbatch
#echo "cuffnorm -p ${threads} ${dir2}/merged.gtf ${dir2}/`echo ${groupedsamples} | sed "s#,#.cuff/abundances.cxb,${dir2}/#g;s# #.cuff/abundances.cxb ${dir2}/#g"`.cuff/abundances.cxb -L ${labels} -o ${dir2}/Normalized || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
echo "cuffnorm -p ${threads} ${dir2}/merged.gtf ${dir2}/`echo ${groupedsamples} | sed "s#,#.cuff/abundances.cxb,${dir2}/#g;s# #.cuff/abundances.cxb ${dir2}/#g"`.cuff/abundances.cxb -L ${labels} -o ${tmp}/Normalized || if [ -f ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err ]; then echo \${SLURM_JOB_NODELIST} >> ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && exit 1; else echo \${SLURM_JOB_NODELIST} > ${dir2}/\${SLURM_JOBID}-${samplename}_${job}.err && scontrol requeue \${SLURM_JOBID} && sleep 42m; fi" >> ${dir2}/${samplename}_${job}.sbatch
echo "cp -rf ${tmp}/Normalized ${dir2}/Normalized" >> ${dir2}/${samplename}_${job}.sbatch

# Cleaning commands
echo "[ -f ${dir2}/assembly_GTF_list.txt ] && rm ${dir2}/assembly_GTF_list.txt" >> ${dir2}/${samplename}_${job}.sbatch
# remove .sbatch
echo "rm ${dir2}/${samplename}_${job}.sbatch" >> ${dir2}/${samplename}_${job}.sbatch

# Queue job
SBcuffnorm=$(sbatch ${dir2}/${samplename}_${job}.sbatch)
echo -e "-- Cuffnorm job queued --"



######################
## Final processing of result files
job="notif"
# We are no longer processing sample files
samplename="`basename ${dir2}`-RNAseq"

# General SLURM parameters
echo '#!/bin/bash' > ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --job-name=${samplename}_${job} --output=${dir2}/${logs}/${samplename}_${job}.out --error=${dir2}/${logs}/${samplename}_${job}.err --open-mode=append" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH `if [ ! -z ${threads} ] && [ ${threads} -gt "1" ]; then echo "--cpus-per-task=1"; else echo "--cpus-per-task=${threads}"; fi`" >> ${dir2}/${samplename}_${job}.sbatch
echo "#SBATCH --time=3:00" >> ${dir2}/${samplename}_${job}.sbatch
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
echo "#SBATCH --dependency=afterok:${SBcuffdiff##* }:${SBcuffnorm##* }`if [ -n \"${SBfqcIDs}\" ]; then echo \"${SBfqcIDs}\"; fi``if [ -n \"${SBmatrix}\" ]; then echo \":${SBmatrix##* }\"; fi``if [ -n \"${SBdeseq2IDs}\" ]; then echo \"${SBdeseq2IDs}\"; fi`" >> ${dir2}/${samplename}_${job}.sbatch

# General commands
echo "mkdir -p ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${customcmd}" ]
then
	echo "${customcmd}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Job specific commands
# Creating folders for archive
echo "[ -f ${dir2}/Results.tar.gz ] && rm ${dir2}/Results.tar.gz" >> ${dir2}/${samplename}_${job}.sbatch
echo "[ -d ${dir2}/Results_`basename ${dir2}` ] && rm -rf ${dir2}/Results_`basename ${dir2}`" >> ${dir2}/${samplename}_${job}.sbatch
echo "mkdir -p ${dir2}/Results_`basename ${dir2}`/" >> ${dir2}/${samplename}_${job}.sbatch
# List all csv, html, pdf and zip files
echo "csvarchive=\`ls ${dir2}/*.csv\`" >> ${dir2}/${samplename}_${job}.sbatch
echo "htmlarchive=\`ls ${dir2}/*.html\`" >> ${dir2}/${samplename}_${job}.sbatch
echo "pdfarchive=\`ls ${dir2}/*.pdf\`" >> ${dir2}/${samplename}_${job}.sbatch
echo "ziparchive=\`ls ${dir2}/*.zip\`" >> ${dir2}/${samplename}_${job}.sbatch
# Create an archive of all csv, html, pdf and zip files for easy download
echo "for i in \${csvarchive}; do cp -rf \$i ${dir2}/Results_`basename ${dir2}`/\`basename \$i\`; done" >> ${dir2}/${samplename}_${job}.sbatch
echo "for i in \${htmlarchive}; do cp -rf \$i ${dir2}/Results_`basename ${dir2}`/\`basename \$i\`; done" >> ${dir2}/${samplename}_${job}.sbatch
echo "for i in \${pdfarchive}; do cp -rf \$i ${dir2}/Results_`basename ${dir2}`/\`basename \$i\`; done" >> ${dir2}/${samplename}_${job}.sbatch
echo "for i in \${ziparchive}; do cp -rf \$i ${dir2}/Results_`basename ${dir2}`/\`basename \$i\`; done" >> ${dir2}/${samplename}_${job}.sbatch
echo "cp -rf ${dir2}/ExpDiff ${dir2}/Results_`basename ${dir2}`/" >> ${dir2}/${samplename}_${job}.sbatch
echo "cp -rf ${dir2}/Normalized ${dir2}/Results_`basename ${dir2}`/" >> ${dir2}/${samplename}_${job}.sbatch
echo "cp -rf ${dir2}/DESeq2 ${dir2}/Results_`basename ${dir2}`/" >> ${dir2}/${samplename}_${job}.sbatch
echo "cp -rf ${dir2}/merged.gtf ${dir2}/Results_`basename ${dir2}`/" >> ${dir2}/${samplename}_${job}.sbatch
echo "tar --remove-files -C ${dir2} -pczf ${dir2}/Results.tar.gz Results_`basename ${dir2}`" >> ${dir2}/${samplename}_${job}.sbatch
if [ -n "${iftttkey}" ] && [ -n "${iftttevent}" ]
then
	# Trigger IFTTT maker channel event when it's ready, nice isn't it?
	echo "curl -X POST -H \"Content-Type: application/json\" -d '{ \"value1\" : \"`basename ${dir2}`\" }' https://maker.ifttt.com/trigger/${iftttevent}/with/key/${iftttkey}" >> ${dir2}/${samplename}_${job}.sbatch
fi

# Cleaning commands
# Remove Temporary directory
echo "rm -rf ${tmp}" >> ${dir2}/${samplename}_${job}.sbatch
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
[ -f ${dir2}/matrices ] && rm ${dir2}/matrices

# That's all folks!
exit
