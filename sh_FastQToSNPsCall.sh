#!/bin/bash
#
# Usage: sh_FastQToSNPsCall.sh </path/to/fastq(.gz)/folder> </path/to/Aligned(.bam)/destination/folder> </path/to/SNPsCalled/folder> [/path/to/config/file.ini]
#
## Description ##
#
# Will call sh_bowtie2_AlignAll.sh to convert fastq to .sam
# then sh_samtools_ProcessSams.sh to convert .sam into sorted and indexed .bam
# Finally will launch sh_gatkSNPcalling.sh to call SNPs with GATK and annotate them with ANNOVAR.
#
##

if [ -z "$1" -o -z "$2" -o -z "$3" ]
then
    echo "Usage: sh_FastQToSNPsCall.sh </path/to/fastq(.gz)/folder> </path/to/Aligned(.bam)/destination/folder> </path/to/SNPsCalled/folder> [/path/to/config/file.ini]"
    exit
fi

if [ -n "$4" ] && [ "${4: -4}" != ".ini" ]
then
    echo "Invalid config file detected. Is it an .ini file?"
    echo "Usage: sh_FastQToSNPsCall.sh </path/to/fastq(.gz)/folder> </path/to/Aligned(.bam)/destination/folder> </path/to/SNPsCalled/folder> [/path/to/config/file.ini]"
    exit
fi

sh_bowtie2_AlignAll.sh "$1" "$2" "$4"
sh_samtools_ProcessSams.sh "$2" "$4"
sh_gatkSNPcalling.sh "$2" "$3" "$4"
