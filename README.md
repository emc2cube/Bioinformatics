Bioinformatics
==============

High throughput sequencing scripts: bowtie2, GATK, etc...


Workflow scripts:
-----------------


## sh_bowtie2_AlignAll.sh

Usage: sh_bowtie2_AlignAll.sh </path/to/fastq(.gz)/folder> </path/to/Aligned(.bam)/destination/folder> [/path/to/config/file.ini]

# Description

This script will process fastq files to .sam files
First it will trim fastqs with Jason script to remove sequences in the adapter region.
It will also generate a merged file of all individual samples to be aligned.
Finally trimmed sequences are aligned to a reference genome using bowtie2


## sh_samtools_ProcessSams.sh

Usage: sh_samtools_ProcessSams.sh </path/to/sam/files/> [/path/to/config/file.ini]

# Description

This script will process .sam files and convert them to .bam files.
Optional: reads can be filtered on their MapQ score to remove poorly aligned reads.
These .bam files are then sorted and indexed.
Intermediate files are deleted to save space.


## sh_gatkSNPcalling.sh

Usage: sh_gatkSNPcalling.sh </path/to/Aligned(.bam)/destination/folder> </path/to/SNPsCalled/folder> [/path/to/config/file.ini]

# Description

This script will process .bam files
Optional: Will first remove duplicate reads (edit options).
This script will, for all samples:
- perform a local realignment around known indels.
- perform a quality score recalibration.
- call SNPs.
- annotate using annovar
- merge csv files and do some cleaning for an easy downloadable file
As this script will be the last one from the workflow, an option exist to send an email when the job is complete.


## sh_FastQToSNPsCall.sh

Usage: sh_FastQToSNPsCall.sh </path/to/fastq(.gz)/folder> </path/to/Aligned(.bam)/destination/folder> </path/to/SNPsCalled/folder> [/path/to/config/file.ini]

# Description

Will call sh_bowtie2_AlignAll.sh to convert fastq to .sam
then sh_samtools_ProcessSams.sh to convert .sam into sorted and indexed .bam
Finally will launch sh_gatkSNPcalling.sh to call SNPs with GATK and annotate them with ANNOVAR.


## sh_FastQToTrimmedBam.sh

Usage: sh_FastQToTrimmedBam.sh </path/to/fastq(.gz)/folder> [/path/to/destination/folder] [/path/to/config/file.ini]

# Description

Will call sh_bowtie2_AlignAll.sh to convert fastq to .sam
then sh_samtools_ProcessSams.sh to convert .sam into sorted and indexed .bam


## sh_filterSNPs.sh

Usage: sh_filterSNPs.sh </path/to/.vcf/files> </path/to/save/filtered.vcf/files> [/path/to/config/file.ini]

# Description

This script will process GATK .vcf files
This script will, for all samples:
- filter SNPs.
- annotate using annovar
- merge csv files and do some cleaning for an easy downloadable file
This function is already included in sh_gatkSNPcalling.sh


## sh_gatkCoverage.sh

Usage: sh_gatkCoverage.sh </path/to/Aligned(.bam)/destination/folder> </path/to/SNPsCalled/folder> [/path/to/config/file.ini]

# Description

This script will process .bam files and compute coverage information.
This function is already included in sh_gatkSNPcalling.sh



Utilities scripts:
------------------


## sh_md5alldir.sh

Usage: sh_md5alldir.sh </path/to/dir/> [-options, -? or --help for help]

# Description

This script will process all sub-directories of the input folders and for each of them
will create a <directory_name>.md5 file if it does not exist yet, or check <directory> files
against the existing <directory_name>.md5 file.

# Options:

-f or --force : even if a <directory>.md5 file is detected, will replace it by a fresh one
and will not check files against it.


## sh_catFiles.sh

Usage: sh_catFiles.sh </path/to/original/fastq(.gz)/folder> </path/to/merged/fastq(.gz)/folder>

# Description

This script will merge all fastq(.gz) from a folder to fastq(.gz) with the same
name in a second folder.


## sh_csvmerge.sh

Usage: sh_csvmerge.sh </path/to/.csv/containing/folder> [/path/to/destination/folder]

# Description

This script will merge .snps.exome_summary.csv file to a new one, adding the sample
name in the first column, and fixing GATK header columns in the resulting ANNOVAR
.snps.exome_summary.csv file.
