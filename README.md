Bioinformatics
==============

High throughput sequencing scripts: bowtie2, GATK, etc...


Workflow scripts:
-----------------


## sh_WES.sh (SLURM compatible)

Usage: sh_WES.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]

# Description

This script will process fastq(.gz) files and align them to a reference genome using bowtie2.
It will then use Picard and GATK following GATK according to June 2016 best practices workflow.
SNPs will then be annotated with ANNOVAR.
Include a failsafe, if a job fails, it will be requeued once in case of a hardware failure.

# Options

Can call trimmomatic, FastQC and compute coverage.
Settings can be modified by using a customized config_WES.ini file.

## sh_RNAseq.sh (SLURM compatible)

Usage: sh_RNAseq.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]

# Description

This script will process fastq(.gz) files and align them to a reference genome using either STAR (recommended), hishat2 or tophat2.
Differential expression will then be computed using cufflinks.
If STAR is used then RSEM will also be used to generate gene read counts, pairwise comparison matrices will be created and DESeq2 analysis will be performed.
Include a failsafe, if a job fails, it will be requeued once in case of a hardware failure.

# Options

Can call trimmomatic and FastQC.
Settings can be modified by using a customized config_RNAseq.ini file.

## sh_bowtie2_AlignAll.sh (deprecated)

Usage: sh_bowtie2_AlignAll.sh </path/to/fastq(.gz)/folder> </path/to/Aligned(.bam)/destination/folder> [/path/to/config/file.ini]

# Description

This script will convert fastq files to bowtie2 aligned .bam files.
Optional: Can call a trimming program (trimmomatic, Trim Galore or your own script).
This script will, for all samples in input folder:
- convert .fastq or .fastq.gz files to .sam.
- align .sam to reference genome.
- convert .sam to .bam.
- Sort and index .bam file.

## sh_gatkSNPcalling.sh (deprecated)

Usage: sh_gatkSNPcalling.sh </path/to/Aligned(.bam)/destination/folder> </path/to/SNPsCalled/folder> [/path/to/config/file.ini]

# Description

This script will process aligned .bam files.
Optional: Will first remove duplicate reads.
This script will, for all samples:
- perform a local realignment around known indels.
- perform a quality score recalibration.
- generate .g.vcf file using HaplotypeCaller.
- Optional: stop here
- perform joint genotyping
- Filter variants using VQSR.
- annotate using annovar.
- do some cleaning on .csv for an easy downloadable file.
- Optional: Can trigger an IFTTT event using the maker channel.


## sh_FastQToSNPsCall.sh (deprecated)

Usage: sh_FastQToSNPsCall.sh </path/to/fastq(.gz)/folder> </path/to/Aligned(.bam)/destination/folder> </path/to/SNPsCalled/folder> [/path/to/config/file.ini]

# Description

Will call sh_bowtie2_AlignAll.sh to convert fastq to aligned .bam and then launch sh_gatkSNPcalling.sh to call SNPs with GATK and annotate them with ANNOVAR.


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
