Bioinformatics
==============

High throughput sequencing scripts: bowtie2, GATK, etc...

## sh_FastQToSNPsCall.sh
usage: sh_FastQToSNPsCall.sh </path/to/fastq(.gz)/folder> </path/to/Aligned(.bam)/destination/folder> </path/to/SNPsCalled/folder>

# Description
Will call sh_bowtie2_AlignAll.sh to convert fastq to .sam then sh_samtools_ProcessSams.sh to convert .sam into sorted and indexed .bam
Finally will launch sh_gatkSNPcalling.sh to call SNPs with GATK and annotate them with ANNOVAR.
