# Bioinformatics pipelines, [SLURM](https://slurm.schedmd.com/overview.html) friendly

![GitHub package.json version](https://img.shields.io/github/package-json/v/emc2cube/Bioinformatics)
![GitHub top language](https://img.shields.io/github/languages/top/emc2cube/Bioinformatics?color=green)
![GitHub](https://img.shields.io/github/license/emc2cube/Bioinformatics?color=yellow)
[![Runs on Sherlock](https://img.shields.io/badge/Runs_on-Sherlock-red)](https://www.sherlock.stanford.edu)

> Set of high throughput sequencing analysis scripts to quickly generate and queue jobs on [SLURM](https://slurm.schedmd.com/overview.html)-based HPC clusters, such as [Stanford's Sherlock](https://www.sherlock.stanford.edu)🕵🏻‍♂️️
>
> Most scripts include some sort of failsafe: if a job fails it will be requeued once. This is useful in case of unexpected node failure.
>
> Currently available pipelines:
> * Whole Exome Sequencing
> * RNA Sequencing
> * CRISPR screens


## sh_WES.sh

This script will process fastq(.gz) files and align them to a reference genome using bowtie2.
It will then use Picard and GATK following the June 2016 best practices workflow.
SNPs will then be annotated using ANNOVAR.

See the [WES.ini](https://github.com/emc2cube/Bioinformatics/blob/master/config_WES.ini) configuration file for all available options and settings.

Options:
* --help : Display help message.
* --version : Display version number.

### Usage

```sh
sh_WES.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]
```


## sh_RNAseq.sh

This script will process fastq(.gz) files and align them to a reference genome using either STAR (recommended), hishat2 or tophat2.
If STAR is used then RSEM will also be used and differential expression will be analyzed using DESeq2.
Differential expression can also be computed using cufflinks (cufflinks is pretty much deprecated, should be avoided unless trying to reproduce old results).

See the [RNAseq.ini](https://github.com/emc2cube/Bioinformatics/blob/master/config_RNAseq.ini) configuration file for all available options and settings.

Options:
* --help : Display help message.
* --version : Display version number.

### Usage

```sh
sh_RNAseq.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]
```


## sh_CRISPR.sh

This script will process the fastq(.gz) files generated in a typical CRISPR screen using either [casTLE](https://bitbucket.org/dmorgens/castle/) or [MAGeCK](https://sourceforge.net/projects/mageck/).
* If using casTLE, a reference file of all the indices will be automatically created using bowtie (NOT bowtie2). It will then analyze the screen and generate basic graphs.
* If using MAGeCK counts, tests, mle and pathway analysis will be performed. It will also run the [R](https://www.r-project.org) package "[MAGeCKFlute](https://bioconductor.org/packages/release/bioc/html/MAGeCKFlute.html)" and in all cases generate basic graphs.

See the [CRISPR.ini](https://github.com/emc2cube/Bioinformatics/blob/master/config_CRISPR.ini) configuration file for all available options and settings.

Options:
* --help : Display help message.
* --version : Display version number.

Dependancies:
[csvkit](https://csvkit.readthedocs.io/en/latest/) should be installed on your system in a location included in your $PATH.


### Usage

```sh
sh_CRISPR.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]
```

### Python 3.6 compatibility

For easy integration along MAGeCK, or any other modern tools, a python 3.6+ compatible version of casTLE is included.
This is based on [casTLE commit 981d6d8](https://bitbucket.org/dmorgens/castle/commits/981d6d877c0fe3ee233e9fd977b13800987a032c) and may not be up to date.
You still need to download the whole [casTLE repository](https://bitbucket.org/dmorgens/castle/) even if you end up switching the scripts with their python 3.6+ compatible version.


## sh_md5alldir.sh

This script will process all sub-directories of the input folders and for each of them will create a <directory_name>.md5 file if it does not exist yet, or check <directory> files against the existing <directory_name>.md5 file.

Options:
* -f or --force : even if there is already a <directory>.md5 file, it will be replaced by a new <directory>.md5 file.
* --help : Display help message.
* --version : Display version number.

### Usage

```sh
sh_md5alldir.sh </path/to/dir/> [OPTIONS]
```


## sh_sha1alldir.sh

This script will process all sub-directories of the input folders and for each of them will create a <directory_name>.sha1 file if it does not exist yet, or check <directory> files against the existing <directory_name>.sha1 file.

Options:
* -f or --force : even if there is already a <directory>.sha1 file, it will be replaced by a new <directory>.sha1 file.
* --help : Display help message.
* --version : Display version number.

### Usage

```sh
sh_sha1alldir.sh </path/to/dir/> [OPTIONS]
```


## sh_ACMGfilter.sh

This script will look for an annovar .snps.exome_summary.csv file and generate a list of all SNPs found in the ACMG guidelines in a new ACMG_genes.csv file.
This file can be directly sent to a clinician for incidental findings reports, if required.

Options:
* --help : Display help message.
* --version : Display version number.

### Usage

```sh
sh_ACMGfilter.sh </path/to/.csv/containing/folder> [/path/to/destination/folder]
```


## sh_mergeFastQ.sh

Simple script to consolidate fragmented .fastq files from different sequencing lanes.
Original files will be backed up in a FastQbackup folder.

Options:
* --help : Display help message.
* --version : Display version number.

### Usage

```sh
sh_mergeFastQ.sh </path/to/fastq(.gz)/folder>
```


## Author(s) contributions

👤 **Julien Couthouis**

*Initial work and releases*

* Linkedin: [@jcouthouis](https://www.linkedin.com/in/jcouthouis/)
* Github: [@emc2cube](https://github.com/emc2cube)


## Show your support

Give a ![GitHub stars](https://img.shields.io/github/stars/emc2cube/Bioinformatics?style=social) if this project helped you!


## License

Copyright © 2019 [Julien Couthouis](https://github.com/emc2cube).

This project is [EUPL-1.2](https://github.com/emc2cube/Bioinformatics/blob/master/LICENSE) licensed.