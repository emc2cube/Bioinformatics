# Bioinformatics pipelines, [SLURM](https://slurm.schedmd.com/overview.html) friendly

![GitHub package.json version](https://img.shields.io/github/package-json/v/emc2cube/Bioinformatics)![GitHub top language](https://img.shields.io/github/languages/top/emc2cube/Bioinformatics?color=green)![GitHub](https://img.shields.io/github/license/emc2cube/Bioinformatics?color=yellow)[![Runs on Sherlock](https://img.shields.io/badge/Runs_on-Sherlock-red)](https://www.sherlock.stanford.edu)

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
It is currently only compatible with GATK version 3.X, and have been used extensively with the latest available version [GATK 3.8.1](https://software.broadinstitute.org/gatk/download/archive)
SNPs will then be annotated using ANNOVAR.

See the [WES.ini](https://github.com/emc2cube/Bioinformatics/blob/master/config_WES.ini) configuration file for all available options and settings.

##### Options:
* ```--help``` : Display help message.
* ```--version``` : Display version number.

##### Dependencies:
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) should be installed on your system in a location included in your $PATH for alignment.
* [Samtools](http://www.htslib.org) should be installed on your system in a location included in your $PATH.
* [Picard](http://broadinstitute.github.io/picard/) should be installed on your system.
* [GATK 3.X](https://software.broadinstitute.org/gatk/download/archive) should be installed on your system for transcript quantification.
* (optional) [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) should be installed on your system for read trimming.
* (optional) [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) should be installed on your system in a location included in your $PATH for quality control.

### Usage

```sh
sh_WES.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]
```


## sh_RNAseq.sh

This script will process fastq(.gz) files and align them to a reference genome using either STAR (recommended), hishat2 or tophat2.
If STAR is used then RSEM will also be used and differential expression will be analyzed using DESeq2.
Differential expression can also be computed using cufflinks (cufflinks is pretty much deprecated, should be avoided unless trying to reproduce old results).
Local Splicing Variation can now be computed using MAJIQ and/or LeafCutter.
If a 4th '[OutputDirName]' argument is provided only the secondary analyses selected in the config file will be queued, 6ft/2m apart, using the already aligned and processed files from a previous run, and results will be saved in a '_OutputDirName' directory.


See the [RNAseq.ini](https://github.com/emc2cube/Bioinformatics/blob/master/config_RNAseq.ini) configuration file for all available options and settings.

##### Options:
* ```--help``` : Display help message.
* ```--version``` : Display version number.

##### Dependencies:
* [STAR](https://github.com/alexdobin/STAR) (recommended) or [tophat2](https://ccb.jhu.edu/software/tophat/index.shtml) / [hisat2](https://daehwankimlab.github.io/hisat2/) should be installed on your system in a location included in your $PATH for alignment.
* [RSEM](https://github.com/deweylab/RSEM) should be installed on your system in a location included in your $PATH for transcript quantification.
* [R](https://www.r-project.org) should be installed on your system in a location included in your $PATH.
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) (recommended) or [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) (deprecated) for differential expression analysis.
* [MAJIQ](https://majiq.biociphers.org) and/or [LeafCutter](https://davidaknowles.github.io/leafcutter/) for local splicing variation detection.
	* This pipeline was created when LeafCutter documentation (and code) was on [commit 249fc26](https://github.com/davidaknowles/leafcutter/tree/249fc26f6e35201fc12fb560d347c9e792e64e5f). Since then documentation and potentially code changed considerably and *may* not be backward compatible despite keeping the same 0.2.9 version number...
* (optional) [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) should be installed on your system for read trimming.
* (optional) [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) should be installed on your system in a location included in your $PATH for quality control.
* (optional) [DESeqAnalysis](https://deseqanalysis.acidgenomics.com) for advanced graph options downstream of DESeq2 analysis.
* (optional) [cummeRbund](https://bioconductor.org/packages/release/bioc/html/cummeRbund.html) for advanced graph options downstream of Cufflinks. Not supported by the pipeline itself but results will be compatible with cummeRbund.

### Usage

```sh
sh_RNAseq.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini] [OutputDirName]
```


## sh_CRISPR.sh

This script will process the fastq(.gz) files generated in a typical CRISPR screen using either [casTLE](https://bitbucket.org/dmorgens/castle/) or [MAGeCK](https://sourceforge.net/projects/mageck/).
* If using casTLE, a reference file of all the indices will be automatically created using bowtie (NOT bowtie2). It will then analyze the screen and generate basic graphs.
* If using MAGeCK counts, tests, mle and pathway analysis will be performed. It will also run the [R](https://www.r-project.org) package "[MAGeCKFlute](https://bioconductor.org/packages/release/bioc/html/MAGeCKFlute.html)" and in all cases generate basic graphs.

See the [CRISPR.ini](https://github.com/emc2cube/Bioinformatics/blob/master/config_CRISPR.ini) configuration file for all available options and settings.

##### Options:
* ```--help``` : Display help message.
* ```--version``` : Display version number.

##### Dependencies:
* [casTLE](https://bitbucket.org/dmorgens/castle/) and/or [MAGeCK](https://sourceforge.net/projects/mageck/).
* [R](https://www.r-project.org) should be installed on your system in a location included in your $PATH.
* [MAGeCKFlute](https://bioconductor.org/packages/release/bioc/html/MAGeCKFlute.html) for downstream analysis if using MAGeCK.
* [csvkit](https://csvkit.readthedocs.io/en/latest/) should be installed on your system in a location included in your $PATH.
* [pathos](https://pypi.org/project/pathos/) should be installed on your system, this will provide [ppft](https://pypi.org/project/ppft/), a fork of [Parallel Python](https://www.parallelpython.com) working with both python2.7 and python3.6.

### Usage

```sh
sh_CRISPR.sh </path/to/fastq(.gz)/folder> </path/to/destination/folder> [/path/to/config/file.ini]
```

### Python 3.6 compatibility

For easy integration along MAGeCK, or any other modern tools, a python 3.6+ compatible version of casTLE is included.
This is based on [casTLE commit 981d6d8](https://bitbucket.org/dmorgens/castle/commits/981d6d877c0fe3ee233e9fd977b13800987a032c) and may not be up to date.
You still need to download the whole [casTLE repository](https://bitbucket.org/dmorgens/castle/) even if you end up switching the scripts with their python 3.6+ compatible version.
If you previously used casTLE with python 2.7 and already have [Parallel Python](https://www.parallelpython.com) installed you will need to uninstall it before installing [pathos](https://pypi.org/project/pathos/)


## sh_md5alldir.sh

This script will process all sub-directories of the input folders and for each of them will create a <directory_name>.md5 file if it does not exist yet, or check <directory> files against the existing <directory_name>.md5 file.

##### Options:
* ```-f``` or ```--force``` : even if there is already a <directory>.md5 file, it will be replaced by a new <directory>.md5 file.
* ```--help``` : Display help message.
* ```--version``` : Display version number.

### Usage

```sh
sh_md5alldir.sh [OPTIONS] </path/to/dir/>
```


## sh_sha1alldir.sh

This script will process all sub-directories of the input folders and for each of them will create a <directory_name>.sha1 file if it does not exist yet, or check <directory> files against the existing <directory_name>.sha1 file.

##### Options:
* ```-f``` or ```--force``` : even if there is already a <directory>.sha1 file, it will be replaced by a new <directory>.sha1 file.
* ```--help``` : Display help message.
* ```--version``` : Display version number.

### Usage

```sh
sh_sha1alldir.sh [OPTIONS] </path/to/dir/>
```


## sh_ACMGfilter.sh

This script will look for an annovar .snps.exome_summary.csv file and generate a list of all SNPs found in the ACMG guidelines in a new ACMG_genes.csv file.
This file can be directly sent to a clinician for incidental findings reports, if required.

##### Options:
* ```--help``` : Display help message.
* ```--version``` : Display version number.

### Usage

```sh
sh_ACMGfilter.sh </path/to/.csv/containing/folder> [/path/to/destination/folder]
```


## sh_mergeFastQ.sh

Simple script to consolidate fragmented .fastq files from different sequencing lanes.
Original files will be backed up in a FastQbackup folder.

##### Options:
* ```--help``` : Display help message.
* ```--version``` : Display version number.

### Usage

```sh
sh_mergeFastQ.sh </path/to/fastq(.gz)/folder>
```


## Author(s) contributions

👤 **Julien Couthouis**

*Initial work and releases*

* Linkedin: [@jcouthouis](https://www.linkedin.com/in/jcouthouis/)
* Github: [@emc2cube](https://github.com/emc2cube)

👤 **Rosa Ma**

*Local Splicing Variation*

* Linkedin: [@rosaxma](https://www.linkedin.com/in/rosaxma/)
* Github: [@rosaxma](https://github.com/rosaxma)


## Show your support

Give a ![GitHub stars](https://img.shields.io/github/stars/emc2cube/Bioinformatics?style=social) if this project helped you!


## License

Copyright © 2019 [Julien Couthouis](https://github.com/emc2cube).

This project is [EUPL-1.2](https://github.com/emc2cube/Bioinformatics/blob/master/LICENSE) licensed.