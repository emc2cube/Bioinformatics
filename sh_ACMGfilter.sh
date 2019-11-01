#!/bin/bash
#
# Usage: sh_ACMGfilter.sh </path/to/.csv/containing/folder> [/path/to/destination/folder]
#
##############################################################
##                      Description                         ##
##############################################################
#
# This script will look for an annovar .snps.exome_summary.csv file and generate a list
# of all SNPs found in the ACMG guidelines in a new ACMG_genes.csv file. This file can be directly
# sent to a clinician for incidental findings reports, if required.
#
##

# Help!
if [ "${1}" == "--help" ] || [ "${2}" == "--help" ] || [ "${3}" == "--help" ]
then
	echo "Usage: $(basename $0) </path/to/.csv/containing/folder> [/path/to/destination/folder]"
	echo ""
	echo "Description"
	echo ""
	echo "This script will look for an annovar .snps.exome_summary.csv file and generate a list"
	echo "of all SNPs found in the ACMG guidelines in a new ACMG_genes.csv file. This file can be directly"
	echo "sent to a clinician for incidental findings reports, if required."
	echo ""
	echo "Options:"
	echo "$(basename $0) --help : Display this help message."
	echo ""
	exit
fi

# Version
if [ "${1}" == "--version" ] || [ "${2}" == "--version" ] || [ "${3}" == "--version" ]
then
	echo "$(basename $0) version 2.0"
	echo "ACMG SF v2.0"
	exit
fi


dir="${1}"
out="${2}"

# Check paths and trailing / in directories
if [ -z "$dir" ]
then
	$(echo "${0} --help")
	exit
fi

if [ ${dir: -1} == "/" ]
then
	dir=${dir%?}
fi

if [ -z ${out} ]
then
	out="."
fi

if [ ${out: -1} == "/" ]
then
	out=${out%?}
fi

[ -f ${out}/filecontent ] && rm ${out}/filecontent

annovarfile=$(ls ${dir}/ | grep _multianno.csv)
fixheaderfile=$(ls ${dir}/ --hide=*.idx | grep .filtered.vcf)
headerfix=$(cat ${dir}/$fixheaderfile | grep CHROM | sed 's/#//g' | sed 's/\t/,/g')

echo "`head -1 ${dir}/${annovarfile} | sed 1s/Otherinfo//g`$headerfix" > ${out}/ACMG_genes.csv

for i in BRCA1 BRCA2 TP53 STK11 MLH1 MSH2 MSH6 PMS2 APC MUTYH VHL MEN1 RET NTRK1 PTEN RB1 SDHD SDHAF2 SDHC SDHB TSC1 TSC2 WT1 NF2 COL3A1 FBN1 TGFBR1 TGFBR2 SMAD3 ACTA2 MYH11 MYBPC3 MYH7 TNNT2 TNNI3 TPM1 MYL3 ACTC1 PRKAG2 GLA MYL2 LMNA RYR2 PKP2 DSP DSC2 TMEM43 DSG2 KCNQ1 KCNH2 SCN5A LDLR APOB PCSK9 RYR1 CACNA1S BMPR1A SMAD4 ATP7B OTC
do
	tail -n +2 ${dir}/${annovarfile} | grep \",${i},\" >> ${out}/ACMG_genes.csv
done

echo "ACMG guidelines gene list generated as ${out}/ACMG_genes.csv"
