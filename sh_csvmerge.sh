#!/bin/bash

dir="$1"
out="$2"

# Check paths and trailing / in directories
if [ -z $dir ]
then
    echo "usage: sh_csvmerge.sh <.csv folder> [destination folder]"
    exit
fi

if [ ${dir: -1} == "/" ]
then
    dir=${dir%?}
fi

if [ -z $out ]
then
    out="."
fi

if [ ${out: -1} == "/" ]
then
    out=${out%?}
fi

[ -f $out/filecontent ] && rm $out/filecontent

csvfiles=`ls $dir/ | grep .snps.exome_summary.csv`
header=`ls $dir | grep .snps.exome_summary.csv | head -1`

echo "Sample,`head -1 $dir/$header`" > $out/All_SNPs_merged.csv

for i in $csvfiles
do
    filename=`echo $i | awk -F. '{print $1}'` 
    tail -n +2 $dir/$i > $out/filecontent
    while read line
    do 
        echo "\"$filename\",$line" >> $out/All_SNPs_merged.csv
        done < $out/filecontent
done

rm $out/filecontent
