#!/bin/bash
#
# Usage: sh_sha1alldir.sh </path/to/dir/> [-options, -? or --help for help]
#
## Description ##
#
# This script will process all sub-directories of the input folders and for each of them
# will create a <directory_name>.sha1 file if it does not exist yet, or check <directory> files
# against the existing <directory_name>.sha1 file.
#
## Options:
#
# -f or --force : even if a <directory>.sha1 file is detected, will replace it by a fresh one
# and will not check files against it.
#
##

dir="$1"

# Check paths and trailing / in directories
if [ -z "$dir" ]
then
	echo "Usage: sh_sha1alldir.sh </path/to/dir/> [-options, -? or --help for help]"
	exit
fi

if [ ${dir: -1} == "/" ]
then
	dir=${dir%?}
fi

echo ""

# Help!
if [ "$1" == "-?" ] || [ "$2" == "-?" ] || [ "$1" == "--help" ] || [ "$2" == "--help" ]
then
	echo "Usage: sh_sha1alldir.sh </path/to/dir/> [-options, -? or --help for help]"
	echo ""
	echo "Description"
	echo ""
	echo "This script will process all sub-directories of the input folders and for each of them"
	echo "will create a <directory_name>.sha1 file if it does not exist yet, or check <directory>"
	echo "files against the existing <directory_name>.sha1 file."
	echo ""
	echo "Options:"
	echo ""
	echo "-f or --force : even if a <directory_name>.sha1 file is detected, will overwrite it by a fresh one."
	echo "-? or --help : This help message..."
	echo ""
	exit
fi

folders=`find ${dir} -mindepth 1 -type d -not -path "*/.*"`

if [ -z "$folders" ]
then
	cd $dir
	if  [ "$2" == "--force" ] || [ "$2" == "-f" ]
	then
		files=`find -maxdepth 1 -type f -not -name ".*" -not -name "*.sha1" | sed "s#./##g"`
		if [ -n "${files}" ]
		then
			echo "Creating sha1 file for ${dir} contents"
			sha1sum $files > $(basename $dir).sha1
		fi
	elif [ -f "$(basename $dir).sha1" ]
	then
		echo "Checking $dir contents"
		sha1sum -c $(basename $dir).sha1
	else
		files=`find -maxdepth 1 -type f -not -name ".*" -not -name "*.sha1" | sed "s#./##g"`
		if [ -n "${files}" ]
		then
			echo "Creating sha1 file for ${dir} contents"
			sha1sum $files > $(basename $dir).sha1
		fi
	fi
else
	for i in $folders
	do
	cd $i
		if  [ "$2" == "--force" ] || [ "$2" == "-f" ]
		then
			files=`find -maxdepth 1 -type f -not -name ".*" -not -name "*.sha1" | sed "s#./##g"`
			if [ -n "${files}" ]
			then
				echo "Creating sha1 file for files in ${i}"
				sha1sum ${files} > $(basename ${i}).sha1
			fi
		elif [ -f "$(basename ${i}).sha1" ]
		then
			echo "Checking ${i} content"
			sha1sum -c $(basename ${i}).sha1
		else
			files=`find -maxdepth 1 -type f -not -name ".*" -not -name "*.sha1" | sed "s#./##g"`
			if [ -n "${files}" ]
			then
				echo "Creating sha1 file for files in $i"
				sha1sum basename ${files} > $(basename ${i}).sha1
			fi
		fi
	done
fi

echo ""
echo "Done"
