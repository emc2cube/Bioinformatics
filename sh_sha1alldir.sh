#!/bin/bash
#
# Usage: sh_sha1alldir.sh </path/to/dir/> [OPTIONS]
#
##############################################################
##                       Description                        ##
##############################################################
#
# This script will process all sub-directories of the input folders and for each of them
# will create a <directory_name>.sha1 file if it does not exist yet, or check <directory>
# files against the existing <directory_name>.sha1 file.
#
## Options:
#
# -f or --force : even if there is already a <directory>.sha1 file, it will be replaced by a new <directory>.sha1 file.
# --help : Display help message.
# --version : Display version number.
#
##

# Help!
if [ "${1}" == "--help" ] || [ "${2}" == "--help" ] || [ "${3}" == "--help" ]
then
	echo "Usage: $(basename $0) </path/to/dir/> [OPTIONS]"
	echo ""
	echo "Description"
	echo ""
	echo "This script will process all sub-directories of the input folders and for each of them"
	echo "will create a <directory_name>.sha1 file if it does not exist yet, or check <directory>"
	echo "files against the existing <directory_name>.sha1 file."
	echo ""
	echo "Options:"
	echo "$(basename $0) -f or --force : even if there is already a <directory>.sha1 file, it will be replaced by a new <directory>.sha1 file."
	echo "$(basename $0) --help : Display this help message."
	echo "$(basename $0) --version : Display version number."
	echo ""
	exit
fi

# Version
if [ "${1}" == "--version" ] || [ "${2}" == "--version" ] || [ "${3}" == "--version" ]
then
	echo "$(basename $0) version 1.0"
	exit
fi

dir="${1}"

# Check paths and trailing / in directories
if [ -z "${dir}" ]
then
	$(echo "${0} --help")
	exit
fi

if [ ${dir: -1} == "/" ]
then
	dir=${dir%?}
fi

echo ""

folders=$(find ${dir} -mindepth 1 -type d -not -path "*/.*")

if [ -z "${folders}" ]
then
	cd ${dir}
	if  [ "${2}" == "--force" ] || [ "${2}" == "-f" ]
	then
		files=$(find -maxdepth 1 -type f -not -name ".*" -not -name "*.sha1" | sed "s#./##g")
		if [ -n "${files}" ]
		then
			echo "Creating sha1 file for ${dir} contents"
			sha1sum ${files} > $(basename ${dir}).sha1
		fi
	elif [ -f "$(basename ${dir}).sha1" ]
	then
		echo "Checking ${dir} contents"
		sha1sum -c $(basename ${dir}).sha1
	else
		files=$(find -maxdepth 1 -type f -not -name ".*" -not -name "*.sha1" | sed "s#./##g")
		if [ -n "${files}" ]
		then
			echo "Creating sha1 file for ${dir} contents"
			sha1sum ${files} > $(basename ${dir}).sha1
		fi
	fi
else
	for i in ${folders}
	do
		cd ${i}
		if  [ "${2}" == "--force" ] || [ "${2}" == "-f" ]
		then
			files=$(find -maxdepth 1 -type f -not -name ".*" -not -name "*.sha1" | sed "s#./##g")
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
			files=$(find -maxdepth 1 -type f -not -name ".*" -not -name "*.sha1" | sed "s#./##g")
			if [ -n "${files}" ]
			then
				echo "Creating sha1 file for files in ${i}"
				sha1sum basename ${files} > $(basename ${i}).sha1
			fi
		fi
	done
fi

echo ""
echo "Done"
