#!/bin/bash
#
# Usage: sh_sha1alldir.sh [OPTIONS] </path/to/dir/>
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
if [ "${1}" = "--help" ]
then
	echo "Usage: $(basename "$0") [OPTIONS] </path/to/dir/>"
	echo ""
	echo "Description"
	echo ""
	echo "This script will process all sub-directories of the input folders and for each of them"
	echo "will create a <directory_name>.sha1 file if it does not exist yet, or check <directory>"
	echo "files against the existing <directory_name>.sha1 file."
	echo ""
	echo "Options:"
	echo "$(basename "$0") -f or --force : even if there is already a <directory>.sha1 file, it will be replaced by a new <directory>.sha1 file."
	echo "$(basename "$0") --help : Display this help message."
	echo "$(basename "$0") --version : Display version number."
	echo ""
	exit
fi

# Version
if [ "${1}" = "--version" ]
then
	echo "$(basename "$0") version 1.0.1"
	exit
fi

if  [ "${1}" = "--force" ] || [ "${1}" = "-f" ]
then
	dir="${2}"
else
	dir="${1}"
fi

# Check paths and trailing / in directories
if [ -z "${dir}" ]
then
	${0} --help
	exit
fi

echo ""

folders=$(find -L "${dir}" -mindepth 1 -type d -not -path "*/.*")

if [ -z "${folders}" ]
then
	cd "${dir}" || exit
	if  [ "${1}" = "--force" ] || [ "${1}" = "-f" ]
	then
		files=$(find -L . -maxdepth 1 -type f -not -name ".*" -not -name "*.sha1" -not -name "*.md5" | sed "s#./##g" | sort -n)
		if [ -n "${files}" ]
		then
			echo "Creating sha1 file for ${dir} contents"
			sha1sum ${files} > "$(basename "${dir}")".sha1
		fi
	elif [ -f "$(basename "${dir}").sha1" ]
	then
		echo "Checking ${dir} contents"
		sha1sum -c "$(basename "${dir}")".sha1
	else
		files=$(find -L . -maxdepth 1 -type f -not -name ".*" -not -name "*.sha1" -not -name "*.md5" | sed "s#./##g" | sort -n)
		if [ -n "${files}" ]
		then
			echo "Creating sha1 file for ${dir} contents"
			sha1sum ${files} > "$(basename "${dir}")".sha1
		fi
	fi
else
	for i in ${folders}
	do
		cd "${i}" || exit
		if  [ "${1}" = "--force" ] || [ "${1}" = "-f" ]
		then
			files=$(find -L . -maxdepth 1 -type f -not -name ".*" -not -name "*.sha1" -not -name "*.md5" | sed "s#./##g" | sort -n)
			if [ -n "${files}" ]
			then
				echo "Creating sha1 file for files in ${i}"
				sha1sum ${files} > "$(basename "${i}")".sha1
			fi
		elif [ -f "$(basename "${i}").sha1" ]
		then
			echo "Checking ${i} content"
			sha1sum -c "$(basename "${i}")".sha1
		else
			files=$(find -L . -maxdepth 1 -type f -not -name ".*" -not -name "*.sha1" -not -name "*.md5" | sed "s#./##g" | sort -n)
			if [ -n "${files}" ]
			then
				echo "Creating sha1 file for files in ${i}"
				sha1sum basename ${files} > "$(basename "${i}")".sha1
			fi
		fi
	done
fi

echo ""
echo "Done"
