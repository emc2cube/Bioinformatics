#!/bin/bash
#
# Usage: sh_md5alldir.sh [OPTIONS] </path/to/dir/>
#
##############################################################
##                       Description                        ##
##############################################################
#
# This script will process all sub-directories of the input folders and for each of them
# will create a <directory_name>.md5 file if it does not exist yet, or check <directory>
# files against the existing <directory_name>.md5 file.
#
## Options:
#
# -f or --force : even if there is already a <directory>.md5 file, it will be replaced by a new <directory>.md5 file.
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
	echo "will create a <directory_name>.md5 file if it does not exist yet, or check <directory>"
	echo "files against the existing <directory_name>.md5 file."
	echo ""
	echo "Options:"
	echo "$(basename "$0") -f or --force : even if there is already a <directory>.md5 file, it will be replaced by a new <directory>.md5 file."
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
		files=$(find -L . -maxdepth 1 -type f -not -name ".*" -not -name "*.md5" -not -name "*.md5" | sed "s#./##g" | sort -n)
		if [ -n "${files}" ]
		then
			echo "Creating md5 file for ${dir} contents"
			md5sum ${files} > "$(basename "${dir}")".md5
		fi
	elif [ -f "$(basename "${dir}").md5" ]
	then
		echo "Checking ${dir} contents"
		md5sum -c "$(basename "${dir}")".md5
	else
		files=$(find -L . -maxdepth 1 -type f -not -name ".*" -not -name "*.md5" -not -name "*.md5" | sed "s#./##g" | sort -n)
		if [ -n "${files}" ]
		then
			echo "Creating md5 file for ${dir} contents"
			md5sum ${files} > "$(basename "${dir}")".md5
		fi
	fi
else
	for i in ${folders}
	do
		cd "${i}" || exit
		if  [ "${1}" = "--force" ] || [ "${1}" = "-f" ]
		then
			files=$(find -L . -maxdepth 1 -type f -not -name ".*" -not -name "*.md5" -not -name "*.md5" | sed "s#./##g" | sort -n)
			if [ -n "${files}" ]
			then
				echo "Creating md5 file for files in ${i}"
				md5sum ${files} > "$(basename "${i}")".md5
			fi
		elif [ -f "$(basename "${i}").md5" ]
		then
			echo "Checking ${i} content"
			md5sum -c "$(basename "${i}")".md5
		else
			files=$(find -L . -maxdepth 1 -type f -not -name ".*" -not -name "*.md5" -not -name "*.md5" | sed "s#./##g" | sort -n)
			if [ -n "${files}" ]
			then
				echo "Creating md5 file for files in ${i}"
				md5sum basename ${files} > "$(basename "${i}")".md5
			fi
		fi
	done
fi

echo ""
echo "Done"
