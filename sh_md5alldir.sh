#!/bin/bash
#
# Usage: sh_md5alldir.sh </path/to/dir/> [OPTIONS]
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
if [ "${1}" == "--help" ] || [ "${2}" == "--help" ] || [ "${3}" == "--help" ]
then
	echo "Usage: $(basename $0) </path/to/dir/> [OPTIONS]"
	echo ""
	echo "Description"
	echo ""
	echo "This script will process all sub-directories of the input folders and for each of them"
	echo "will create a <directory_name>.md5 file if it does not exist yet, or check <directory>"
	echo "files against the existing <directory_name>.md5 file."
	echo ""
	echo "Options:"
	echo "$(basename $0) -f or --force : even if there is already a <directory>.md5 file, it will be replaced by a new <directory>.md5 file."
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

folders=$(find $(readlink -f ${dir}) -mindepth 1 -type d -not -path "*/.*")

if [ -z "${folders}" ]
then
	if [ -f "${dir}/$(basename $dir).md5" ] && [ "$2" != "--force" -o "$2" != "-f" ]
	then
		echo "Checking ${dir} content"
		cd ${dir}
		md5sum -c $(basename ${dir}).md5 2>/dev/null
	else
		echo "Creating md5 file for ${dir} content"
		md5sum ${dir}/* > ${dir}/$(basename $dir).md5 2>/dev/null
	fi
else
	for i in ${folders}
	do
		if [ -f "${i}/$(basename $i).md5" ] && [ "${2}" != "--force" -o "${2}" != "-f" ]
		then
			echo "Checking ${i} content"
			cd ${i}
			md5sum -c $(basename ${i}).md5 2>/dev/null
		else
			echo "Creating md5 file for files in ${i}"
			files=$(find $(readlink -f ${i}) -type f -not -name ".*")
			md5sum ${files} > ${i}/$(basename $i).md5 2>/dev/null
		fi
	done
fi

echo ""
echo "Done"
