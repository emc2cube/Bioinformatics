#!/bin/bash

#
# Goes into each subdirectory and creates a <directoryname>-checksum.md5 file
# with checksums for all files within (incl subdir or subdir)
#
# Steve Allison -- https://www.nooblet.org/
#

# Temporary directory
TMPDIR="/dev/shm/checksum-create"
# Debug mode
DEBUG=0

# Grab debug flag -d
if [ "$1" = "-d" ]; then
	DEBUG=1;
	shift;
fi

# Give group write access by default
umask 022

# Will be set to 1 if  checksum file is created/modified
EDIT=0

# Check we have directory
if [ "$1" = "" ]; then
	echo "Syntax: $0 <dir>"; exit 1;
fi

# Check we have md5sum
if ! which md5sum &>/dev/null; then
	echo "Error: md5sum not found"
	exit 1;
fi

# Only outputs if DEBUG=1, also stores output for later use
debug() {
	echo "$*" >> ${TMPDIR}/output
	if [ $DEBUG -eq 1 ]; then
		echo "$*"
	fi
}

# Header of checksum file
_header() {
	echo "# Checksum created on $(date -d "@${1}")"
	echo "# by $USER@$(hostname --fqdn) ($0)"
	echo "# for $PWD"
	echo
}

# Checksum file has been edited
_edited() {
	if [ ! -e "${TMPDIR}/edited" ]; then touch ${TMPDIR}/edited; fi
	EDIT=1
}

# Parse data line from existing checksum file
_data() {
	dir="${PWD}"
	# ; data creation modified size md5 filename
	if [ "$7" = "" ]; then
		_edited
		debug "   - nodata: ${file}"
		return
	fi
	shift; shift
	# Parse data line
	creation="$1"; shift
	modified="$1"; shift
	size="$1"; shift
	md5="$1"; shift
	file="$*"
	# If file no longer exists, skip
	if [ ! -e "${file}" ]; then
		_edited
		debug "   - missing: ${file}"
		return
	# If date is different, compute checksum
	elif [ "$(date -r "${file}" +%s)" -ne "${modified}" ]; then
		debug "   - mtime: ${file} # $(date -r "${file}" +%s) # ${modified}"
		_addfile ${creation} "${file}"
	# If size is different, compute checksum
	elif [ "$(stat -c%s "${file}")" -ne "${size}" ]; then
		debug "   - size: ${file} # $(stat -c%s "${file}") # $size"
		_addfile ${creation} "${file}"
	# Else just add the same line again
	else
		echo "# data $creation $modified $size $md5 $file" >> ${TMPDIR}/md5file
		echo "${md5} *${file}" >> ${TMPDIR}/md5file
	fi
}

# Use md5sum to create checksum and add data line
_addfile() {
	dir="${PWD}"
	# We have now edited the checksum file
	_edited
	# Creation date given on commandline
	creation=$1; shift
	file=$*
	modified=$(date -r "${file}" +%s)
	if [ "$creation" -eq 0 ]; then
		creation=${modified}
	fi
	# This is where the magic happens
	md5="$(md5sum -b $file|awk '{print $1}')"
	# Get filesize
	size="$(stat -c%s "${file}")"
	# Output data to checksum file
	echo "# data $creation $modified $size $md5 $file" >> ${TMPDIR}/md5file
	echo "${md5} *${file}" >> ${TMPDIR}/md5file
}

# Clean workspace
rm -rf ${TMPDIR}
mkdir -p ${TMPDIR}

# Lets begin
debug " ! Looping $1 .."

# Get list of directories
FILES=$(find $1 -maxdepth 1 -mindepth 1 -type d)
while read dir; do
	debug " o ${dir}"
	begin="$(date +%s)"
	# Current checksum file
	md5file="$(basename ${dir})-checksums.md5"
	cd "${dir}"
	# If checksum file exists, we need to parse data lines
	if [ -e "${md5file}" ]; then
		OUTPUT=$(grep "^# data " ${md5file} 2>&1)
	while read data; do if [ "$data" != "" ]; then _data $data; fi; done << EOF
${OUTPUT}
EOF
	fi
	# Now we find any new files
	OUTPUT=$(find -type f|sed -e 's/^\.\///' -e '/\-checksums\.md5$/d')
	while read file; do
		if [ "$file" = "" ]; then continue; fi
		# Don't want to add checksum.md5 to the md5file
		if [[ "$file" =~ -checksums.md5$ ]]; then continue; fi
		# Does the file already exist in the md5file?
		if ! grep "^# data [0-9]\+ [0-9]\+ [0-9]\+ [0-9a-f]\+ $file\$" ${TMPDIR}/md5file &>/dev/null; then
			# Doesn't exist, so compute the checksum
			debug "   - newfile: ${file}"
			_addfile 0 "$file"
		fi
	done << EOF
${OUTPUT}
EOF
	# If checksum file has been edited, copy data
	if [ -e ${TMPDIR}/edited ]; then
		_header $begin > ${md5file}
		cat ${TMPDIR}/md5file >> ${md5file}
		echo -e "\n# Time taken: $(($(date +%s) - $begin)) secs" >> ${md5file}
	fi
	cd ..
	# Clean working directory
	rm -rf ${TMPDIR}/md5file
	rm -rf ${TMPDIR}/edited
done << EOF
${FILES}
EOF
# If debug is disabled, and the files were edited, output the changed
if [ "$DEBUG" -eq 0 ] && [ "$EDIT" -eq 1 ]; then
	cat ${TMPDIR}/output
fi
# Purge working directory
rm -rf ${TMPDIR}
