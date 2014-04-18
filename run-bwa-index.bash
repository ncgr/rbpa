#! /bin/bash

###############################################################################
#                               run-bwa-index.bash
#       		Shell for launching bwa index
#                               February 4, 2014
#                               NCGR: www.ncgr.org
###############################################################################

. /home/analysis/$USER/.bash_profile

if [ $# -ne 1 ];then

	cat << _TEXT
	
	usage:$0 <ref>

_TEXT
exit 1
fi

ref=$1

if [ ! -f $ref ];then
	echo "could not find $ref.  please provide a complete path"
	exit 1
fi

bwa index $ref
