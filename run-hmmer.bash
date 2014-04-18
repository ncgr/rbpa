#!/bin/bash

###############################################################################
#                               run-hmmer.bash
#       		Shell for launching hmmscan
#                               February 4, 2014
#                               NCGR: www.ncgr.org
###############################################################################

if [ $# -ne 4 ]; then
	echo "usage:$0 <in> <db> <out> <cpu>"
	exit 1
fi

input=$1
db=$2
out=$3
cpu=$4

#mkdir -p hmmscan;cd hmmscan
						      #dir label    hmm
hmmscan -o $out.hits --tblout $out.hits.tbl --cpu $cpu --notextw $db $input 2> hmm3_$out.err
