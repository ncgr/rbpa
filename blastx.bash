#!/bin/bash

###############################################################################
#                               blastx.bash
#       		Shell for launching blastx
#                               February 4, 2014
#                               NCGR: www.ncgr.org
###############################################################################

if [ $# -lt 7 ];then
	echo "usage:$0 <query> <ref> <out> <threads> <eval> <outfmt> <gencode>"
	exit 1
fi

query=$1
ref=$2
out=$3
threads=$4
e=$5
fmt=$6
code=$7

blastx \
-query $query \
-db $ref \
-out $out \
-evalue $e \
-outfmt $fmt \
-query_gencode $code \
-num_threads $threads
