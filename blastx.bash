#!/bin/bash

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

/sw/compbio/ncbi_blast+/2.2.28/bin/blastx \
-query $query \
-db $ref \
-out $out \
-evalue $e \
-outfmt $fmt \
-query_gencode $code \
-num_threads $threads
