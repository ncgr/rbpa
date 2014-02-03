#!/bin/bash

if [ $# -lt 2 ];then
	echo "usage: $0 <reference> <fa>"
	exit 1
fi

. /home/analysis/ctc/dev/abyss_1.3.3/paths.bash

scripath=`dirname $0`
reference=$1

mkdir -p uniref; cd uniref;mkdir -p tmp

path=$PWD

cat $2 | $scripath/seqs_chunkify -n 100

$scripath/run-blast.bash $scripath/blastx.bash $path/chunks $reference $path/tmp
