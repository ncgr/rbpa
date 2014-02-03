#!/bin/bash

if [[ $# < 4 ]]; then
	echo "useage:$0 alignCmd inputDir ref outDir"
	exit 1
fi
cmd=$1
ind=$2
ref=$3
outd=$4

cd $ind;

for f in *;do
	if [ ! -e "${outd}/$f" ];then
		qsub -N ctc_blastx \
		-pe smp 2 \
		-b n \
		-q normal.q \
		-o /home/align/logs/ \
		-e /home/align/logs/ \
		-j y \
		$cmd ${ind}/$f $ref ${outd}/$f
	fi
done
