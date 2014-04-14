#! /bin/bash

###############################################################################
#                               run-bwa-single.bash
#                       Shell for realigning single-end reads
#                               February 4, 2014
#                               NCGR: www.ncgr.org
###############################################################################

. /home/analysis/ctc/rbpa/paths.bash

if [ $# -lt 3 ];then

cat << USAGE

	usage:$0 <in_fastq> <ref> <out_prefix> <threads> [mismatch]

USAGE

exit 1
fi

input=$1
ref=$2
output_file_prefix=$3
threads=$4
mm=0.04

if [ ! -z $5 ]; then
	mm=$5
fi

bwa aln -n $mm -t $threads -f ${output_file_prefix}.0.sai $ref <( cat $input )

bwa samse -f ${output_file_prefix}.0.int $ref ${output_file_prefix}.0.sai <( cat $input )

samtools view -bT $ref ${output_file_prefix}.0.int | samtools sort - ${output_file_prefix}.0

