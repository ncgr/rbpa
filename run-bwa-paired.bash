#! /bin/bash

. /home/analysis/$USER/.bash_profile

if [ $# -lt 4 ];then

	cat << _TEXT
	
	usage:$0 <r1> <r2> <ref> <outprefix> <threads> [miss match]

_TEXT
exit 1
fi

input1=$1
input2=$2
ref=$3
output_file_prefix=$4
threads=$5
mism=0.04

if [ ! -z $6 ];then
	mism=$6
fi


bwa aln -n $mism -t $threads -f ${output_file_prefix}.1.sai $ref <( cat $input1 )

bwa aln -n $mism -t $threads -f ${output_file_prefix}.2.sai $ref <( cat $input2 )

bwa sampe -f ${output_file_prefix}.int $ref ${output_file_prefix}.1.sai ${output_file_prefix}.2.sai <( cat $input1 ) <( cat $input2 )

samtools view -bT $ref ${output_file_prefix}.int | samtools sort - ${output_file_prefix}
