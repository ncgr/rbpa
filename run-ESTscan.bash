#!/bin/bash


#############command line args check#####################

if [ $# -lt 2 ]; then
	echo "usage:$0 <in> <name.smat> <out> [filter_length]"
	exit 1
fi

in=$1
smat=$2
out=$3
len=30


if [ ! -z $4 ]; then
	len=$4
fi

if [ -d estscan ]; then
	echo "estscan directory exists.  please remove before running"
	exit 1
fi

mkdir estscan

#############Begin peptide prediction###################

	estscan -M $smat -l $len -t estscan/${out}_protein_motifs_raw.faa $in > estscan/${out}_motif_seqs.fna 

echo "Parameter File=$smat" > estscan/smat.txt
