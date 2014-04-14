#! /bin/bash

#placeholder for header EDIT

if [ $# -lt 3 ]; then

cat << _TEXT

	usage:$0 <kmer file> <"fastq file(s)"> <cores> [coverage]

	Use full paths.

	Run with kmer file, examples contained in kmer_lists.

	For multiple fastq files, or when using *, please use 
	"" around file list.

	Cores must be set with respect to your SGE cluster environment.	

_TEXT

exit 1
fi

kmer=$1
data=$2
cores=$3

if [ ! -z $4 ]; then
	cov=$4
else
	cov=5
fi

spcs=`pwd | awk -F/ '{print $NF}'`

mkdir -p abyss;cd abyss

if [ ! -f $kmer ]; then

cat << _TEXT

	ERROR: kmer file not found, "$kmer", please create kmer file.
	Examples can be found in the kmer_lists directory in BPA.

_TEXT

exit 1
fi

for f in $data; do
        if [ ! -f $f ]; then
                cat << _TEXT

        FAILED to locate fastq data in "$f" please check data input.

_TEXT

		exit 1
        fi
done

cp /dev/null ../properties.txt
echo "kmer file = $kmer" >> ../properties.txt
echo "cores = $cores" >> ../properties.txt
echo "input file(s) = `ls -1 $data | tr '\n' ' '`" >> ../properties.txt

path=$(dirname `which manage_cluster.bash`)
cp -f $kmer ./kmer.txt
kmers=`cat kmer.txt`

for i in $kmers; do
	let coretest=$(find $i -name "core*"|wc -l)
	let filetest=$(find $i -name "abyssrun-3.fa"|wc -l)
	label=$$'_'$i'_'$spcs
	echo $label
	if [[ "$coretest" -ne "0" || "$filetest" -ne "1" || ! -s $i/abyssrun-3.fa ]]; then
		if [ -d $i ]; then 
			rm -rf $i 2>/dev/null
		fi 
		mkdir -p $i
		( echo $i; cd $i; manage_cluster.bash $i "$data" $label $cores $cov $path > ${label}.out 2>${label}.err )
	fi
done

