#! /bin/bash

if [ $# -lt 3 ];then

#usage:$0 <-r "read file(s)" | -b "bam file"> <-R reference.fa> <-k kmer> [-m maximum mismatch bwa] [-p min pairs] [-o output prefix] [-s seed length] [-t threads]

cat << _TEXT

        usage:$0 <-r "read file(s)"> <-R reference.fa> <-k kmer> [-m maximum mismatch bwa] [-p min pairs] [-o output prefix] [-s seed length] [-t threads]

_TEXT

  exit 1

fi

i=0

while getopts ":r:m:p:o:s:t:k:R:b:" fl; do
        case $fl in
                r)
                        reads=$OPTARG
                        for c in $reads;do
                                if [ ! -f $c ];then
                                        echo "could not locate file \"$c\", please provide a complete path"
                                        exit 1
                                else
                                        if [ $i -eq 0 ];then
                                                input=$c
                                                let i+=1
                                        else
                                                input="$input $c"
                                        fi
                                fi
                        done
                        ;;
		b)
			bam=$OPTARG
			;;
                m)
                        mmatch=$OPTARG
                        ;;
                p)
                        mpairs=$OPTARG
                        ;;
		o)
			pref=$OPTARG
			;;
		s)
			seed=$OPTARG
			;;
		t)
			threads=$OPTARG
			;;
		k)
			kmer=$OPTARG
			;;
		R)
			ref=$OPTARG
			;;
                \?)
                        echo "$fl is not a valid option"
                        exit 1
                        ;;
                :)
                        echo "$fl requires an argument"
                        exit 1
                        ;;
        esac
done

mkdir -p scaffold;cd scaffold

mkdir -p aln ref;cd ref

if [ ! -f $ref ]; then
	echo "could not locate $ref, please provide a complete path"
	exit 1
fi

if [ -z $mmatch ]; then
	mmatch=0.04
fi

if [ -z $mpairs ]; then
	mpairs=5
fi

if [ -z $pref ]; then
	pref="bpa_out"
fi

if [ -z $seed ]; then
	seed=200
fi

if [ -z $threads ]; then
	threads=2
fi

#if [ -z $bam ];then

awk '{if(/^>/){a++;print ">"a}else{print}}' $ref > scaffolds_in.fna

bwa index scaffolds_in.fna

cd ../aln

bwa aln -t $threads -n $mmatch -f $pref.sai ../ref/scaffolds_in.fna <( cat $input )
bwa samse -f $pref.int ../ref/scaffolds_in.fna $pref.sai <( cat $input )
cd ..

cat aln/$pref.int | abyss-fixmate -h scaffold.hist \
|sort -T /tmp/ -snk3 -k4 \
|DistanceEst --dot -v -k$kmer -j$threads -s$seed -n$mpairs -o scaffold.dist.dot scaffold.hist
abyss-todot -v -e ref/scaffolds_in.fna scaffold.dist.dot > gmerge.dist.dot
abyss-scaffold -v -n$mpairs -s$seed -k$kmer -o scaffold.paths -g scaffold.adj gmerge.dist.dot
MergeContigs -v -k$kmer -o scaffolds.fna ref/scaffolds_in.fna gmerge.dist.dot scaffold.paths
#else
#	samtools view -h $bam | abyss-fixmate -h scaffold.hist \
#         |sort -T /tmp/ -snk3 -k4 \
#         |DistanceEst --dot -v -k$kmer -j$threads -s$seed -n$mpairs -o scaffold.dist.dot scaffold.hist
#         abyss-todot -v -e ref/scaffolds_in.fna scaffold.dist.dot > gmerge.dist.dot
#         abyss-scaffold -v -n$mpairs -s$seed -k$kmer -o scaffold.paths -g scaffold.adj gmerge.dist.dot
#         MergeContigs -v -k$kmer -o scaffolds.fna ref/scaffolds_in.fna gmerge.dist.dot scaffold.paths
#fi
