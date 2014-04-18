#!/bin/bash -e

###############################################################################
#                               run-gapclose.bash
#       		Shell for launching GapCloser
#                               February 4, 2014
#                               NCGR: www.ncgr.org
###############################################################################

if [ $# -lt 5 ]; then
        echo "usage:$0 <r1.fq> <r2.fq> <ref.fa> <max_length> <threads> [insert_size]"
        exit 1
fi

r1=$1
r2=$2
ref=$3
len=$4
threads=$5
ins=$6

if [[ -f scaffold/scaffold.hist && -s scaffold/scaffold.hist ]]; then
	ins=`cat scaffold/scaffold.hist | awk '{if($2 >= b){a=$1;b=$2}}END{print a}'`
fi

if [ -z $ins ]; then
	echo "if not using pipeline scaffolder, insert size must be provided"
	exit 1
fi

#jname=n_$$_gapclose_${species}
#echo ${jname}

mkdir -p gapclose;cd gapclose

#################creates gapclose config file#######################

echo -e "#maximal read length\n""max_rd_len=$len\n""#2x$len\n""[LIB]\n""avg_ins=$ins\n""reverse_seq=0\n""asm_flags=3\n""rank=1" > gapclose.conf.header

if [ `echo $r1 | tr ' ' '\n' | wc -l` -eq `echo $r2 | tr ' ' '\n' | wc -l` ]; then
	paste <( echo $r1 | tr ' ' '\n' | sort ) <( echo $r2 | tr ' ' '\n' | sort ) | while read l; do
		i=0
		for r in `echo $l`; do
			if [ ! -f $r ]; then
				echo "$r does not exist, please check path"
				exit 1
			fi
			if [ $i -eq 0 ]; then
				echo "q1=$r" >> gapclose.conf.header
				let i+=1
			else
				echo "q2=$r" >> gapclose.conf.header
			fi
		done
	done
else
	echo "mate 1 and mate 2 did not contain the same number of files"
	exit 1
fi

#################end of gapclose.conf.header#######################

cat $ref |  awk '{if(/^>/){a+=1;if(a==1){print ">BPA_"a}else{printf "\n>BPA_"a"\n"}}else{printf $0}}' > gapclose_in.fna

#################lauch jobs########################################

GapCloser -a gapclose_in.fna -b gapclose.conf.header -t $threads -o gapclosed.fna > log_gapclose.out

##################END#####################
