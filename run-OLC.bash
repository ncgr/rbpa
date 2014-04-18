#!/bin/bash

###############################################################################
#                               run-OLC.bash
#               Shell for launching CD-HIT-EST and CAP3 OLC
#                               February 4, 2014
#                               NCGR: www.ncgr.org
###############################################################################

if [ $# -lt 1 ];then
	
cat << _TEXT

        usage:$0 <-f "fasta file(s)"> [-c "cap3 args"] [-d "cd-hit-est args"] 

        For multiple fasta files, or when using *, please use
        "" around file list.  If not using the perl wrapper, please
	provide all arguments in quotations.

_TEXT

exit 1
fi

fcheck=0
dcheck=0
i=0

while getopts ":f:c:d:" fl; do
	case $fl in
		f)
			files=$OPTARG
			for c in $files;do
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

			fcheck=1
			;;
		c)
			c3args=$OPTARG
			;;
		d)
			ceargs=$OPTARG
			dcheck=1
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

if [ $fcheck -eq 0 ];then
	echo "-f required, see usage"
	exit 1
fi

mkdir -p overlap
cd overlap

if [ -f pool_uncondensed.fna ];then rm pool_uncondensed.fna pool_final.fna;fi

cat $input | awk '{if(/^>/){a++;gsub(">","");print ">"a" "$0}else{print}}' > pool_uncondensed.fna

if [ $dcheck -eq 0 ];then
	mv -f pool_uncondensed.fna pool_final.fna
else	
	echo "$ceargs"
	cd-hit-est -i pool_uncondensed.fna -o pool_final.fna -M 0 $ceargs
fi
echo $c3args
cap3 pool_final.fna $c3args
cat pool_final.fna.cap.contigs pool_final.fna.cap.singlets > cap3_final.fna

