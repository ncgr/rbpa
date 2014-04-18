#! /bin/bash
#$ -S /bin/bash
#$ -m e
#$ -cwd

###############################################################################
#                               run-abyss-post.bash
#       		Shell for completing unipaths
#                               February 4, 2014
#                               NCGR: www.ncgr.org
###############################################################################

if [ z$1 = z ]; then
  exit 1
fi

kmer=$1
path=$2

export PATH=$path:$PATH

AdjList -k${kmer} -m30 abyssrun-1.fa >abyssrun-1.adj
abyss-filtergraph  -k${kmer} -g abyssrun-2.adj abyssrun-1.adj >abyssrun-1.path
PopBubbles -j8 -v -k${kmer} -p0.9 -g abyssrun-3.adj abyssrun-1.fa abyssrun-2.adj > abyssrun-2.path
MergeContigs -k${kmer} -o abyssrun-3.fa abyssrun-1.fa abyssrun-2.adj abyssrun-2.path
awk '!/^>/ {x[">" $1]=1; next} {getline s} $1 in x {print $0 "\n" s}' \
                abyssrun-2.path abyssrun-1.fa >abyssrun-indel.fa
ln -sf abyssrun-3.fa abyssrun-unitigs.fa

