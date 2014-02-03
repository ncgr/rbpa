#!/bin/bash

###############################################################################
#				manage_cluster.bash			
#	Shell for launching openMPI enabled ABySS and completing unipaths
# 				September 5, 2013
###############################################################################

kmer=$1
data=$2
label=$3
cores=$4
cov=$5
script=$6
echo $script
path=$(dirname `which ABYSS-P`)
mpi=$(dirname `which mpirun`)
lib=$(dirname `which mpirun | sed -e 's/bin//'`)

#run the unitigs stage
name1=abyss_unitigs_${label}_$$
qsub -pe orte $cores -N ${name1} $script/run-abyss-unitigs.bash $kmer "$data" $cov $path $mpi $lib

#the -hold_jid may reference multiple job patterns, according to SGE type definitions.
name2=abyss_post_${label}_$$
qsub -N ${name2} -hold_jid ${name1} $script/run-abyss-post.bash $kmer $path

