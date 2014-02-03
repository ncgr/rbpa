#! /bin/bash
#$ -S /bin/bash
#$ -m e
#$ -cwd

if [ z$1 = z ]; then
  exit 1
fi

kmer=$1
data=$2
cov=$3
path=$4
mpi=$5
lib=$6
openmpihome="$lib"
#abysshome=$ABYSS

export LD_LIBRARY_PATH=${openmpihome}/lib:${LD_LIBRARY_PATH}

export OMPI_MCA_plm_rsh_disable_qrsh=1
export OMPI_MCA_btl_tcp_if_include=eth0

echo "$lib" > lib.txt

cat $PE_HOSTFILE | cut -f1 -d' ' > hosts.dat

$mpi/mpirun \
 --prefix ${openmpihome} \
 -v \
 -np $NSLOTS \
 --hostfile hosts.dat \
 $path/ABYSS-P \
 --out=abyssrun-1.fa \
 --kmer=${kmer} \
 --graph=se-graph.dot \
 --snp=popped-bubbles.fa \
 --coverage-hist=k-mer_coverage.hist \
 --coverage=$cov \
 --verbose \
 $data \

