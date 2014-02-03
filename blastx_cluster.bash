#!/bin/bash

/sw/compbio/ncbi_blast+/2.2.28/bin/blastx \
-query $1 \
-db $2 \
-out $3 \
-evalue 1e-5 \
-outfmt 7 \
-query_gencode 11 \
-num_threads 2
