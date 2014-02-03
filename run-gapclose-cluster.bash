#! /bin/bash

. /home/analysis/ctc/dev/abyss_1.3.3/paths.bash

GapCloser -a gapclose_in.fna -b gapclose.conf.header -t 15 -o gapclosed.fna > log_gapclose.out

