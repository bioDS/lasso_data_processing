#!/bin/sh

# Run TargetScan's base script (targetscan_60.pl) on all seed files.
#
# Author: Fabian Schmich (fabian.schmich@bsse.ethz.ch)
#

threads=16

SEEDS_DIR="/home/kieran/work/infx_lasso_data/seeds"
TS_UTRSEQ="./data/hg19_ucsc_3p.txt"
TSCAN_60="/home/kieran/work/tscan60/targetscan_60/targetscan_60.pl"
OUT_DIR="/home/kieran/work/infx_lasso_data/tscan60"

for s in $SEEDS_DIR/*.txt
do
	bname=$(basename $s)
	if [ ! -f $OUT_DIR/tscan.$bname ]; then
		((i=i%threads)); ((i++==0)) && wait
		{ perl $TSCAN_60 $s $TS_UTRSEQ $OUT_DIR/tscan.$bname || true; } &
	fi
done
