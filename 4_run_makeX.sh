#!/bin/sh

# Run makeX.py on the output directory of predicted context+ scores for siRNA off-targets.
#
# Author: Fabian Schmich (fabian.schmich@bsse.ethz.ch)
#

CPS_DIR="/home/kieran/work/infx_lasso_data/targetpredict"
SIRNAS="./siRNAs.txt"
GENES="./genes.txt"

OUT="/home/kieran/work/infx_lasso_data/xmatrix_vaccinia"

#python2 makeX.py --cps $CPS_DIR --sirnas $SIRNAS --genes $GENES --out $OUT
python2 makeX.py --cps $CPS_DIR --sirnas $SIRNAS --genes $GENES --out $OUT
