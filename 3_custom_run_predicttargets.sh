#!/bin/bash

export OUT_DIR="/home/kieran/work/infx_lasso_data/targetpredict"

parallel -j 16 -a 3_custom_run_predicttargets.arguments python /home/kieran/work/infx_lasso_data/predicttargets.py --out $OUT_DIR {=uq=}
