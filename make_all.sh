#!/bin/bash

./make_data_files.R
cut -d' ' -f 2 3_custom_run_predicttargets.arguments > siRNAs.txt
rg --text --no-filename --only-matching "^[1-9]\d*" targetpredict/ | sort -n | uniq > genes.txt
