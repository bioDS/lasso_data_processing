#!/usr/bin/env bash

run_sets_on_methods() {
  for set in $sets; do
    echo "benchmarking collection $set"
    root_dir="whinter_pint_comparison/$set/"
    directories="$root_dir $root_dir/summaries"
    for dir in directories; do
      if ! [ -d $dir ]; then
        mkdir -p whinter_pint_comparison
      fi
    done

    bench_sets_dir="~/work/data/simulated_rerun/$set/"
    for bench_set in `bash -c "cd $bench_sets_dir; ls *rds"`; do
      echo "benchmarking $bench_set"
      set_dir="$root_dir/summaries/${bench_set/\.rds/}"
      if ! [ -d $set_dir ]; then
        mkdir -p $set_dir
      fi
      if ! [ -f $set_dir/all.rds ]; then
        taskset -c $taskset_threads bash -c "./process_rds_whinter.R $bench_sets_dir/$bench_set $set_dir $methods $num_features $depth $use_cores" &> $set_dir/log.txt
      else
        echo "skipping $set_dir"
      fi
    done
  done
}

#set="1k_whinter_working_data"
#set="simulated_small_data"

# to reproduce current paper figures:
## three-way comparison
num_features=5000
sets="3way"
methods="all"
depth="3"
taskset_threads="0-95"
use_cores="48"
run_sets_on_methods

depth=2
num_features=1000
sets="simulated_small_data_sample"
methods="all"
taskset_threads="0"
use_cores="1"
run_sets_on_methods

depth=2
num_features=5000
sets="8k_only"
methods="all"
taskset_threads="0-95"
use_cores="48"
run_sets_on_methods

## new wide comparison
num_features=5000
sets="wide_only_10k"
methods="all"
taskset_threads="0-95"
use_cores="48"
run_sets_on_methods


# testing extra things
# num_features=60
# sets="simulated_small_data_sample"
# methods="noglint"
# taskset -c 0 run_sets_on_methods
