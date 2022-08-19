count=1
for cores in $(seq 0 7); do
  taskset -c 0-$cores ./pint_only_time.R &> "time_$cores.txt"
  mv time.rds "time_large_$count.rds"
  count=$((count+1))
done
for cores in $(seq 16 23); do
  taskset -c 0-7,16-$cores ./pint_only_time.R &> "time_$cores.txt"
  mv time.rds "time_large_$count.rds"
  count=$((count+1))
done
for cores in $(seq 8 15); do
  taskset -c 0-7,16-23,8-$cores ./pint_only_time.R &> "time_$cores.txt"
  mv time.rds "time_large_$count.rds"
  count=$((count+1))
done
for cores in $(seq 24 31); do
  taskset -c 0-$cores ./pint_only_time.R &> "time_$cores.txt"
  mv time.rds "time_large_$count.rds"
  count=$((count+1))
done

#taskset -c 0-0 ./pint_only_time.R &> "time_0.txt"
#mv time.rds "time_small_0.rds"
#taskset -c 0,1 ./pint_only_time.R &> "time_small_1.txt"
#mv time.rds "time_small_1.rds"
#taskset -c 0,1,2 ./pint_only_time.R &> "time_small_2.txt"
#mv time.rds "time_small_2.rds"
#taskset -c 0,1,2,3 ./pint_only_time.R &> "time_small_3.txt"
#mv time.rds "time_small_3.rds"
#taskset -c 0-7 ./pint_only_time.R &> "time_small_7.txt"
#mv time.rds "time_small_7.rds"
