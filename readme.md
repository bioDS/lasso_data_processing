# Benchmarking and data processing scripts.

N.B. The infectX and antibiotic data are not included.

## Reproducing simulations

A compiled version of WHInter should be placed in `./WHInter`.

```
Rscript regenerate_data.R
bash bench_pint_whinter.sh
Rscript collect_bench_stats.R
Rscript plot.R
```

## InfectX

Requires the package `normInfectX`.

```
Rscript process_ifx.R
```

## Antibiotic data

Assumes a table of SNVs and MIC is present in "binary_table_T_anonymised.csv".

```
./process_bacteria_data.R
./process_antibio_ifx.R
```
