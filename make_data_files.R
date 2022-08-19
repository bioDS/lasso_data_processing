#!/usr/bin/Rscript

library(normInfectX)
library(dplyr)

run_predicttargets_script = "3_custom_run_predicttargets.arguments"

#write("OUT_DIR=\"/home/kieran/work/infx_lasso_data/targetpredict\"", file=run_predicttargets_script)
cat("", file=run_predicttargets_script)

rows = function(x) lapply(seq_len(nrow(x)), function(i) lapply(x,"[",i))

mock_data =  mock %>% filter(WellType != "CONTROL") %>% filter(!is.na(ID)) %>% select(Catalog_number, Sequence_antisense_5_3) %>% mutate(Catalog_number = as.character(Catalog_number))
vaccinia_data = vaccinia %>% filter(WellType != "CONTROL") %>% filter(!is.na(ID)) %>% select(Catalog_number, Sequence_antisense_5_3) %>% mutate(Catalog_number = as.character(Catalog_number))


t = union(mock_data, vaccinia_data)

t$Catalog_number = make.unique(t$Catalog_number)

for (entry in rows(t)) {
	sirna_id =  as.character(entry$Catalog_number)
	sequences = strsplit(as.character(entry$Sequence_antisense_5_3), ' ')
	for (s in sequences) {
		#print(sprintf("python predicttargets.py --id %s --seq %s --out $OUT_DIR", sirna_id, s))
		write(sprintf("--id %s --seq %s", sirna_id, s), file=run_predicttargets_script, append=TRUE)
	}
}
