#!/usr/bin/env Rscript
require(Matrix)
require(dplyr)
require(normInfectX)

ds = clean(pathogen="vaccinia", killers=FALSE, controls=FALSE, quality="bad")[["vaccinia"]] %>% filter(!is.na(eCount_oCells))
lethals = (ds %>% filter(eCount_oCells == 0))

barcodes = (vaccinia %>% filter(Catalog_number %in% lethals$Catalog_number))$Barcode

print("lethal genes are contained within barcodes:")
table(as.character(barcodes))

for (barcode in unique(barcodes)) {
	print("summarissing cell counts in barcode")
	print(as.character(barcode))

	print(summary((vaccinia %>% filter(Barcode == barcode))$eCount_oCells))
}
