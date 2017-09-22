library(devtools)

#install_github("timjeske/USBseq", auth_token = "856ab1ec38c789a1a0ac30cdfbcd8ccf1c6f224f")

setwd("~")
install("USBseq")

library(USBseq)
in_dir <- system.file("extdata", package = "USBseq")
countTable <- createCountTableFromFastQs(in_dir)

#write.table(countTable, "~/tmp/AllCounts_R.tsv", col.names=T, quote=F, sep="\t", row.names=T)

map <- createMap(countTable)

phenofile <- system.file("extdata", "condition.tsv", package = "USBseq")
countDataFilt <- filterLowExp(countTable, phenofile)
