library(devtools)

#install_github("timjeske/USBseq", auth_token = "856ab1ec38c789a1a0ac30cdfbcd8ccf1c6f224f")

setwd("~")
install("USBseq")
library(USBseq)

in_dir <- system.file("extdata", package = "USBseq")
phenofile <- system.file("extdata", "condition.tsv", package = "USBseq")
phenoInfo <- read.table(phenofile, header=T, row.names=1, check.names=FALSE)

countTable <- createCountTableFromFastQs(in_dir)
countDataFilt <- filterLowExp(countTable, phenoInfo)
#write.table(countDataFilt, "~/tmp/AllCounts_filtered.tsv", col.names=T, quote=F, sep="\t", row.names=T)

#do differential expression analysis
countData <- countDataFilt
colData <- phenoInfo
map <- createMap(countTable)
design <- ~ condition
filter_val <- 'IHWPval'
filter_threshold <- 0.05
out_dir <- "~/tmp/"
source("~/Dropbox/Scripts/DiffExpression/DESeq2.R")
sigResults <- read.table("~/tmp/DESeq2_sig_output.tsv", head=T, row.names=1)

#run blast
blast_exec <-...
blast_db <-...
ncores <- 2
blastResult <- runBlast(blast_exec, blast_db, ncores, sigResults, map)
