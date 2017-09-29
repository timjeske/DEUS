library(devtools)

# install package
#install_github("timjeske/USBseq", auth_token = "856ab1ec38c789a1a0ac30cdfbcd8ccf1c6f224f")
setwd("~")
install("USBseq")
library(USBseq)

# set input and output
in_dir <- system.file("extdata", package = "USBseq")
phenofile <- system.file("extdata", "condition.tsv", package = "USBseq")
phenoInfo <- read.table(phenofile, header=T, row.names=1, check.names=FALSE)
out_dir <- "/storageNGS/ngs4/projects/other/Mouse_beckers_huypens/smallRNA/Pipeline/2017_09_28_Rpackage"
classes <- c("piR-mmu","ENSMUST","tRNA","mmu-miR","retro")

# create count table
countTable <- createCountTableFromFastQs(in_dir)
countDataFilt <- filterLowExp(countTable, phenoInfo)
write.table(countDataFilt, paste(out_dir,"AllCounts_filtered.tsv",sep="/"), col.names=T, quote=F, sep="\t", row.names=T)

# do differential expression analysis
countData <- countDataFilt
colData <- phenoInfo
map <- createMap(countTable)
design <- ~ condition
filter_val <- 'IHWPval'
filter_threshold <- 0.05
source("~/USBseq/DESeq2.R")
sigResults <- read.table(paste(out_dir,"DESeq2_sig_output.tsv",sep="/"), head=T, row.names=1)

# run blast
blast_exec <- "/storageNGS/ngs1/software/ncbi-blast-2.6.0+/bin/blastn"
blast_db <-"/storageNGS/ngs1/software/ncbi-blast-2.6.0+/blastdb/MouseDB2.fa"
ncores <- 2
blastResult <- runBlast(blast_exec, blast_db, ncores, sigResults, map)
write.table(blastResult, paste(out_dir, "Sig_sequences.blastn.tsv",sep="/"), col.names=T, quote=F, sep="\t", row.names=F)

# merge results
summary <- mergeResults(sigResults, blastResult, map)
summary <- addCountsOfFeatureClasses(summary, classes)
