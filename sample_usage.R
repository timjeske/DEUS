library(devtools)

# install package
install_github("timjeske/USBseq", auth_token = "856ab1ec38c789a1a0ac30cdfbcd8ccf1c6f224f")
#setwd("~")
#install("USBseq")
library(USBseq)

# defaults may be changed
out_dir <- "~/USBseq/inst/extdata/results"
blast_ncores <- 2

# leave paths empty if not yet installed, intermediate results will be used for testing purposes
blast_exec <- "/storageNGS/ngs1/software/ncbi-blast-2.6.0+/bin/blastn"
cd_hit <- "/storageNGS/ngs1/software/cdhit/cd-hit-est"

# load data delivered with the package
in_dir <- system.file("extdata", package = "USBseq")
phenofile <- system.file("extdata", "condition.tsv", package = "USBseq")
blast_db <-system.file("extdata", "blastdb/DASHR_subset.fa", package = "USBseq")
blast_int <- system.file("extdata", "results/Sig_sequences.blastn.tsv", package="USBseq")
clust_int <- system.file("extdata", "results/clustResult.tsv", package="USBseq")

# create and filter count table, create sequence to sequenceID map
phenoInfo <- read.table(phenofile, header=T, row.names=1, check.names=FALSE)
countTable <- createCountTableFromFastQs(in_dir, phenoInfo=phenoInfo)
countTable <- filterLowExp(countTable, phenoInfo)
write.table(countTable, paste(out_dir,"AllCounts_filtered.tsv",sep="/"), col.names=T, quote=F, sep="\t", row.names=T)

# run differential expression analysis
design <- ~ condition
deResults <- runDESeq2(countTable, phenoInfo, design, map, out_dir)
sigResults <- deResults$deResult
sigResults <- sigResults[!is.na(sigResults$IHWPval) & sigResults$IHWPval < 0.05,]
map <- createMap(sigResults)
sigSeqFasta <- sequencesAsFasta(sigResults,map)

# get count stats
countStats <- getConditionCountStats(deResults$normCounts, phenoInfo)

# run blast
if(blast_exec == "") {
  blastResult <- read.table(blast_int, header=T, sep="\t")
} else {
  blastResult <- runBlast(blast_exec, blast_db, blast_ncores, sigSeqFasta)
}

# run clustering
if(cd_hit == "") {
  clustResult <- read.table(clust_int, header=T, sep="\t")
} else {
  clustResult<-runClustering(cd_hit, sigSeqFasta, out_dir, 0.9, 0.9, 9, map)
}

# merge results
classes <- c("tRNA","[Hh]sa","^U")
summary <- mergeResults(sigResults, countStats, blastResult, clustResult, map)
summary <- addCountsOfFeatureClasses(summary, classes)
writeSummaryFiles(summary, out_dir)

