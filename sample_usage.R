library(devtools)

# install package
#install_github("timjeske/USBseq", auth_token = "856ab1ec38c789a1a0ac30cdfbcd8ccf1c6f224f")
setwd("~")
install("USBseq")
library(USBseq)

# set input and output
in_dir <- system.file("extdata", package = "USBseq")
out_dir <- "/storageNGS/ngs4/projects/sncRNA_USB/pipeline/2018_04_09_sample_usage"
phenofile <- system.file("extdata", "condition_test.tsv", package = "USBseq")
phenoInfo <- read.table(phenofile, header=T, row.names=1, check.names=FALSE)

# create and filter count table, create sequence to sequenceID map
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
blast_exec <- "/storageNGS/ngs1/software/ncbi-blast-2.6.0+/bin/blastn"
blast_db <-"/storageNGS/ngs1/software/ncbi-blast-2.6.0+/blastdb/MouseDB_piRNABank.fa"
ncores <- 2
blastResult <- runBlast(blast_exec, blast_db, ncores, sigSeqFasta, identity = 95)
write.table(blastResult, paste(out_dir, "Sig_sequences.blastn.tsv",sep="/"), col.names=T, quote=F, sep="\t", row.names=F)

# run clustering
cd_hit <- "/storageNGS/ngs1/software/cdhit/cd-hit-est"
clustResult<-runClustering(cd_hit, sigSeqFasta, out_dir, 0.9, 0.9, 9, map)

# merge results
classes <- c("mmu_piR","ENSMUST","tRNA","mmu-miR","retro")
summary <- mergeResults(sigResults,countStats, blastResult, clustResult, map)
summary <- addCountsOfFeatureClasses(summary, classes)
writeSummaryFiles(summary,out_dir)
