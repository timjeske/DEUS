library(devtools)

# install package
install_github("timjeske/USBseq", auth_token = "856ab1ec38c789a1a0ac30cdfbcd8ccf1c6f224f")
setwd("~")
#install("USBseq")
library(USBseq)

# set input and output
in_dir <- system.file("extdata", package = "USBseq")
out_dir <- "/storageNGS/ngs4/projects/sncRNA_USB/pipeline/2018_04_09_sample_usage_getStats"
phenofile <- system.file("extdata", "condition_test.tsv", package = "USBseq")
phenoInfo <- read.table(phenofile, header=T, row.names=1, check.names=FALSE)

# create and filter count table, create sequence to sequenceID map
countTable <- createCountTableFromFastQs(in_dir, phenoInfo=phenoInfo)
countTable <- filterLowExp(countTable, phenoInfo)
map <- createMap(countTable)
sigSeqFasta <- sequencesAsFasta(countTable,map)

# run blast
blast_exec <- "/storageNGS/ngs1/software/ncbi-blast-2.6.0+/bin/blastn"
blast_db <-"/storageNGS/ngs1/software/ncbi-blast-2.6.0+/blastdb/MouseDB_piRNABank.fa"
ncores <- 16
blastResult <- runBlast(blast_exec, blast_db, ncores, sigSeqFasta, identity = 95)
write.table(blastResult, paste(out_dir, "Sig_sequences.blastn.tsv",sep="/"), col.names=T, quote=F, sep="\t", row.names=F)

# merge results
classes <- c("mmu_piR","ENSMUST","tRNA","mmu-miR","retro")
summary <- mergeResults(blastResult=blastResult,map=map)
summary <- addCountsOfFeatureClasses(summary, classes)
writeSummaryFiles(summary,out_dir)

generateSummaryPlots(summary,classes, out_dir)

res <- getNoBlastHitFraction(summary, deResults$normCounts)
write.table(res,paste(out_dir,"NA_fraction_per_sample.tsv",sep="/"),row.names = T, col.names = T, quote=F, sep="\t")

