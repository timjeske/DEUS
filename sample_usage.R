library(devtools)

# install package
#install_github("timjeske/USBseq", auth_token = "856ab1ec38c789a1a0ac30cdfbcd8ccf1c6f224f")
setwd("~")
install("USBseq")
library(USBseq)

# set input and output
in_dir <- system.file("extdata", package = "USBseq")
out_dir <- "/storageNGS/ngs4/projects/sncRNA_USB/pipeline/2018_01_18_summaryTableStats"
phenofile <- system.file("extdata", "condition.tsv", package = "USBseq")
phenoInfo <- read.table(phenofile, header=T, row.names=1, check.names=FALSE)

# create and filter count table, create sequence to sequenceID map
countTable <- createCountTableFromFastQs(in_dir, phenoInfo=phenoInfo)
map <- createMap(countTable)
countDataFilt <- filterLowExp(countTable, phenoInfo)
write.table(countDataFilt, paste(out_dir,"AllCounts_filtered.tsv",sep="/"), col.names=T, quote=F, sep="\t", row.names=T)

# run differential expression analysis
design <- ~ condition
deResults <- runDESeq2(countDataFilt, phenoInfo, design, map, out_dir)
sigResults <- deResults$deResult
sigResults <- sigResults[!is.na(sigResults$IHWPval) & sigResults$IHWPval < 0.05,]
sigSeqFasta <- sequencesAsFasta(sigResults,map)

# get count stats
countStats <- getConditionCountStats(deResults$normCounts, phenoInfo)

# run blast
blast_exec <- "/storageNGS/ngs1/software/ncbi-blast-2.6.0+/bin/blastn"
blast_db <-"/storageNGS/ngs1/software/ncbi-blast-2.6.0+/blastdb/MouseDB2.fa"
ncores <- 2
blastResult <- runBlast(blast_exec, blast_db, ncores, sigSeqFasta)
write.table(blastResult, paste(out_dir, "Sig_sequences.blastn.tsv",sep="/"), col.names=T, quote=F, sep="\t", row.names=F)

# run clustering
cd_hit <- "/storageNGS/ngs1/software/cdhit/cd-hit-est"
seq_fasta <- paste(out_dir,"sig_sequences.fa",sep="/")
write.table(sigSeqFasta,seq_fasta,quote = F,row.names = F,col.names = F)
clustResult<-runClustering(cd_hit,seq_fasta,out_dir,0.9,0.9,9,map)

# merge results
classes <- c("piR-mmu","ENSMUST","tRNA","mmu-miR","retro")
summary <- mergeResults(sigResults,countStats, blastResult, clustResult, map)
summary <- addCountsOfFeatureClasses(summary, classes)
writeSummaryFiles(summary,out_dir)

#Remove tmp files
deleteTmp(out_dir)

generateSummaryStats(summary, phenoInfo, classes)
