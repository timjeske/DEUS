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
out_dir <- "/storageNGS/ngs4/projects/other/Mouse_beckers_huypens/smallRNA/Pipeline/2017_11_16_checkRefactoring"
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
deResults <- runDESeq2(countData, colData, design, map, out_dir)
sigResults <- deResults$deResult
sigResults <- sigResults[!is.na(sigResults$IHWPval) & sigResults$IHWPval < 0.05,]

# get count stats
countStats <- getCountStats(deResults$normCounts, colData)

# run blast
blast_exec <- "/storageNGS/ngs1/software/ncbi-blast-2.6.0+/bin/blastn"
blast_db <-"/storageNGS/ngs1/software/ncbi-blast-2.6.0+/blastdb/MouseDB2.fa"
ncores <- 2
blastResult <- runBlast(blast_exec, blast_db, ncores, sigResults, map,out_dir)
write.table(blastResult, paste(out_dir, "Sig_sequences.blastn.tsv",sep="/"), col.names=T, quote=F, sep="\t", row.names=F)

cd_hit <- "/storageNGS/ngs1/software/cdhit/cd-hit-est"
seq_fasta <- paste(out_dir,"sig_sequences.fa",sep="/")
write.table(sequencesAsFasta(sigResults,map),seq_fasta,quote = F,row.names = F,col.names = F)
clustResult<-runClustering(cd_hit,seq_fasta,out_dir,0.9,0.9,9,map)

# merge results
summary_err <- mergeResults(sigResults,map=map)
summary_blast <- mergeResults(sigResults,blastResult=blastResult,map=map)
summary_clust <- mergeResults(sigResults,clustResult=clustResult,map=map)

summary <- mergeResults(sigResults,countStats, blastResult, clustResult, map)
summary <- addCountsOfFeatureClasses(summary, classes)
summary_err <- addCountsOfFeatureClasses(summary_clust, classes)
writeSummaryFiles(summary,out_dir)

#Remove tmp files
deleteTmp(out_dir)

