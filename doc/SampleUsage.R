## ----gh-installation, message=FALSE--------------------------------------
if (!require("devtools")) install.packages("devtools", repos='http://cran.us.r-project.org')
devtools::install_github("timjeske/DEUS")
library(DEUS)

## ---- eval=FALSE---------------------------------------------------------
#  if (!require("devtools")) install.packages("devtools", repos='http://cran.us.r-project.org')
#  devtools::install_local("./DEUS-master")
#  library(DEUS)

## ---- eval=TRUE----------------------------------------------------------
in_dir <- system.file("extdata", package = "DEUS")
out_dir <- "~/tmp"
phenofile <- system.file("extdata", "condition.tsv", package = "DEUS")
pheno_info <- read.table(phenofile, header=T, row.names=1, check.names=FALSE)
pheno_info

## ---- eval=TRUE, message=FALSE-------------------------------------------
count_table <- createCountTableFromFastQs(in_dir, pheno_info=pheno_info)
count_table_filtered <- filterLowExp(count_table, pheno_info)
head(count_table_filtered, n=2)

## ---- eval=TRUE----------------------------------------------------------
map <- createMap(count_table_filtered)
head(map, n=2)

## ---- eval=TRUE, message=FALSE-------------------------------------------
design <- ~ condition
de_results <- runDESeq2(count_table_filtered, pheno_info, design, out_dir=out_dir)
sig_results <- de_results$de_result
sig_results <- sig_results[!is.na(sig_results$IHWPvalue) & sig_results$IHWPvalue < 0.05,]
head(sig_results, n=2)

## ---- eval=FALSE---------------------------------------------------------
#  blast_exec <- "/your/path/to/ncbi-blast-x.x.x/bin/blastn"
#  blast_db <-system.file("extdata", "blastdb/DASHR_subset.fa", package = "DEUS")
#  blast_ncores <- 2
#  sig_seq_fasta <- sequencesAsFasta(sig_results,map)
#  blast_result <- runBlast(blast_exec, blast_db, blast_ncores, sig_seq_fasta)

## ---- eval=TRUE----------------------------------------------------------
blast_int <- system.file("extdata", "results/blast_result_sig_sequences.tsv", package="DEUS")
blast_result <- read.table(blast_int, header=T, sep="\t")
head(blast_result)

## ---- eval=TRUE----------------------------------------------------------
count_stats <- getConditionCountStats(de_results$norm_counts, pheno_info)
head(count_stats, n=2)

## ---- eval=FALSE---------------------------------------------------------
#  cd_hit <- "/your/path/to/cdhit/cd-hit-est"
#  sig_seq_fasta <- sequencesAsFasta(sig_results,map)
#  clust_result <- runClustering(cd_hit, sig_seq_fasta, out_dir, identity_cutoff=0.9, length_cutoff=0.9, wordlength=9, map)

## ---- eval=TRUE----------------------------------------------------------
clust_int <- system.file("extdata", "results/clust_result_sig_sequences.tsv", package="DEUS")
clust_result <- read.table(clust_int, header=T, sep="\t")
rownames(clust_result) <- clust_result$SequenceID
head(clust_result)

## ---- eval=TRUE,message=FALSE--------------------------------------------
classes <- c("tRNA","[Hh]sa","^U")
summary <- mergeResults(sig_results, count_stats, blast_result, clust_result, map)
summary <- addCountsOfFeatureClasses(summary, classes)
writeSummaryFiles(summary, out_dir)
head(summary, n=2)

## ---- eval=TRUE----------------------------------------------------------
getNoBlastHitFraction(summary, de_results$norm_counts)

## ---- eval=TRUE----------------------------------------------------------
map <- createMap(count_table)

## ---- eval=FALSE---------------------------------------------------------
#  cd_hit <- "/your/path/to/cdhit/cd-hit-est"
#  all_seq_fasta <- sequencesAsFasta(count_table,map)
#  clust_result <- runClustering(cd_hit, all_seq_fasta, out_dir, identity_cutoff=0.9, length_cutoff=0.9, wordlength=9, map)
#  cl_counts <- mergeAndAggregate(map, count_table, clust_result)

## ---- eval=TRUE,message=FALSE--------------------------------------------
clust_int <- system.file("extdata", "results/clust_result_clust_sequences.tsv", package="DEUS")
clust_result <- read.table(clust_int, header=T, sep="\t")
rownames(clust_result) <- clust_result$SequenceID
cl_counts <- mergeAndAggregate(map, count_table, clust_result)
head(cl_counts, n=2)

## ---- eval=TRUE, message=FALSE-------------------------------------------
design <- ~ condition
cl_de_results <- runDESeq2(cl_counts, pheno_info, design, out_dir = out_dir, prefix = "Cluster")
cl_sig_results <- cl_de_results$de_result
head(cl_sig_results, n=2)

## ---- eval=TRUE----------------------------------------------------------
sig_results <- mergeSingleAndClusterResults(cl_sig_results,clust_result,sig_results,map)
sig_results <- sig_results[((!is.na(sig_results$IHWPval) & sig_results$IHWPval < 0.05) | (!is.na(sig_results$Cl_IHWPval) & sig_results$Cl_IHWPval < 0.05)),]
head(sig_results, n=2)

## ---- eval=FALSE---------------------------------------------------------
#  blast_exec <- "/your/path/to/ncbi-blast-x.x.x/bin/blastn"
#  blast_db <-system.file("extdata", "blastdb/DASHR_subset.fa", package = "DEUS")
#  blast_ncores <- 2
#  sig_seq_fasta <- sequencesAsFasta(sig_results,map)
#  blast_result <- runBlast(blast_exec, blast_db, blast_ncores, sig_seq_fasta)

## ---- eval=TRUE----------------------------------------------------------
blast_int <- system.file("extdata", "results/blast_result_clust_sequences.tsv", package="DEUS")
blast_result <- read.table(blast_int, header=T, sep="\t")
head(blast_result)

## ---- eval=TRUE, message=FALSE-------------------------------------------
# merge results
classes <- c("tRNA","[Hh]sa","^U")
summary <- mergeResults(sig_results, count_stats, blast_result, clust_result, map)
summary <- addCountsOfFeatureClasses(summary, classes)
writeSummaryFiles(summary, out_dir)
head(summary, n=2)

## ---- eval=TRUE----------------------------------------------------------
printClusterSummary(summary)

