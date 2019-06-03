#*  Copyright (C) 2018 the DEUS contributors.
#*  Website: https://github.com/timjeske/DEUS
#*
#*  This file is part of the DEUS R package.
#*
#*  The DEUS R package is free software: you can redistribute it and/or modify
#*  it under the terms of the GNU General Public License as published by
#*  the Free Software Foundation, either version 3 of the License, or
#*  (at your option) any later version.
#*
#*  This program is distributed in the hope that it will be useful,
#*  but WITHOUT ANY WARRANTY; without even the implied warranty of
#*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#*  GNU General Public License for more details.
#*
#*  You should have received a copy of the GNU General Public License
#*  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# install package
if (!require("devtools")) install.packages("devtools", repos='http://cran.us.r-project.org')
devtools::install_github("timjeske/DEUS")
library(DEUS)

# defaults may be changed
out_dir <- tempdir()
blast_ncores <- 2

# leave paths empty if not yet installed, intermediate results will be used for testing purposes
blast_exec <- "/data/Software/ncbi-blast-2.6.0+/bin/blastn"
#blast_exec <- ""
cd_hit <- "/data/Software/cdhit/cd-hit-est"
#cd_hit <- ""

# load data delivered with the package
in_dir <- system.file("extdata", package = "DEUS")
phenofile <- system.file("extdata", "condition.tsv", package = "DEUS")
blast_db <-system.file("extdata", "blastdb/DASHR_subset.fa", package = "DEUS")
blast_int <- system.file("extdata", "results/blast_result_sig_sequences.tsv", package="DEUS")
clust_int <- system.file("extdata", "results/clust_result_sig_sequences.tsv", package="DEUS")

# create and filter count table, create sequence to sequenceID map
pheno_info <- read.table(phenofile, header=T, row.names=1, check.names=FALSE)
count_table <- createCountTableFromFastQs(in_dir, pheno_info=pheno_info)
count_table <- filterLowExp(count_table, pheno_info)
map <- createMap(count_table)
write.table(count_table, paste(out_dir,"AllCounts_filtered.tsv",sep="/"), col.names=T, quote=F, sep="\t", row.names=T)

# run differential expression analysis
design <- ~ condition
de_results <- runDESeq2(count_table, pheno_info, design, out_dir=out_dir)
sig_results <- de_results$de_result
sig_results <- sig_results[!is.na(sig_results$IHWPvalue) & sig_results$IHWPvalue < 0.05,]
sig_seq_fasta <- sequencesAsFasta(sig_results,map)

# get count stats
count_stats <- getConditionCountStats(de_results$norm_counts, pheno_info)

# run blast
if(file.exists(blast_exec)) {
  blast_result <- runBlast(blast_exec, blast_db, blast_ncores, sig_seq_fasta)

  # write intermediate result file
  git_path <- paste(getwd(),"inst/extdata/results", sep="/")
  if(file.exists(git_path)) {
    blast_result <- blast_result[c("SequenceID", "Annotation", "Length", "BlastEvalue")]
    write.table(blast_result, paste(git_path, "blast_result_sig_sequences.tsv", sep="/"), col.names=T, sep="\t", quote=F, row.names=F)
  }
} else {
  blast_result <- read.table(blast_int, header=T, sep="\t")
}

# run clustering
if(file.exists(cd_hit)) {
  clust_result <- runClustering(cd_hit, sig_seq_fasta, out_dir, identity_cutoff=0.9, length_cutoff=0.9, wordlength=9, map)

  # write intermediate result file
  git_path <- paste(getwd(),"inst/extdata/results", sep="/")
  if(file.exists(git_path)) {
    write.table(clust_result, paste(git_path, "clust_result_sig_sequences.tsv", sep="/"), col.names=T, sep="\t", quote=F, row.names=F)
  }
} else {
  clust_result <- read.table(clust_int, header=T, sep="\t")
  rownames(clust_result) <- clust_result$SequenceID
}

# merge results
classes <- c("tRNA","[Hh]sa","^U")
summary <- mergeResults(sig_results, count_stats, blast_result, clust_result, map)
summary <- addCountsOfFeatureClasses(summary, classes)
writeSummaryFiles(summary, out_dir)

# compute NA fraction
notAnnotated <- getNoBlastHitFraction(summary, de_results$norm_counts)
print(notAnnotated)

