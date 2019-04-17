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
out_dir <- "~/tmp"
blast_ncores <- 2

# leave paths empty if not yet installed, intermediate results will be used for testing purposes
#blast_exec <- "/storageNGS/ngs1/software/ncbi-blast-2.6.0+/bin/blastn"
blast_exec <- "/data/Software/ncbi-blast-2.8.1+/bin/blastn"
#cd_hit <- "/storageNGS/ngs1/software/cdhit/cd-hit-est"
cd_hit <- "/data/Software/cdhit/cd-hit-est"

# load data delivered with the package
in_dir <- system.file("extdata", package = "DEUS")
phenofile <- system.file("extdata", "condition.tsv", package = "DEUS")
blast_db <-system.file("extdata", "blastdb/DASHR_subset.fa", package = "DEUS")
blast_int <- system.file("extdata", "results/Sig_sequences.blastn.tsv", package="DEUS")
clust_int <- system.file("extdata", "results/clustResult.tsv", package="DEUS")

# create and filter count table, create sequence to sequenceID map
pheno_info <- read.table(phenofile, header=T, row.names=1, check.names=FALSE)
count_table <- createCountTableFromFastQs(in_dir, pheno_info=pheno_info)

#Create map with all sequences
map <- createMap(count_table)

#<<<<<<<<
#Cluster Approach first

#cl_map <- map
#Get all sequences as one fasta file
allSeqFasta <- sequencesAsFasta(count_table,map)
#Cluster them
clustResult<-runClustering(cd_hit, allSeqFasta, out_dir, 0.9, 0.9, 9, map)
#Aggregate counts by cluster
cl_counts <- mergeAndAggregate(map,count_table,clustResult)

# run differential expression analysis on clusters
design <- ~ condition
cl_deResults <- runDESeq2(cl_counts, pheno_info, design, out_dir = out_dir)
cl_sigResults <- cl_deResults$deResult
#>>>>>>>>>

#Continue with regular pipeline

#Filter lowexpressed but only for single seq approach
count_table <- filterLowExp(count_table, pheno_info)
write.table(count_table, paste(out_dir,"AllCounts_filtered.tsv",sep="/"), col.names=T, quote=F, sep="\t", row.names=T)

# run differential expression analysis
deResults <- runDESeq2(count_table, pheno_info, design, out_dir=out_dir)
sigResults <- deResults$deResult

# get count stats
countStats <- getConditionCountStats(deResults$normCounts, pheno_info)

###########################
#Merge sigResults and cl_sigResults
sigResults <- mergeSingleAndClusterResults(cl_sigResults,clustResult,sigResults,map)
##########################
sigResults <- sigResults[((!is.na(sigResults$IHWPval) & sigResults$IHWPval < 0.05) | (!is.na(sigResults$Cl_IHWPval) & sigResults$Cl_IHWPval < 0.05)),]

#Get sequences for blast
sigSeqFasta <- sequencesAsFasta(sigResults,map)
# run blast
if(file.exists(blast_exec)) {
  blastResult <- runBlast(blast_exec, blast_db, blast_ncores, sigSeqFasta)
} else {
  blastResult <- read.table(blast_int, header=T, sep="\t")
}

# merge results
classes <- c("tRNA","[Hh]sa","^U")
summary <- mergeResults(sigResults, countStats, blastResult, clustResult, map)
summary <- addCountsOfFeatureClasses(summary, classes)
writeSummaryFiles(summary, out_dir)
