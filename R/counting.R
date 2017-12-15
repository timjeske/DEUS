#' Function to create sequence-based count tables
#'
#' This function creates count tables summarizing counts for each observed read in a set of fastQ files.
#' Input files are taken from the given input directory and are selected based on the given fastQ suffix.
#' @param in_dir input directory containing fastq files to be analyzed
#' @param fq_suffix suffix of your fastq files, default is ".fq.gz"
#' @param phenoInfo by default count table headers are the basenames of the input fastq files. If a file is given with sample names in first column and a sample column (=fastq path), the sample names will be used as headers in the count table.
#' @keywords read counting
#' @export
#' @examples
#' in_dir <- system.file("extdata", package = "USBseq")
#' countTable <- createCountTableFromFastQs(in_dir)

createCountTableFromFastQs <- function(in_dir,fq_suffix = ".fq.gz",phenoInfo=NULL) {
  file.list <- list.files(in_dir)
  fastqs <- file.list[endsWith(file.list,fq_suffix)]
  fastqs <- paste(in_dir,fastqs,sep="/")

  print(paste("Counting all reads in",fastqs[1]))
  content <- ShortRead::readFastq(fastqs[1])
  counts <- ShortRead::tables(content,n=length(content))
  countTable <- as.data.frame(counts$top)

  names(countTable) <- c(basename(fastqs[1]))
  if(!is.null(phenoInfo) && "sample" %in% colnames(phenoInfo)) {
    phenoInfo$sampleName = rownames(phenoInfo)
    phenoInfo$baseName = basename(as.character(phenoInfo$sample))
    names(countTable) <- c(as.character(phenoInfo$sampleName[basename(fastqs[1])==phenoInfo$baseName]))
  }

  countTable$Read <- row.names(countTable)
  fastqs <- fastqs[-1]
  for(f in fastqs) {
    print(paste("Counting all reads in",f))
    content <- ShortRead::readFastq(f)
    counts <- ShortRead::tables(content,n=length(content))
    counts <- as.data.frame(counts$top)
    names(counts) <- c(basename(f))
    if(!is.null(phenoInfo) && "sample" %in% colnames(phenoInfo)) {
      names(counts) <- c(as.character(phenoInfo$sampleName[basename(f)==phenoInfo$baseName]))
    }
    counts$Read <- row.names(counts)
    countTable <- plyr::join(countTable, counts, type="full")
  }
  row.names(countTable) <- countTable$Read
  countTable$Read <- NULL
  countTable[is.na(countTable)] <- 0
  return(countTable)
}
