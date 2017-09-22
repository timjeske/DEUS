#' Function to create sequence-based count tables
#'
#' This function creates count tables summarizing counts for each observed read in a set of fastQ files.
#' Input files are taken from the given input directory and are selected based on the given fastQ suffix.
#' @param in_dir input directory containing fastq files to be analyzed
#' @param fq_suffix suffix of your fastq files, default is ".fq.gz"
#' @keywords
#' @export
#' @examples
#' createCountTableFromFastQs("./my_experiment/fastQs",".trimmed.fastq.gz")
#' createCountTableFromFastQs("./my_experiment/fastQs") #use if fastq files end with ".fq.gz"

createCountTableFromFastQs <- function(in_dir,fq_suffix = ".fq.gz") {
  file.list <- list.files(in_dir)
  fastqs <- file.list[endsWith(file.list,fq_suffix)]
  fastqs <- paste(in_dir,fastqs,sep="/")

  content <- ShortRead::readFastq(fastqs[1])
  counts <- ShortRead::tables(content,n=length(content))
  countTable <- as.data.frame(counts$top)
  names(countTable) <- c(basename(fastqs[1]))
  countTable$Read <- row.names(countTable)
  fastqs <- fastqs[-1]
  for(f in fastqs) {
    print(paste("Counting all reads in",f))
    content <- ShortRead::readFastq(f)
    counts <- ShortRead::tables(content,n=length(content))
    counts <- as.data.frame(counts$top)
    names(counts) <- c(basename(f))
    counts$Read <- row.names(counts)
    countTable <- plyr::join(countTable, counts, type="full")
  }
  row.names(countTable) <- countTable$Read
  countTable$Read <- NULL
  countTable[is.na(countTable)] <- 0
  return(countTable)
}
