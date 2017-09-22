#' Function to blast significant sequences
#'
#' @param blast_exec
#' @param blast_db
#' @param ncores
#' @param sigResults
#' @param map
#' @keywords filtering
#' @export
#' @examples
#' runBlast(...)

runBlast <- function(blast_exec, blast_db, ncores, sigResults, map){
  sig_sequences <- as.vector(rbind(paste(">",map[row.names(sigResults),1],sep=""),row.names(sigResults)))
  blast.f6 <- c('qseqid', 'qlen', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart' ,'send', 'evalue', 'bitscore')
  command=paste("-db",blast_db, "-num_threads", ncores , "-perc_identity 100 -strand plus -task blastn-short -qcov_hsp_perc 100 -max_target_seqs 2000 -outfmt", sprintf(" '6 %s'",paste(collapse=" ",blast.f6)),sep=" ")
  blast.out <- system2(blast_exec,command, input=sig_sequences, stdout=TRUE)
  blast.out.df <- `names<-`(read.table(quote="",sep='\t',textConnection(blast.out)),blast.f6)
  return(blast.out.df)
}
