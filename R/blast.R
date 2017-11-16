#' Function to blast significant sequences
#'
#' @param blast_exec path to your installation of 'blastn'
#' @param blast_db path to your local blastdb
#' @param ncores number of cores used for blasting
#' @param sigResults data frame with significant sequences as row names
#' @param map data frame with sequences as row names and some identifier for each sequence in the first column
#' @param outdir output folder
#' @keywords blastn
#' @export
#' @examples

runBlast <- function(blast_exec, blast_db, ncores, sigResults, map, out_dir){
  sig_sequences <- sequencesAsFasta(sigResults,map)
  blast.f6 <- c('qseqid', 'qlen', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart' ,'send', 'evalue', 'bitscore')
  command=paste("-db",blast_db, "-num_threads", ncores , "-perc_identity 100 -strand plus -task blastn-short -qcov_hsp_perc 100 -max_target_seqs 2000 -outfmt", sprintf(" '6 %s'",paste(collapse=" ",blast.f6)),sep=" ")
  blast.out <- system2(blast_exec,command, input=sig_sequences, stdout=TRUE)
  blast.out.df <- `names<-`(read.table(quote="",sep='\t',textConnection(blast.out)),blast.f6)
  return(blast.out.df)
}
