#' Function to execute CD-HIT Clustering
#'
#' @param cdhit_path Path to cd-hit-est executable
#' @param sequences File including the sequences in FASTA format
#' @param out_dir Output folder
#' @param identityCutoff Sequence identity cutoff used for clustering
#' @param lengthCutoff Length difference cutoff
#' @param wordlength CD-Hit word length
#' @param optional Optional execution parameters
#' @keywords clustering
#' @export
#' @examples

# run clustering
runClustering <- function(cdhit_path, sequences, out_dir, identityCutoff, lengthCutoff, wordlength, optional="") {
  outfile <- paste(out_dir,"DESeq2_sig_sequences.cd_hit",sep="/")
  command <- paste("-i",sequences,"-c",identityCutoff,"-s",lengthCutoff,"-n",wordlength,"-g 1",optional,"-o",outfile,sep=" ")
  cluster.out <- system2(cdhit_path,command, stdout=TRUE)

  fileName <- "/storageNGS/ngs4/projects/other/Mouse_beckers_huypens/smallRNA/Pipeline/2017_10_12_checkCluster/DESeq2_sig_sequences.cd_hit.clstr"
  out=readChar(fileName, file.info(fileName)$size)
  li=strsplit(out,">Cluster ")
  li[[1]]=li[[1]][-1]
  d=data.frame(unlist(li))
  outdir=paste(out_dir, "Clusters",sep="/")
  dir.create(outdir, showWarnings = FALSE)
  apply(d,1,extractSequences,map=map)
}





extractSequences<-function(entry,map){
  cl_id<-unlist(strsplit(entry,"\n"))[1]
  outfile=paste(outdir,cl_id,sep="/")
  splitentry=unlist(strsplit(entry,","))
  m=regexpr("seq_[0-9]+",splitentry,perl=TRUE)
  o=c(regmatches(splitentry, m))
  seq=row.names(map)[map[,1]%in%o]
  out=paste(">",o,"\n",seq,sep="")
  write.table(out,outfile,row.names=F,col.names=F,quote=F)
}
