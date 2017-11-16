#' Function to execute CD-HIT Clustering
#'
#' @param cdhit_path Path to cd-hit-est executable
#' @param sequences File including the sequences in FASTA format
#' @param out_dir Output folder
#' @param identityCutoff Sequence identity cutoff used for clustering
#' @param lengthCutoff Length difference cutoff
#' @param wordlength CD-Hit word length
#' @param map data frame with sequences as row names and some identifier for each sequence in the first column
#' @param optional Optional execution parameters
#' @keywords clustering
#' @export
#' @examples
runClustering <- function(cdhit_path, sequences, out_dir, identityCutoff, lengthCutoff, wordlength, map, optional="") {
  clusters <- paste(out_dir,"DESeq2_sig_sequences.cd_hit",sep="/")
  command <- paste("-i",sequences,"-c",identityCutoff,"-s",lengthCutoff,"-n",wordlength,"-g 1",optional,"-o",clusters,sep=" ")
  cluster.out <- system2(cdhit_path,command, stdout=TRUE)
  cl.dir<-paste(out_dir, "Clusters",sep="/")
  dir.create(cl.dir, showWarnings = FALSE)
  out.df <- processClusters(map,out_dir)
  out.df["sequences"] <- NULL
  names(out.df)<-c("SequenceID","ClusterID")
  return(out.df)
}

#' Processes CD-Hit Outfile
#'
#' @param map data frame with sequences as row names and some identifier for each sequence in the first column
#' @param out_dir Output folder
#' @keywords clustering
#' @export
#' @examples
processClusters <- function(map,out_dir) {
  cl.dir <- paste(out_dir, "Clusters",sep="/")
  clusters <- paste(out_dir,"DESeq2_sig_sequences.cd_hit",sep="/")
  cluster.file <- paste(clusters,".clstr",sep="")
  con <- file(cluster.file, "r")
  out.df <- data.frame()

  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if(startsWith(line,">")){
      cl_id <- unlist(strsplit(line,"\\s+",perl = T))[2]
    }else{
      seqIDsIndex <- regexpr("seq_[0-9]+",line,perl=TRUE)
      qseqid <- c(regmatches(line, seqIDsIndex))
      out.df <- rbind.data.frame(out.df,cbind(qseqid,cl_id))
    }
  }
  close(con)

  map$sequences<-row.names(map)
  out.df<-merge(x = out.df, y = map, by.x = "qseqid",by.y="seq_id", all.x = TRUE)
  by(out.df,out.df$cl_id,printCluster,cl.dir=cl.dir)
  return(out.df)
}

#' Helper function for generating individual cluster fastas, not meant for external use
#'
#' @param data Subset of overall data provided within the processClusters function
#' @param cl.dir Cluster output folder
#' @keywords clustering
#' @export
#' @examples
printCluster<-function(data,cl.dir){
  cl_id<-data[1,2]
  outfile <- paste(cl.dir,"/Cluster",cl_id,".fa",sep="")
  out<-paste(">",data[,1],"\n",data[,3],sep="")
  write.table(out,outfile,col.names=F,row.names=F,quote=F)
}
