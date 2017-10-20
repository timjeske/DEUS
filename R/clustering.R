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

# run clustering
runClustering <- function(cdhit_path, sequences, out_dir, identityCutoff, lengthCutoff, wordlength, map, optional="") {
  clusters <- paste(out_dir,"DESeq2_sig_sequences.cd_hit",sep="/")
  command <- paste("-i",sequences,"-c",identityCutoff,"-s",lengthCutoff,"-n",wordlength,"-g 1",optional,"-o",clusters,sep=" ")
  cluster.out <- system2(cdhit_path,command, stdout=TRUE)

  #Read in outfile and prepare for further processing
#  cluster.file <- paste(clusters,".clstr",sep="")
#  out<-readChar(cluster.file, file.info(cluster.file)$size)
#  cl.list<-strsplit(out,">Cluster ")
#  #Remove first empty element
#  cl.list[[1]]<-cl.list[[1]][-1]
#  cl.df<-data.frame(unlist(cl.list))
  cl.dir<-paste(out_dir, "Clusters",sep="/")
  dir.create(cl.dir, showWarnings = FALSE)
  out.df <- processClusters(map,out_dir)

#  out.df<-(apply(cl.df,1,extractSequences,map=map,cl.dir=cl.dir))
#  out.df<-do.call("rbind",out.df)
  out.df["sequences"] <- NULL
  names(out.df)[2]<-"ClusterID"
  return(out.df)
}

extractSequences<-function(entry,map,cl.dir){
  cl_id<-unlist(strsplit(entry,"\n"))[1]
  outfile<-paste(cl.dir,cl_id,sep="/")
  splitentry<-unlist(strsplit(entry,","))
  seqIDsIndex<-regexpr("seq_[0-9]+",splitentry,perl=TRUE)
  qseqid<-c(regmatches(splitentry, seqIDsIndex))
  sequences<-row.names(map)[map[,1]%in%qseqid]
  out=paste(">",qseqid,"\n",sequences,sep="")
  write.table(out,outfile,row.names=F,col.names=F,quote=F)
  out.df<-as.data.frame(cbind(qseqid,sequences,cl_id))
  return(out.df)
}

processClusters <- function(map,out_dir) {
  cl.dir <- paste(out_dir, "Clusters",sep="/")
  clusters <- paste(out_dir,"DESeq2_sig_sequences.cd_hit",sep="/")
  cluster.file <- paste(clusters,".clstr",sep="")
  con <- file(cluster.file, "r")
  outfile <- -1
  out.df <- data.frame()

  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if(startsWith(line,">")){
      if(!(outfile==-1)){
        close(out.con)
      }
      cl_id <- unlist(strsplit(line,"\\s+",perl = T))[2]
      outfile <- paste(cl.dir,"/Cluster",cl_id,".fa",sep="")
      out.con <- file(outfile,"w")
    }else{
      seqIDsIndex <- regexpr("seq_[0-9]+",line,perl=TRUE)
      qseqid <- c(regmatches(line, seqIDsIndex))
      sequences <- row.names(map)[map[,1]%in%qseqid]
      write(paste(">",qseqid,sep=""),out.con)
      write(sequences,out.con)
      out.df <- rbind.data.frame(out.df,cbind(qseqid,sequences,cl_id))
    }
  }
  close(con)
  close(out.con)
  return(out.df)
}

