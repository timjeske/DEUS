% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{runClustering}
\alias{runClustering}
\title{Run CD-HIT clustering}
\usage{
runClustering(cdhit_path, sequences, out_dir, identity_cutoff,
  length_cutoff, wordlength, map, write_fastas = FALSE, optional = "")
}
\arguments{
\item{cdhit_path}{Path to cd-hit-est executable}

\item{sequences}{Vector of sequences in FASTA style generated by the \link[DEUS]{sequencesAsFasta}}

\item{out_dir}{Directory to save output files of clustering}

\item{identity_cutoff}{Sequence identity cutoff used for clustering}

\item{length_cutoff}{Length difference cutoff}

\item{wordlength}{CD-Hit word length}

\item{map}{A data frame with sequences as row names and sequence identifiers in first column.
Can be generated by \link[DEUS]{createMap}}

\item{write_fastas}{Boolean that indicates whether a fasta file will be generated for each cluster}

\item{optional}{Optional execution parameters}
}
\value{
A data frame with the columns 'SequenceID' and 'ClusterID' assigning each sequence to a cluster of similar sequences via their identifiers.
Additionally, a file CD-HIT.fa, CD-HIT.fa.clstr and a folder Clusters is generated in the given output directory.
The CD-HIT.fa file is the FASTA file of all cluster representatives.
The CD-HIT.fa.clustr file lists all identified clusters and the assigned sequence identifiers together with the percentage of overlapping sequence with the cluster representative.
In the Clusters directory there is a FASTA file for each cluster.
}
\description{
Sequentially, each sequency is either assinged to an existing cluster or is classified as a new
cluster representative if no matching cluster can be found.
}
