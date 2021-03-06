% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{processClusters}
\alias{processClusters}
\title{Processes CD-HIT result file}
\usage{
processClusters(map, clusters, out_dir, write_fastas)
}
\arguments{
\item{map}{A data frame with sequences as row names and sequence identifiers in first column.
Can be generated by \link[DEUS]{createMap}}

\item{clusters}{The path to CD-HIT.fa}

\item{out_dir}{Directory of CD-HIT result file and where generated files will be saved}

\item{write_fastas}{Boolean that indicates whether a fasta file will be generated for each cluster}
}
\value{
A data frame with the columns 'qseqid', 'cl_id' and 'sequences' containing the sequence identifier, the sequence and the assigned cluster identifier.
}
\description{
This function is internally called by \link[DEUS]{runClustering}
}
