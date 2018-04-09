#' Function to run differential expression analysis based on DESeq2
#'
#' @param countData raw counts of sequences
#' @param phenoData data frame with sample names as rownames and assigned condition in first column
#' @param design design formular for DE analysis
#' @param map data frame with sequences as rownames and sequence ID in first column
#' @param out_dir directory for sample distance map and MA plot
#' @keywords differential expression analysis
#' @export
#' @examples

runDESeq2 <- function(countData, phenoData, design, map, out_dir) {
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = phenoData, design = design)
  dds <- dds[rowMeans(DESeq2::counts(dds)) > 1, ]
  dds <- DESeq2::DESeq(dds, betaPrior=TRUE)
  res <- DESeq2::results(dds)
  res <- res[, c(ncol(res)-1, ncol(res), (1:(ncol(res)-2)))]
  res$"IHWPval"=IHW::adj_pvalues(IHW::ihw(res$pvalue~res$baseMean, data=res,alpha=0.1))

  # plot sample distance
  rld <- DESeq2::rlog(dds, blind=FALSE)
  plotSampleDistanceMap(rld,out_dir)

  # plot PCA
  plotPCA(rld,phenoData,out_dir)

  # plot MA
  pdf(paste(out_dir,"DESeq2_MAplot_shrunken.pdf",sep="/"), onefile=FALSE)
  DESeq2::plotMA(res, main = "DESeq2", ylim=c(-4,4))
  dev.off()

  counts_norm<-DESeq2::counts(dds, normalized=TRUE)
  newList <- list("normCounts" = counts_norm, "deResult" = as.data.frame(res))
  return(newList)
}

#' Function to plot sample distance map
#'
#' @param rld regularized log transformed counts
#' @param out_dir directory for sample distance map
#' @keywords sample distance map
#' @export
#' @examples

plotSampleDistanceMap <- function(rld, out_dir) {
  pdf(paste(out_dir, "DESeq2_sample_dist.pdf", sep="/"), onefile=FALSE)
  sampleDists <- dist(t(SummarizedExperiment::assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(rld)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  dev.off()
}

#' Function to plot PCA
#'
#' @param rld regularized log transformed counts
#' @param phenoData data frame with sample names as rownames and assigned condition in first column
#' @param out_dir directory for sample distance map
#' @keywords sample distance map
#' @export
#' @examples

plotPCA <- function(rld, phenoData, out_dir) {
  myPhenoData <- phenoData[, !(names(phenoData) %in% c("sample")), drop=FALSE]
  pdf(paste(out_dir,"DESeq2_PCA.pdf", sep="/"), onefile=FALSE)
  if(ncol(myPhenoData) > 1) {
    data <- plotPCA(rld, intgroup=colnames(myPhenoData), returnData=TRUE)
    percentVar <- round(100 * attr(data, "percentVar"))
    plt <- ggplot(data, aes(PC1, PC2, color=myPhenoData[,1], shape=myPhenoData[,2])) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      labs(color=colnames(myPhenoData)[1]) +
      labs(shape=colnames(myPhenoData)[2])
  } else {
    plt <- plotPCA(rld, intgroup=colnames(myPhenoData))
  }
  print(plt)
  dev.off()
}
