#!/usr/bin/Rscript

# call: Rscript DESeq2.R counts.tsv conditions.tsv "~ condition" "(pvalue|padj|IHWPval)" 0.05 map out_dir [gene_list]
# or set: countData, colData, design, filter_val, filter_threshold and map
# map: first col is entity id, second col is description for plots

# CRAN packages
for (package in c("scales","ggplot2","pheatmap","RColorBrewer","tools")) {
    if (!require(package, character.only=T, quietly=T)) {
        install.packages(package, repos = "http://cran.us.r-project.org")
    }
    suppressMessages(library(package, character.only=T))
}

# Bioconductor packages
source("http://bioconductor.org/biocLite.R")

for (package in c("DESeq2","IHW")) {
    if (!require(package, character.only=T, quietly=T)) {
        biocLite(package)
    }
    suppressMessages(library(package, character.only=T))
}

# get mean and sd for normalized counts for each condition
getCountStats<-function(countData,colData){

        pheno <- colData
        #Get Group Means
        groups = unique(pheno$condition)
        for(type in groups){
                cols  = row.names(pheno)[which(pheno$condition==type)]
                subset= countData[,cols]
                countData=cbind(countData,rowMeans(subset),data.frame(apply(subset,1,sd)))
                names(countData)[ncol(countData)-1]=paste("NormCounts",type,"Mean",sep="_")
                names(countData)[ncol(countData)]=paste("NormCounts",type,"Sd",sep="_")
        }
        index=length(groups)*2
        #Return only Mean & Sd columns
        return(countData[,c((ncol(countData)-index+1):ncol(countData))])
}


getMatrix<-function(all_list,DE_list,gene_list) {
	overlap_DE <- length(intersect(DE_list,gene_list))

        DE_IN <- overlap_DE
        DE_nIN <- length(DE_list) - overlap_DE
        nDE_IN <- length(intersect(gene_list,all_list)) - overlap_DE
        nDE_nIN <- length(all_list) - DE_IN - DE_nIN - nDE_IN
        x <- matrix(c(DE_IN, DE_nIN, nDE_IN, nDE_nIN),2,2)
        dimnames(x) <- list( c("in list", "not in list"), c("DE", "not DE"))
        return(x)
}


args<-commandArgs(TRUE)
if( length(args) > 6 ) {
	countData <- read.table(args[1], header=T, row.names=1, check.names=FALSE)
	colData <- read.table(args[2], header=T, row.names=1, check.names=FALSE)
	design <- as.formula(args[3])
	filter_val <- args[4]
	filter_threshold <- as.numeric(args[5])
	map <- read.table(args[6], row.names=1)
	out_dir <- args[7]
}

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = design)
dds <- dds[rowMeans(counts(dds)) > 1, ]
dds <- DESeq(dds, betaPrior=TRUE)
df <- as.data.frame(colData(dds)[,c("condition")])
names(df) <- c("condition")
row.names(df) <- row.names(colData)

#Get normalized read counts
dds_norm <- estimateSizeFactors(dds)
counts_norm<-counts(dds_norm, normalized=TRUE)
write.table(counts_norm, paste(out_dir, "DESeq2_counts_norm.tsv", sep="/") , sep="\t", quote=F)
dds_norm<-0

#regularized log transformation
rld <- rlog(dds, blind=FALSE)

#most expressed heatmap
pdf(paste(out_dir, "DESeq2_most_expressed_heatmap.pdf", sep="/") ,onefile=FALSE)
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:50]
tmp <- assay(rld)[select,]
row.names(tmp) <- as.vector(map[row.names(tmp),1])
pheatmap(tmp,  annotation_col=df, fontsize_row=7)
dev.off()

#sample-to-sample distance
pdf(paste(out_dir, "DESeq2_sample_dist.pdf", sep="/"), onefile=FALSE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(rld)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
clustering_distance_rows=sampleDists,
clustering_distance_cols=sampleDists,
col=colors)
dev.off()

res <- results(dds)

#Run IHW Method
res$"IHWPval"=adj_pvalues(ihw(res$pvalue~res$baseMean, data=res,alpha=0.1))

print_res <- cbind(res,getCountStats(counts_norm,colData))
write.table(print_res[order(res$pvalue),], paste(out_dir,"DESeq2_output.tsv",sep="/"), sep="\t", quote=F)
sig_res <- print_res[!is.na(print_res[[filter_val]]) & print_res[[filter_val]] < filter_threshold,]
sig_res <- sig_res[order(sig_res$pvalue),]
write.table(sig_res, paste(out_dir,"DESeq2_sig_output.tsv",sep="/"), sep="\t", quote=F)
write.table(sig_res[sig_res$log2FoldChange < 0, ], paste(out_dir,"DESeq2_sig_down.tsv",sep="/"), sep="\t", quote=F)
write.table(sig_res[sig_res$log2FoldChange > 0, ], paste(out_dir,"DESeq2_sig_up.tsv",sep="/"), sep="\t", quote=F)

#PCA-plot
pdf(paste(out_dir,"DESeq2_PCA.pdf", sep="/"), onefile=FALSE)
if(ncol(colData) > 1) {
	data <- plotPCA(rld, intgroup=colnames(colData), returnData=TRUE)
	percentVar <- round(100 * attr(data, "percentVar"))
	plt <- ggplot(data, aes(PC1, PC2, color=colData[,1], shape=colData[,2])) +
	geom_point(size=3) +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance"))
} else {
	plt <- plotPCA(rld, intgroup=colnames(colData))
}
print(plt)
dev.off()

#MA-plot
pdf(paste(out_dir,"DESeq2_MAplot_shrunken.pdf",sep="/"), onefile=FALSE)
DESeq2::plotMA(res, main = "DESeq2")
dev.off()

#heatmap of significance
n <- sum(!is.na(res[[filter_val]]) & res[[filter_val]]<filter_threshold)
if(n > 1) {
	if(n > 50) {
		n <- 50
	}
	pdf(paste(out_dir,"DESeq2_significant_heatmap.pdf",sep="/"), onefile=FALSE)
	select <- order(res$pvalue)[1:n]
	tmp <- assay(rld)[select,]
	row.names(tmp) <- as.vector(map[row.names(tmp),1])
	pheatmap(tmp,  annotation_col=df, fontsize_row=7)
	dev.off()

	pdf(paste(out_dir,"DESeq2_significant_log_z-scores_heatmap.pdf",sep="/"), onefile=FALSE)
        mycounts <- assay(rld)[select,]
	row.names(mycounts) <- as.vector(map[row.names(mycounts),1])
        row_sd <- apply(mycounts,1,sd)
        z_matrix <- (mycounts-rowMeans(mycounts))/row_sd
        pheatmap(z_matrix, annotation_col=df, fontsize_row=7)
        dev.off()
}

#additional plots if gene list is given as fourth input

if(length(args)>7) {

	geneLists <- strsplit(args[8], ",")[[1]]

	for(geneList in geneLists) {
		print(geneList)
		my_name <- file_path_sans_ext(basename(geneList))

		#special for InduRag analysis, not harmful for other gene lists
		my_name <- strsplit(my_name, "_SUM4.genes")[[1]]
		dir.create(file.path(getwd(),my_name))
		sub_out_dir <- paste(out_dir, my_name, sep="/")
        	geneList <- read.table(geneList, header=F, row.names=1)
       		overlap_ids <- length(intersect(row.names(counts_norm),row.names(geneList)))
		overlap_names <- length(intersect(as.vector(map[row.names(counts_norm),1]),row.names(geneList)))

		overlap <- overlap_ids
		gene_ids <- row.names(geneList)
		if(overlap_ids < overlap_names) {
			overlap <- overlap_names
			gene_ids <- row.names(map)[map[,1] %in% row.names(geneList)]
		}

		#print enrichment stuff
		sink(paste(sub_out_dir, paste("DESeq2_enrichment_",my_name,".txt",sep=""),sep="/"))
		cat(sprintf("Length of your list\t%d", length(row.names(geneList))),"\n")
		cat(sprintf("Overlap with count data\t%d", overlap),"\n")
		cat(sprintf("DE genes\t%d", nrow(sig_res)),"\n")

		genes_down <- row.names(sig_res[sig_res$log2FoldChange < 0, ])
		genes_up <- row.names(sig_res[sig_res$log2FoldChange > 0, ])
		cat(sprintf("DE genes down\t%d", length(genes_down)), "\n")
		cat(sprintf("DE genes up\t%d", length(genes_up)), "\n")

		overlap_DE <- length(intersect(row.names(sig_res),gene_ids))
		cat(sprintf("List overlap with DE genes\t%d", overlap_DE),"\n")

		#cat("2x2 table all genes","\n")
		#x <- getMatrix(row.names(counts_norm),row.names(sig_res),gene_ids)
		#print(x)
		#print(fisher.test(x))

		cat("","\n")
		cat("2x2 table down-regulated genes","\n")
		x <- getMatrix(row.names(counts_norm),genes_down,gene_ids)
		print(x)
		fres <- fisher.test(x)
		cat(sprintf("p-value\t%.4g", fres$p.value),"\n")
		cat(sprintf("OR\t%.2f", as.numeric(fres$estimate)),"\n")

		cat("","\n")
		cat("2x2 table up-regulated genes","\n")
        	x <- getMatrix(row.names(counts_norm),genes_up,gene_ids)
        	print(x)
        	fres <- fisher.test(x)
		cat(sprintf("p-value\t%.4g", fres$p.value),"\n")
                cat(sprintf("OR\t%.2f", as.numeric(fres$estimate)),"\n")
		sink()

        	#normalized gene counts
        	select <- row.names(counts_norm) %in% intersect(row.names(counts_norm), gene_ids)
		tmp <- counts_norm[select,]
		row.names(tmp) <- as.vector(map[row.names(tmp),1])
        	write.table(tmp, paste(sub_out_dir,paste("DESeq2_all_counts_",my_name,".tsv",sep=""), sep="/"), dec="," , sep="\t", quote=F)

		ov_down <- intersect(gene_ids, genes_down)
		select <- row.names(counts_norm) %in% intersect(row.names(counts_norm), ov_down)
                tmp <- counts_norm[select,]
                row.names(tmp) <- as.vector(map[row.names(tmp),1])
                write.table(tmp, paste(sub_out_dir,paste("DESeq2_down_counts_",my_name,".tsv",sep=""), sep="/"), dec=",", sep="\t", quote=F)

		ov_up <- intersect(gene_ids, genes_up)
		select <- row.names(counts_norm) %in% intersect(row.names(counts_norm), ov_up)
                tmp <- counts_norm[select,]
                row.names(tmp) <- as.vector(map[row.names(tmp),1])
                write.table(tmp, paste(sub_out_dir,paste("DESeq2_up_counts_",my_name,".tsv",sep=""), sep="/"), dec=",", sep="\t", quote=F)

		#MAplot
		select <- row.names(res) %in% intersect(row.names(res), gene_ids)
		pdf(paste(sub_out_dir, paste("DESeq2_MAplot_",my_name,".pdf",sep=""), sep="/"),onefile=FALSE)
		DESeq2::plotMA(res[select,], main="DESeq2")
        	dev.off()

	}

	#heatmap of genes
        #pdf(paste(out_dir, "DESeq2_gene_list_heatmap.pdf", sep="/") ,onefile=FALSE)
        #select <- row.names(dds) %in% intersect(row.names(dds), gene_ids)
	#tmp <- assay(rld)[select,]
	##row.names(tmp) <- as.vector(map[row.names(tmp),1])
        #pheatmap(tmp,  annotation_col=df, fontsize_row=5)
        #dev.off()

        #pdf(paste(out_dir, "DESeq2_gene_list_z-scores_heatmap.pdf", sep="/"), onefile=FALSE)
        #mycounts <- assay(rld)[select,]
	#row.names(mycounts) <- as.vector(map[row.names(mycounts),1])
        #row_sd <- apply(mycounts,1,sd)
        #z_matrix <- (mycounts-rowMeans(mycounts))/row_sd
        #pheatmap(z_matrix, annotation_col=df, fontsize_row=5)
        #dev.off()
}

#res <- cbind(res,getCountStats(counts_norm,colData))
#write.table(res[order(res$pvalue),], paste(out_dir,"DESeq2_output.tsv",sep="/"), sep="\t", quote=F)
#tmp_res <- res[!is.na(res[[filter_val]]) & res[[filter_val]] < filter_threshold,]
#tmp_res <- tmp_res[order(tmp_res$pvalue),]
#write.table(tmp_res, paste(out_dir,"DESeq2_sig_output.tsv",sep="/"), sep="\t", quote=F)
#write.table(tmp_res[tmp_res$log2FoldChange < 0, ], paste(out_dir,"DESeq2_sig_down.tsv",sep="/"), sep="\t", quote=F)
#write.table(tmp_res[tmp_res$log2FoldChange > 0, ], paste(out_dir,"DESeq2_sig_up.tsv",sep="/"), sep="\t", quote=F)

