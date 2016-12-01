
# Note - I cut out sample GRC321-ONT22 (Cytosine Arabinoside 1uM treated) because it appears as outlier currently


set_dir <- "~/MEA_Data/RNA_Pilot/RNAseq/gene_counts/"
set_key <- "~/MEA_Data/RNA_Pilot/RNAseq/gene_counts/sample_info.csv"


## Load packages
library("DESeq2")
library("ggplot2")

## Set directory
setwd(set_dir)


## Prepare for cross-species comparisons; Get human orthologs of rat gene IDs
library(biomaRt)
listMarts(host="www.ensembl.org")
human <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
rat <- useMart("ENSEMBL_MART_ENSEMBL", dataset="rnorvegicus_gene_ensembl", host="www.ensembl.org")

# Find all rat gene orthologs of human gene IDs
hum_rat_ortho <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "rnorvegicus_homolog_ensembl_gene", "rnorvegicus_homolog_associated_gene_name"), mart=human)
hum_rat_ortho <- subset(hum_rat_ortho, rnorvegicus_homolog_associated_gene_name!="") # remove non-homologs

## Separate biomart query to get types of gene annotations
hum_gene_types <- getBM(attributes=c("ensembl_gene_id", "gene_biotype"), mart=human)
hum_pc_genes <- subset(hum_gene_types, gene_biotype=="protein_coding")
rm(hum_gene_types)

# Merge lists, filter to only protein-coding genes
hum_rat_ortho <- merge(hum_rat_ortho, hum_pc_genes, by="ensembl_gene_id") # 21,160 genes
rm(hum_pc_genes)

# Get human entrezIDs to match with Ensembl IDs
#human_entrez <- getBM(attribute=c("ensembl_gene_id", "entrezgene"), mart=human)
#human_entrez <- na.omit(human_entrez) # remove rows without corresponding names
#hum_rat_ortho2 <- merge(hum_rat_ortho, human_entrez, by="ensembl_gene_id") # merge to only entrezgene matching entries





## Read in all similar gene expression data from set directory and assemble into one counts table
read_expression_table <- function(set_dir, pattern="*_geneCounts_ENSEMBL.txt"){
  
  fileList <- list.files(path=set_dir, pattern=pattern)
  
  for (geneCounts in fileList){
    data <- read.delim(geneCounts, header=F, stringsAsFactors=FALSE, skip=2, col.names=c("Geneid","Chr","Start","End","Strand","Length","Counts"))
    data_rd <- data[,c(1,7)] # take gene name and counts columns only
    sample_name <- strsplit(geneCounts, split="-")[[1]][1]

    
    if (exists("all_data")) {
      all_data[,sample_name] <- data_rd[,2]  # append column from this sample with sample name
    
    } else {
      all_data <- data_rd
      names(all_data)[2] <- sample_name
    }
  }  
  all_data
}


## Filter expression counts to exclude information-poor genes
exp_table <- read_expression_table(set_dir)
countData <- exp_table[apply(exp_table[,2:ncol(exp_table)], 1, function(x) sum(x > 10)) / (ncol(exp_table)-1) >= 0.1,] # select genes (rows) which have at least greater than 10 reads for at least 10% of the samples
rownames(countData) <- countData[,1] # set rownames from first column
countData[,1] <- NULL # remove first column from df

## Get sample key
## This should be separate .csv file with experiment summary info like sample name, treatment, concentration, sequencing run, etc.
colData <- read.csv(set_key, row.names=1)
colData[,"Concentration"] <- as.factor(colData[,"Concentration"])
# all(rownames(colData)==colnames(countData))  # to check for consistent sample names; this should be all true

## Create new factor that combines compound and concentration information
colData[,"Group"] <- factor(paste0(colData[,"Treatment"], "_", colData[,"Concentration"])) 


## Initiate DESeq object, comparison formula, and specify control references
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design= ~ Group)
dds$Group <- relevel(dds$Group, ref="Control_0") # this specifies the reference level of treatment, which we labeled as "Control"
dds$Treatment <- relevel(dds$Treatment, ref="Control") 

## Run DESeq algorithm
dds <- DESeq(dds)

## Get results table for each contrast
for (i in unique(colData$Group)) {
  if (i != "Control_0") {
    res <- results(dds, contrast=c("Group", i, "Control_0"))
    print(paste0("Comparing ", i, " vs. ", "Controls"))
    print(summary(res))
    
    write.table(res[order(res$log2FoldChange, decreasing=T),], paste0("../", i, "_vs_control_DESeq.txt"), quote=FALSE, sep="\t")
    
    pdf(paste0("../", i, " vs. ", "controls", ".pdf"))
    plotMA(res, ylim=c(-2,2), main=i, cex.axis=1.4, cex.lab=1.6) # generate MA-plot
    dev.off()
    
    res <- data.frame(res)
    res[,"rnorvegicus_homolog_ensembl_gene"] <- row.names(res)
    res_ids <- merge(res, hum_rat_ortho, by="rnorvegicus_homolog_ensembl_gene") # associate rat Ensembl IDs with human Ensembl IDs
    res_sub <- data.frame(res_ids[,"external_gene_name"], res_ids[,"log2FoldChange"]) 
    names(res_sub) <- c("gene","value")
    res_means <- aggregate(value ~ gene, res_sub, FUN=mean) # aggregate to take the average of multiple fold-changes associated with one gene ID
    res_ord <- res_means[order(res_means$value, decreasing=T),] # order by fold-change only
    write.table(res_ord, paste0("../", i, "_vs_control_GSEA.rnk"), quote=FALSE, sep="\t", row.names=FALSE)  
  }
}



## regularized log transformation
rld <- rlog(dds, blind=TRUE)  # use blind=TRUE for unbiased look at sample clustering


## look at sample distances among transformed values
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$Treatment, rld$Concentration, sep="_")
colnames(sampleDistMatrix) <- paste(rld$Treatment, rld$Concentration, sep="_")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf("../sampleDist_heatmap.pdf")
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
dev.off()


## PCA plot
data <- plotPCA(rld, intgroup=c("Treatment","Concentration"), returnData=TRUE)
percentVar <- round(100*attr(data, "percentVar"))
pdf("../sample_PCA.pdf")
ggplot(data, aes(PC1, PC2, color=Treatment, shape=Concentration)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
dev.off()


## Heatmap of highest variance genes
RowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

pdf("../heatmap_top50_variance_genes.pdf")
heatmap(assay(rld)[head(order(RowVar(assay(rld)), decreasing=TRUE), 50),], margins=c(8,14))
dev.off()


## Other QC plots
pdf("../dispersion_plot.pdf")
plotDispEsts(dds) # dispersion plot
dev.off()

pdf("../Cooks_Dist_boxplots.pdf")
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, main="Cook's distances") # boxplots of Cook's distances per sample
dev.off()






## regularized log transformation for downstream use
#rld <- rlog(dds, blind=FALSE)  # use blind=FALSE for downstream use of normalized data
#norm_genes <- assay(rlog)













