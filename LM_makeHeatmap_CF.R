
# Load required packages
library(colorRamps)
library(gplots)



# Copy LM data from excel
data <- read.table("clipboard", sep="\t", header=TRUE, check.names=FALSE)

# Rename rows
row.names(data) <- c("MFR","BR","ISI","%SiB","BD","IBI","#AE","#ABE","#NS","NSP","NSD","%SiNS","NS-ISI","NSDsd","#SiNS","r")

# Merge with "q2" data frame from AUC-based EC50 heatmap construction
q_lm <- merge(q2[,"seq", drop=FALSE], t(data), by="row.names")
q_lm <- q_lm[order(q_lm$seq),]
q_lm <- q_lm[,!"seq"]
row.names(q_lm) <- q_lm[,"Row.names"]
q_lm[,"seq"] <- NULL
q_lm[,"Row.names"] <- NULL





# Heatmap with binary (significant vs. not) colors
ramp<-colorRampPalette(c("red","gray20"),space="rgb")

data <- q_lm
data[data < .05] <- 0
data[data > .05] <- 1

pdf("linearmodel_FDR.05Sig_heatmap.pdf", width=11)
heatmap.2(t(data), Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=ramp(200), margins=c(16,14), trace="none", density.info="none", colsep=seq(1:dim(data)[1]), rowsep=seq(1:dim(data)[1]), symkey=FALSE, keysize=1, cexCol=1, cexRow=1.1)
dev.off()


## Note - need to reload dataset to change these values again
data <- q_lm
data[data < .01] <- 0
data[data > .01] <- 1

pdf("linearmodel_FDR.01Sig_heatmap.pdf", width=11)
heatmap.2(t(data), Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=ramp(200), margins=c(16,14), trace="none", density.info="none", colsep=seq(1:dim(data)[1]), rowsep=seq(1:dim(data)[1]), symkey=FALSE, keysize=1, cexCol=1, cexRow=1.1)
dev.off()







# Heatmap with range of colors corresponding to log(p-value)
ramp<-colorRampPalette(c("cyan","cornflowerblue","gray20"),space="rgb")
col_breaks = c(seq(-50,-10.01,length=100), # for red
               seq(-10,-3,length=100), # for purple
               seq(-2.99,0,length=1)) # for black

pdf("linearmodel_range_heatmap.pdf", width=11)
heatmap.2(log(as.matrix(data)), Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=ramp(200), margins=c(16,14), trace="none", density.info="none", key=TRUE, colsep=seq(1:length(unique(names(data)))), rowsep=seq(1:length(unique(names(data)))), breaks=col_breaks, symkey=FALSE, key.xlab="log(p-value)", keysize=1)
dev.off()

