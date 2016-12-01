library(ggplot2)
library(gplots)


setwd("~/MEA_Data/PIP3_AIM2-3-DNT/TippingPointOutputs/ReplicatePoints_Span0.80_7-11-16Update/")
data <-  read.table("clipboard", sep="\t", header=TRUE) # Read in tipping point, cytotox, and AUC EC50 concentrations

# Eliminate unused parameters
#data <- subset(data, !(ontogeny %in% c("mean.isis_AUC", "mean.dur_AUC", "mean.IBIs_AUC", "ns.peak.m_AUC", "ns.durn.m_AUC", "ns.mean.insis_AUC", "ns.durn.sd_AUC", "ns.mean.spikes.in.ns_AUC","cv.time_AUC","cv.network_AUC")))


data_ord <- data[order(data[,"MinCritConc"]),] # Order by ascending critical concentration
data_ord[,"Compound"] <- factor(data_ord[,"Compound"], levels=rev(data_ord[,"Compound"])) # specify the factor order to prevent ggplot from reordering compounds

# If we want to manually calculate confidence intervals
#data_ord[,"EstError"] <- apply(data.frame(a=(data_ord[,"ConcAbove"] - data_ord[,"MinCritConc"]), b=(data_ord[,"MinCritConc"] - data_ord[,"ConcBelow"])), 1, FUN=max)
#data_ord[,"LowerError"] <- data_ord[,"MinCritConc"] - data_ord[,"EstError"]
#data_ord[,"UpperError"] <- data_ord[,"MinCritConc"] + data_ord[,"EstError"]


data_ord[is.na(data_ord)] <- 100 # Optional - set all undetermined values to 100 uM


# Generate comparison plot vs. minimums
vis <- ggplot(data=data_ord, aes(y=Compound, x=MinCritConc)) + 
geom_point(size=2) + 
geom_point(aes(x=CytotoxMin), col="orange", size=2.5, shape=17) +
geom_point(aes(x=MinAUC), col="blue", size=2, shape=15) +
geom_errorbarh(aes(y=Compound, xmax=CI_hi, xmin=CI_low), alpha=0.5, height=0.2) +
scale_x_log10(name="Concentration", breaks=c(0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)) +
ggtitle("Critical Conc. vs. Cytotox and LEC")

pdf("CritConc_vs_CytotoxandLEC_comparison.pdf", width=8)
vis
dev.off()



# Generate comparison plot vs. AUC median conc.
data_ord[,"MedAUC"] <- apply(data_ord[,6:13], 1, function(x) median(x, na.rm=TRUE)) # calculate median of AUC values, ignoring NAs
vis <- ggplot(data=data_ord, aes(y=Compound, x=MinCritConc)) + 
  geom_point(size=2) + 
  geom_point(aes(x=CytotoxMin), col="orange", size=2.5, shape=17) +
  geom_point(aes(x=MedAUC), col="blue", size=2, shape=15) +
  geom_point(aes(x=MinAUC), col="purple", size=2, shape=15) +
  geom_errorbarh(aes(y=Compound, xmax=CI_hi, xmin=CI_low), alpha=0.5, height=0.2) +
  scale_x_log10(name="Concentration", breaks=c(0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)) +
  ggtitle("Critical Conc. vs. Cytotox and MedianC")

pdf("CritConc_vs_CytotoxandAUCMedian_comparison.pdf", width=8)
vis
dev.off()



# Generate comparison plot vs. both AUC median and minimum
vis <- ggplot(data=data_ord, aes(y=Compound, x=MinCritConc)) + 
  geom_point(size=2) + 
  geom_point(aes(x=MedAUC), col="blue", size=2, shape=15) +
  geom_point(aes(x=MinAUC), col="purple", size=2, shape=15) +
  geom_point(aes(x=CytotoxMin), col="orange", size=2.5, shape=17) +
  geom_errorbarh(aes(y=Compound, xmax=CI_hi, xmin=CI_low), alpha=0.5, height=0.2) +
  scale_x_log10(name="Concentration", breaks=c(0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)) +
  ggtitle("Critical Conc. vs. Cytotox and MedianC")

pdf("CritConc_vs_DoubleAUC_comparison.pdf", width=8)
vis
dev.off()






## Generate heatmap of endpoints used
# First subset data to just AUC endpoints and set NA to 100
q <- data_ord[,6:13]
row.names(q) <- data_ord[,"Compound"]
q[is.na(q)] <- 100

ramp<-colorRampPalette(c("cyan","cornflowerblue","gray20"),space="rgb")
col_breaks = c(seq(-5,1,length=100), 
               seq(1.01,4.5,length=100), 
               seq(4.51,5,length=1))

pdf("NetworkParamAUC_Endpoints.pdf", width=6)
heatmap.2(log(as.matrix(q)), Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=ramp(200), margins=c(16,18), trace="none", density.info="none", key=TRUE, colsep=seq(1:8), rowsep=seq(1:35), breaks=col_breaks, symkey=FALSE, key.xlab="log(EC50)", keysize=1)
dev.off()









