# AUC summary statistics and plotting
# Chris Frank - March 2016


# Set parameters
set_directory <- "~/MEA_Data/PIP3_AIM2-3-DNT/TEST3/"
set_ec50_table <- "ec50_allOntogeny_withCytotox.txt" # Table of EC50 values for network parameters and cytotox assays
set_output <- "selectivity"


# Load required packages
library(colorRamps)
library(gplots)
library(car)



## Read in ec50 table and rearrange it
setwd(set_directory)
all_ec50 <- read.delim(set_ec50_table,sep="\t", stringsAsFactors=FALSE)

# Remove cv.time, cv.network parameters
all_ec50 <- subset(all_ec50, ontogeny!="cv.time_AUC" & ontogeny!="cv.network_AUC")

# Split by compound
sp <- split(all_ec50, all_ec50[,"compound"])

# Take the column that has ec50 values only and make new data.frame that can be coerced to numerical matrix
q <- data.frame()
for (i in 1:length(sp)) {
  q <- rbind(q, sp[[i]][,3])
}
row.names(q) <- names(sp)
names(q) <- unique(all_ec50[,"ontogeny"])

names(q) <- c("MFR","BR","ISI","%SiB","BD","IBI","#AE","#ABE","#NS","NSP","NSD","%SiNS","NS-ISI","NSDsd","#SiNS","r","ABcytotox","LDHcytotox") 










if (set_output == "std_heatmap") {
  # Set all NAs to 100 for dark gray color in heatmap
  q[is.na(q)] <- 100
  
  # set the color scheme
  ramp<-colorRampPalette(c("red","purple","gray20"),space="rgb")
  col_breaks = c(seq(-5,1,length=100), # for red
                 seq(1.01,4.5,length=100), # for purple
                 seq(4.51,5,length=1)) # for black
  
  
  
  pdf("ec50_heatmap.pdf", width=11)
  heatmap.2(log(t(q)), Colv=TRUE, Rowv=FALSE, scale="none", dendrogram="row", col=ramp(200), margins=c(16,14), trace="none", density.info="none", key=TRUE, colsep=seq(1:length(unique(names(sp)))), rowsep=seq(1:length(unique(names(sp)))), breaks=col_breaks, symkey=FALSE, key.xlab="log(EC50)", keysize=1)
  dev.off()
  
  png("ec50_heatmap.png", width=4000, height=4000, pointsize=90)
  heatmap.2(log(t(q)), Colv=TRUE, Rowv=FALSE, scale="none", dendrogram="row", col=ramp(200), margins=c(16,14), trace="none", density.info="none", key=TRUE, colsep=seq(1:length(unique(names(sp)))), rowsep=seq(1:length(unique(names(sp)))), breaks=col_breaks, symkey=FALSE, key.xlab="log(EC50)", keysize=1.75)
  dev.off()

  
  # try a different color scheme
  ramp<-colorRampPalette(c("cyan","cornflowerblue","gray20"),space="rgb")
  col_breaks = c(seq(-5,1,length=100), 
                 seq(1.01,4.5,length=100), 
                 seq(4.51,5,length=1))
  
  pdf("ec50_heatmap_withCytotox_altColor1.pdf", width=11)
  heatmap.2(log(t(q)), Colv=TRUE, Rowv=FALSE, scale="none", dendrogram="col", col=ramp(200), margins=c(16,14), trace="none", density.info="none", key=TRUE, colsep=seq(1:length(unique(names(sp)))), rowsep=seq(1:length(unique(names(sp)))), breaks=col_breaks, symkey=FALSE, key.xlab="log(EC50)", keysize=1)
  dev.off()
  
  ## Adding a column color code
  # Import excel formatted data on DNT evidence
  dnt <- read.table("clipboard", sep="\t", header=FALSE, row.names=1, strip.white=TRUE)
  names(dnt) <- c("evidence")
  
  # Map values to colors
  dnt_cols <- recode(dnt[,"evidence"], "0='forestgreen'; 1='gainsboro'; 2='firebrick3'; 3='firebrick3'; 4='darkorange'")
  
  pdf("ec50_heatmap_withCytotox_DNTevidence.pdf", width=11)
  heatmap.2(log(t(q)), Colv=TRUE, Rowv=FALSE, scale="none", dendrogram="col", col=ramp(200), margins=c(16,14), trace="none", density.info="none", key=TRUE, colsep=seq(1:length(unique(names(sp)))), rowsep=seq(1:length(unique(names(sp)))), breaks=col_breaks, symkey=FALSE, key.xlab="log(EC50)", keysize=1, ColSideColors=dnt_cols, cexCol=1.1, cexRow=1.1)
  dev.off()
  
}




if (set_output == "cytotox_ordered_heatmap") {
  ## Order heatmap by cytotoxicity minimum EC50 value
  q[,"min"] <- apply(q[,17:18], 1, min)
  q <- q[order(q[,"min"]),]
  
  pdf("ec50_heatmap_CytotoxOrdered.pdf", width=11)
  heatmap.2(log(t(q)), Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=ramp(200), margins=c(16,14), trace="none", density.info="none", key=TRUE, colsep=seq(1:length(unique(names(sp)))), rowsep=seq(1:length(unique(names(sp)))), breaks=col_breaks, symkey=FALSE, key.xlab="log(EC50)", keysize=1)
  heatmap.2(log(t(q[27:70,])), Colv=TRUE, Rowv=FALSE, scale="none", dendrogram="none", col=ramp(200), margins=c(16,14), trace="none", density.info="none", key=TRUE, colsep=seq(1:length(unique(names(sp)))), rowsep=seq(1:length(unique(names(sp)))), breaks=col_breaks, symkey=FALSE, key.xlab="log(EC50)", keysize=1)
  dev.off()
  
  # Add cytotox data but order heatmap according to default dendrogram
  q <- q[,c(1:16,19)]
  names(q)[17] <- "Viability"
  
  pdf("ec50_heatmap_withCytotox.pdf", width=11)
  heatmap.2(log(t(q)), Colv=TRUE, Rowv=FALSE, scale="none", dendrogram="col", col=ramp(200), margins=c(16,14), trace="none", density.info="none", key=TRUE, colsep=seq(1:length(unique(names(sp)))), rowsep=seq(1:length(unique(names(sp)))), breaks=col_breaks, symkey=FALSE, key.xlab="log(EC50)", keysize=1)
  dev.off()
}






if (set_output == "selectivity_heatmap.pdf") {
  ## Calculate log-ratio of network parameter EC50 value to cytotoxicity EC50
  ## This might require first setting undefined viability EC50 values to maxmimum tested concentrations.
  q[,"min"] <- apply(q[,17:18], 1, min)
  q <- q[order(q[,"min"]),]
  
  p <- log(q/q[,19])
  p <- p[,1:16] # remove cytotox data
  p[p > 0] <- 0 # set values greater than zero to zero (cases of network parameter EC50 greater than cytotox EC50)
  
  
  # set new color scheme for different range
  ramp<-colorRampPalette(c("cyan","cornflowerblue","gray20"),space="rgb")
  col_breaks = c(seq(-7,-2,length=100), # for cyan
                 seq(-1.99,-.01,length=100), # for cornflowerblue
                 seq(0,1,length=1)) # for black
  
  pdf("ec50_heatmap_selectivity_v1.pdf", width=11)
  heatmap.2(t(p), Colv=TRUE, Rowv=FALSE, scale="none", dendrogram="col", col=ramp(200), margins=c(16,14), trace="none", density.info="none", key=TRUE, colsep=seq(1:length(unique(names(sp)))), rowsep=seq(1:length(unique(names(sp)))), breaks=col_breaks, symkey=FALSE, key.xlab="log(Network EC50 / Cytotox EC50)", keysize=1)
  dev.off()
  
  
  ## Establish column clustering by manually calling hierarchical clustering (Not working currently)
  colDistance = dist(log(q[,1:16]), method = "euclidean")
  colCluster = hclust(colDistance, method = "complete")
  colDend = as.dendrogram(colCluster)
  #colDend = reorder(colDend, t(q[,1:16]))
  
  heatmap.2(t(p), Colv=colDend, Rowv=FALSE, scale="none", dendrogram="none", col=ramp(200), margins=c(16,14), trace="none", density.info="none", key=TRUE, colsep=seq(1:length(unique(names(sp)))), rowsep=seq(1:length(unique(names(sp)))), breaks=col_breaks, symkey=FALSE, key.xlab="log(Network EC50 / Cytotox EC50)", keysize=1)
    
}







## Plotting relationship between network parameter effects vs. cytotoxicity
sp_ntwk_vs_tox <- function(q, set_parameter) {
  
  # q is table of EC50 values with each row a compound and each column a network parameter. 
  # the last two columns are assumed to contain cytotoxicity EC50 values
  
  q[is.na(q)] <- 30 # Setting all NA values to maximum tested dose of 30 uM
  q[,"min"] <- apply(q[,(ncol(q)-1):ncol(q)], 1, min) # Find minimum of the two cytotox assay values
  q_t <- log(30/q) # This is the same as -(log(q/30)) ; log of fold-concentration difference between EC50 value and maximum tested dose
  
  # Plotting transformed EC50 of cytotoxicity vs. transformed EC50 of network parameters
  plot(q_t[,"min"], q_t[,set_parameter], pch=16, ylab="Network Parameter: -log(EC50/max dose)", main=set_parameter, xlab="Viability: -log(EC50/max dose)", ylim=c(-1,8), xlim=c(-1,8), bty="n")
  abline(0,1, lwd=2, col="gray20", lty=2)
  text(q_t[,"min"], q_t[,set_parameter], labels=row.names(q), cex= 0.7, pos=3, offset=0.3)
  
  
  ## Calculate selectivity scores for each compound
  euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2)) # Euclidean distance function
  midpoint <- function(x1, x2) c((x1[1] + x2[1])/2, (x1[2] + x2[2])/2) # Find midpoint function
  selscores <- list()
  
  # Loop through each compound
  for (i in unique(row.names(q))){
    
    poi <- c(q_t[i,"min"], q_t[i, set_parameter]) # get coordinates of point of interest
    mid <- midpoint(c(q_t[i,"min"], q_t[i,"min"]), c(q_t[i, set_parameter], q_t[i, set_parameter])) # find midpoint on y=x line
    distA <- euc_dist(poi, mid) # calculate distance from point of interest to midpoint
    
    # Reverse the sign if cytotox value is greater than network parameter value
    if (q_t[i,"min"] > q_t[i,set_parameter]) {
      distA <- -distA
    }
       
    #distB <- (q_t[i,set_parameter] - q_t[i,"min"]) # Find raw difference between network parameter and viability values
    
    selscores[[i]] <- distA
    #selscores[[i]]$distB <- distB
  }
   
  sltv_scores <- as.data.frame(do.call(rbind, selscores)) # make data frame
  class(sltv_scores[,1]) <- "numeric"
  
  sltv_scores
}




if (set_output == "selectivity") {
  
  selectivity <- list()
  pdf("Network_vs_Cytotox_Scatter.pdf")
  for (i in colnames(q)) {
    selectivity[[i]] <- sp_ntwk_vs_tox(q, i)
  }
  dev.off()
  
  # Produce table of selectivity results
  sltv_df <- data.frame(lapply(selectivity, cbind)) # Form data frame
  names(sltv_df) <- c(names(q)) # Rename columns
  sltv_df <- sltv_df[,1:(ncol(sltv_df)-2)] # Remove cytotox columns
  sltv_df[,"mean"] <- apply(sltv_df, 1, mean) # Find mean for each compound
  sltv_df[,"median"] <- apply(sltv_df, 1, median) # Find median for each compound
  sltv_df <- sltv_df[order(sltv_df[,"mean"], decreasing=TRUE),] # Sort results by mean value
  write.table(sltv_df, file="selectivity_scores.csv", sep=",") # Write out table of results
  
  # Produce table of potency results
  np_table <- q[,1:(ncol(q)-2)] # Subset EC50 values
  np_table[is.na(np_table)] <- 30 # Set undetermined EC50 values (NA) to max dose tested (30)
  np_table <- log(30/np_table) # transform values to new scale
  np_table[,"mean"] <- apply(np_table, 1, mean)
  np_table[,"median"] <- apply(np_table, 1, median)
  np_table <- np_table[order(np_table[,"mean"], decreasing=TRUE),]
  write.table(np_table, file="networkparam_scores.csv", sep=",")
  
  # Generating summary rankings plot
  pdf("Potency_vs_Selectivity_RankingsSummary.pdf")
  np_table[,"rank"] <- c(1:nrow(np_table))
  sltv_df[,"rank"] <- c(1:nrow(sltv_df))
  comb_rank <- merge(np_table, sltv_df, by="row.names")
  plot(comb_rank$rank.x, comb_rank$rank.y, pch=16, xlim=c(-2,80), ylim=c(-2,80), xlab="Potency rank", ylab="Selectivity rank", main="Ranking summary", cex.lab=1.2, cex.main=1.2, cex.axis=1.2, bty="n")
  text(comb_rank$rank.x, comb_rank$rank.y, labels=comb_rank[,"Row.names"], cex=0.6, pos=3, offset=0.3)
  dev.off()
  
}






## How to recover rearranged order of heatmap.2 hierarchical clustering
#q[is.na(q)] <- 100
#h <- heatmap.2(log(t(q)), Colv=TRUE, Rowv=FALSE, scale="none", dendrogram="col")
#q2 <- t(t(q)[rev(h$rowInd), h$colInd])
#q2 <- data.frame(q2, seq=seq_len(nrow(q2))) # add column of row order

# Merge with another data frame to reorder that data
#q_dnt <- merge(q2, dnt, by="row.names")

# Re-establish dendrogram order
#q_dnt <- q_dnt[order(q_dnt$seq),]

# Map values to colors
#dnt_cols <- recode(q_dnt[,"evidence"], "0='green'; 1='gray20'; 2='orange'; 3='red'; 4='yellow'")







