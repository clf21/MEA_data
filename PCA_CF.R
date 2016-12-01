## Principal components analysis (PCA) to evaluate parameter importance and redundancy
## Chris Frank - April 2016


set_directory <- "~/MEA_Data/PIP3_AIM2-3-DNT/" # Set directory 
set_data_file1 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/SA2_allData_11-17-2014.csv" # Name of csv file with burst parameters per sample
set_data_file2 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/AllData_SA3_NoDNTRef.csv"
set_data_file3 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/DNTRef_all.csv"
set_mi_data_file <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/MutualInformation/update_10-31-16/all_mi.csv"


# Load packages
library(ggplot2)
library(ggfortify)



## Read in full burst parameter dataset
setwd(set_directory)
all_data1 <- read.delim(set_data_file1, sep=",", stringsAsFactors=FALSE) #Reading in full data set ; assuming it is a csv table
all_data2 <- read.delim(set_data_file2, sep=",", stringsAsFactors=FALSE) # Reading in second data set
all_data3 <- read.delim(set_data_file3, sep=",", stringsAsFactors=FALSE) # Reading in third data set
all_data <- rbind(all_data1, all_data2, all_data3)
rm(all_data1, all_data2, all_data3)

bis_rows <- grep("12_01_", all_data$file.name, fixed=TRUE) #index all Bicuculline-treated wells
all_data <- all_data[- bis_rows,] #Remove all bic-treated wells
all_data <- subset(all_data, DIV %in% c(5,7,9,12)) # Remove sparse DIV2 data
all_data <- unique(all_data) # eliminate duplicate data


## Read in mutual information parameter and attach to full dataset
mi_data <- read.delim(set_mi_data_file, sep=",", stringsAsFactors=FALSE)
mi_data <- subset(mi_data, DIV %in% c(5,7,9,12)) # Remove sparse DIV2 data
mi_data <- unique(mi_data) # eliminate duplicate data
all_data <- merge(all_data, mi_data, all.x=TRUE) # merge two data frames on common columns (plateID, date, well, treatment, etc.)
all_data <- all_data[, c(1:23, 27, 24:26)] # rearrange data frame to accomodate additional parameter
names(all_data)[[24]] <- "mi"



## Missing data values are an issue
## Option 1 - impute missing data as mean of controls on same DIV (see below)

## Option 2 - impute missing data as zero
#all_data[is.na(all_data)] <- 0 # Replace all NAs with zeros - This may be undesirable for MEA parameters that are derived from other parameters.

## Option 3 - remove all rows with missing data (removes ~29% of data)
#all_data <- na.omit(all_data)



## Option 4 - first replace all NA values with mean of defined values for same compound + dose + DIV (This represents ~1/3rd of NAs)
## Then replace remaining NA values with zeros (~2/3rds NAs)
impute_na <- function(all_data) {
  # Split data by DIV, treatment, and dose
  trt_split <- split(all_data, interaction(all_data[,"DIV"], all_data[,"trt"], all_data[,"dose"], drop=TRUE))

  cmplt_data <- list() #initialize list

  # Loop through each treatment
  for (i in trt_split) {
    
    # Loop through each network parameter
    for (j in names(all_data[,8:(ncol(all_data)-1)])) {
      
      # Loop through each row
      for (k in 1:nrow(i)) {
        
        # If value is NA, replace with mean of determined values
        if (is.na(i[k,j]) && length(na.omit(i[,j])) > 0) {
          i[k,j] <- median(na.omit(i[,j]))
        }
        
        # If no value to replace with, set to zero
        if (is.na(i[k,j]) && length(na.omit(i[,j]))==0) {
          i[k,j] <- 0
        } 
      }
    }
    
    # Add corrected data to list
    cmplt_data[[length(cmplt_data)+1]] <- i
  }

  cmplt_data <- do.call(rbind, cmplt_data) #Re-form one table of values
  cmplt_data
}












## Examining PCA of data to determine important network parameters
pdf("PCA_results.pdf")
sink("PCA_results.txt")
for (i in c(5,7,9,12)) {
  
  div_data <- subset(all_data, DIV==i)
  
  names(div_data) <- c("date","Plate.SN","DIV","well","trt","dose","units","MFR","BR","ISI","%SiB","BD","IBI","#AE","#ABE","#NS","NSP","NSD","%SiNS","NS-ISI","NSDsd","#SiNS","R","MI","CVT","CVN","file.name")
  
  #autoplot(prcomp(div_data[,8:24], scale=TRUE), data=div_data, colour="dose", loadings=TRUE, loadings.label=TRUE)
  
  # Perform the PCA
  pca <- prcomp(div_data[,8:24], scale=TRUE)
  
  # Get factor loadings
  eig <- (pca$sdev)^2
  variance <- eig*100/sum(eig)
  cumvar <- cumsum(variance)
  
  # Make scree plot (cumulative variance explained by PCs)
  plot(cumvar, pch=16, ylim=c(0,100), ylab="% Variance explained", xlab="# PCs", main=paste("DIV ", i, sep=""), cex.lab=1.4, cex.axis=1.2)
  lines(cumvar, lwd=2)
  
  # Make unit circle plot with variable loadings (rotation matrix)
  #theta <- seq(0,2*pi,length.out = 100)
  #circle <- data.frame(x = cos(theta), y = sin(theta))
  #p <- ggplot(circle, aes(x,y)) + geom_path()
  
  #loadings <- data.frame(pca$rotation, .names = row.names(pca$rotation))
  #p <- p + geom_text(data=loadings, mapping=aes(x = PC1, y = PC2, label = .names, colour = .names)) + coord_fixed(ratio=1) + labs(x = "PC1", y = "PC2")
  #print(p)  
  
  # Print variable loadings on first 8 principal components
  print(paste("DIV",i,sep=""))
  print(summary(pca))
  print(abs(pca$rotation[,1:8]))
  
  # Convert variable loadings to percent of total
  aload <- abs(pca$rotation[,1:8])
  prop <- (sweep(aload, 2, colSums(aload), "/"))
  print(prop)
  
  barplot(prop[,1:6], beside=T, col=palette(rainbow(n=17)), ylab="Proportion of PC", main=paste("DIV ", i, sep=""), ylim=c(0,max(prop[,1:6])+.1))
  legend("topleft", legend=names(div_data[,8:24]), fill=palette(rainbow(n=17)), ncol=4)
  
  # Make biplot with variable loadings on PC1 vs. PC2
  biplot(pca, col=c("blue","black"), xlabs=rep("o", nrow(div_data)), cex=c(0.6,1.2), cex.axis=1.2, cex.lab=1.4, main=paste("DIV ", i, sep=""))
  
}
dev.off()
sink()







## Look at overall correlation between parameters:
library(colorRamps)
library(gplots)
library(corrplot)

# Change names of dataframe
names(all_data) <- c("date","Plate.SN","DIV","well","trt","dose","units","MFR","BR","ISI","%SiB","BD","IBI","#AE","#ABE","#NS","NSP","NSD","%SiNS","NS-ISI","NSDsd","#SiNS","R","MI","CVT","CVN","file.name")

# Produce heatmap of correlation matrix
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white", "cyan", "#007FFF", "blue", "#00007F"))

pdf("ParameterCorrelationMap.pdf")
corrplot(cor(scale(all_data[,8:24], center=TRUE, scale=TRUE), method="pearson"), method="ellipse", type="upper", order="hclust", col=rev(col1(200)), outline=TRUE, tl.col="black")
dev.off()

## We may want to consider how inclusion of imputed NA values changes the parameter relationships
## Note - read in fresh data set before using this.
pdf("ParameterCorrelationMap_NAomitted.pdf")
corrplot(cor(scale(na.omit(all_data[,8:24]), center=TRUE, scale=TRUE), method="pearson"), method="ellipse", type="upper", order="hclust", col=rev(col1(200)), outline=TRUE, tl.col="black")
dev.off()



## Plotting correlation matrix at different timepoints
# Note - the order parameter has been changed to "original" to keep all the same for comparison
pdf("ParameterCorrelationMap_NAomit_DIV12.pdf")
corrplot(cor(scale(na.omit(subset(all_data, DIV==12)[,8:24]), center=TRUE, scale=TRUE), method="pearson"), method="ellipse", type="upper", order="original", col=rev(col1(200)), outline=TRUE, tl.col="black")
dev.off()

pdf("ParameterCorrelationMap_NAomit_DIV9.pdf")
corrplot(cor(scale(na.omit(subset(all_data, DIV==9)[,8:24]), center=TRUE, scale=TRUE), method="pearson"), method="ellipse", type="upper", order="original", col=rev(col1(200)), outline=TRUE, tl.col="black")
dev.off()

pdf("ParameterCorrelationMap_NAomit_DIV7.pdf")
corrplot(cor(scale(na.omit(subset(all_data, DIV==7)[,8:24]), center=TRUE, scale=TRUE), method="pearson"), method="ellipse", type="upper", order="original", col=rev(col1(200)), outline=TRUE, tl.col="black")
dev.off()

pdf("ParameterCorrelationMap_NAomit_DIV5.pdf")
corrplot(cor(scale(na.omit(subset(all_data, DIV==5)[,8:24]), center=TRUE, scale=TRUE), method="pearson"), method="ellipse", type="upper", order="original", col=rev(col1(200)), outline=TRUE, tl.col="black")
dev.off()

pdf("ParameterCorrelationMap_NAomit_AllDIV.pdf")
corrplot(cor(scale(na.omit(all_data[,8:24]), center=TRUE, scale=TRUE), method="pearson"), method="ellipse", type="upper", order="original", col=rev(col1(200)), outline=TRUE, tl.col="black")
dev.off()



#heatmap.2(cor(scale(all_data[,8:24], center=TRUE, scale=TRUE), method="pearson"), col=ramp(200), scale="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, margins=c(8,7), key=TRUE, density.info="none", colsep=1:16, rowsep=1:16, key.xlab="Correlation", keysize=1, trace="none", symkey=FALSE)


## Getting out ordered list of correlation coefficients
cor_mat <- cor(scale(na.omit(all_data[,8:24]), center=TRUE, scale=TRUE), method="pearson")
cor_list <- as.data.frame(as.table(cor_mat))
cor_ord <- cor_list[order(abs(cor_list[,"Freq"]), decreasing=TRUE),]
write.table(cor_ord, "PearsonCorr_Coefficients_List.txt", quote=FALSE, sep="\t", row.name=F)




## Imputing control values for 'NA' marked parameters
#impute_na <- function(all_data) {
#  
#  # Split data set by DIV
#  div_split <- split(all_data, all_data[,"DIV"]) 
#  
#  cmplt_data <- list() #initialize list
#  
#  # Loop through each DIV
#  for (i in div_split) {
#    
#    # Loop through each ontogeny parameter
#    for (j in names(all_data[,8:(ncol(all_data)-1)])) {
#      
#      cntrls <- na.omit(subset(i, trt=="Control")[,j]) # Get vector of control values for that parameter
#      i[,j][is.na(i[,j])] <- median(cntrls) # Replace NA values with median of control values      
#    } 
#    
#    cmplt_data[[length(cmplt_data)+1]] <- i
#  }
#  
#  rm(i,j)
#  cmplt_data <- do.call(rbind, cmplt_data) #Re-form one table of values
#  
#  cmplt_data  
#}