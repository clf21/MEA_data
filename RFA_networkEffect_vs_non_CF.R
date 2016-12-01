
set.seed(333)

set_directory <- "~/MEA_Data/PIP3_AIM2-3-DNT/" # Set directory 
set_data_file1 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/SA2_allData_11-17-2014.csv" # Name of csv file with burst parameters per sample
set_data_file2 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/AllData_SA3_NoDNTRef.csv"
set_data_file3 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/DNTRef_all.csv"
set_mi_data_file <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/MutualInformation/update_10-31-16/all_mi.csv"
set_ec50_file <- "~/MEA_Data/PIP3_AIM2-3-DNT/TEST5/combined_ec50_allOntogeny.txt"


# Load required libraries, initialize random number generation
library(randomForest)
library(car)
library(ROCR)


## First, a bunch of data formatting is needed ##

## Read in and clean up EC50 table
ec50_table <- read.delim(set_ec50_file, header=TRUE, stringsAsFactors=FALSE)
ec50_table <- subset(ec50_table, !ontogeny %in% c("cv.time_AUC", "cv.network_AUC")) # Remove CV.time and CV.network unused parameters
ec50_table[is.na(ec50_table)] <- 10000 # convert all NA values to arbitrary high value


## Determine minimum EC50 per compound
min_EC50 <- aggregate(ec50_table[,"ec"], by=list(ec50_table[,"compound"]), FUN=min)
names(min_EC50) <- c("trt","cutoff_dose")


## Determine minimum cytotoxicity assay EC50 per compound
cytotox_EC50 <- subset(ec50_table, ontogeny %in% c("ABcytotox", "LDHcytotox"))
min_ctx_EC50 <- aggregate(cytotox_EC50[,"ec"], by=list(cytotox_EC50[,"compound"]), FUN=min)
names(min_ctx_EC50) <- c("trt","cytotox_dose")


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



## Label dataset with above or below minimum EC50
all_data_plus <- merge(all_data, min_EC50, by="trt")
all_data_plus[,"effect"] <-(all_data_plus[,"dose"] > all_data_plus[,"cutoff_dose"]) # marks with True when tested dose exceeds min. NP EC50


## Use cytotox EC50 min. to remove cytotoxic treatments to focus on specific network effects
all_data_plus2 <- merge(all_data_plus, min_ctx_EC50, by="trt")
all_data_plus2[,"cytotoxic"] <- (all_data_plus2[,"dose"] > all_data_plus2[,"cytotox_dose"]) # marks with True when tested dose exceeds cytotox EC50
all_data_plus <- subset(all_data_plus2, cytotoxic==FALSE) # replace all_data_plus with reduced dataset
all_data_plus <- all_data_plus[,1:(ncol(all_data_plus)-2)] # remove cytotox info
rm(all_data_plus2, all_data)





## Function to normalize data by % of control well median at DIV12 on same plate 
normalize_by_plate <- function(all_data) {
  
  # Split data set by plate (date as well because plateIDs get reused occasionally)
  per_plate_split <- split(all_data, interaction(all_data[,"date"], all_data[,"Plate.SN"], drop=TRUE))
  
  # Normalize each plate by percent of untreated control median at DIV12 for that parameter and combine into single data frame
  norm_plates <- list()
  
  # Loop through each plate
  for (i in per_plate_split) {
    
    # Loop through each ontogeny parameter
    for (j in names(all_data[,8:(ncol(all_data)-3)])) {
      i[,j] <- ((i[,j]) / median(subset(i, dose==0 & DIV==12)[,j]))*100
      
      #if(median(subset(i, dose==0 & DIV==12)[,j]) == 0){ 
      #  print(paste("Median for plate", i[1,"Plate.SN"], j, "is zero! - Normalization will produce infinite values"))
      #}
    }
    
    norm_plates[[length(norm_plates)+1]] <- i
  }
  
  rm(i,j)
  norm_plates <- do.call(rbind, norm_plates) #Re-form one table of values
  
  norm_plates  
}









## Function to calculate accuracy and important features from trees
# Input is dataframe with first column of treated vs. untreated classes, 
# Output is prediction accuracy and average importance of each feature, based on its Mean Decrease in Gini
boost.tree.all <- function(data.df, num.iter=500) {
  ret <- NULL
  
  # set first column as classification key
  classes <- as.factor(data.df[,1])
  
  #Set size of training set to 2/3 of total data
  n.row <- nrow(data.df)
  Ntrain <- round(n.row*2/3)
  
  # initiate % correct to 0 and importance of each input variable
  correct <- 0 
  imp <- rep(0, (ncol(data.df)-1))
  
  
  #Perform 500 trials with different training sets
  roc <- list()
  n <- 1
  for (i in 1:num.iter) {
    #print(i) #to see what loop we're on
    
    # randomly sample Ntrain observations from data
    traindat <- sample(1:nrow(data.df), Ntrain)
    
    # use the rest as a test data set
    testdat <- data.df[-traindat,]
    
    # random forest call using all parameters to predict class
    # na.omit is being used to handle cases of undetermined values
    tree.out <- randomForest(as.factor(classes) ~ ., data=data.df[traindat,], importance=TRUE, ntree=500)
    
    # using this decision tree see how well the classes were predicted.
    print(tree.out) # Uncomment to look at confusion matrix
    # predict held-out data classes
    tree.pred <- predict(tree.out, testdat, type="class")
    # print(tree.pred)
    
    # compare to actual predictions ; these functions use ROCR package to plot ROC curves
    tree.pred.prob <- predict(tree.out, testdat, type="prob") # output probabilities of binary classifier
    pred_eval <- prediction(tree.pred.prob[,2], testdat$classes) # match prediction with true labels
    roc[[n]] <- performance(pred_eval, "tpr", "fpr") # evaluate performance with ROC curve ; add to list to look at all curves
    n <- n + 1 # increment n
    #plot(roc)
    
    # correct is the percent correctly predicted + correct before
    correct <- sum(tree.pred==testdat$classes)/(n.row-Ntrain) + correct
    
    # for each variable the sum of the decrease in gini is computed across all nodes that used that variable in decision
    # the larger the number the more important that variable is in classification.
    imp <- importance(tree.out)[,"MeanDecreaseGini"] + imp
  }
  
  plot(roc[[1]]) # plot first ROC curve
  for (curve in roc){
    plot(curve, add=TRUE) # plot all iterations ROC curves
  }
  
  avg.imp <- imp/num.iter # get average importance over number of iterations performed
  ret$accuracy <- correct/num.iter # get average % classes correct
  ret$factors <- avg.imp
  ret

}















## Function to prepare data for RFA, run it, and display results
runRFA <- function(norm_data) {
  
  # loop through each DIV and run RFA
  rfa_out <- list()
  j <- 1
  for (i in c(5,7,9,12)) {
    
    data_sub <- norm_data[norm_data[,"DIV"]==i,] # select data from one DIV
    
    rfa_input <- data_sub[,c(ncol(norm_data), 8:(ncol(norm_data)-5))] # reformat for RFA input
    colnames(rfa_input)[1] <- "classes"
    rfa_input[is.na(rfa_input)] <- 0 # set all undefined values to 0. RFA does not know how to deal with undefined predictor variables
    
    rfa_out[[j]] <- boost.tree.all(rfa_input, num.iter=10) # perform RFA, iterate 100 times
    
    print(paste("DIV", i, sep=""))
    print(rfa_out[[j]]) # report results for that DIV
    j <- j+1
  }
  
  # combine DIV outputs
  DIV_results <- data.frame(cbind(rfa_out[[1]]$factors, rfa_out[[2]]$factors, rfa_out[[3]]$factors, rfa_out[[4]]$factors), stringsAsFactors=FALSE)  
  names(DIV_results) <- c("DIV5","DIV7","DIV9","DIV12")
  
  # calculate mean accuracy across DIVs
  mean_accuracy <- mean(c(rfa_out[[1]]$accuracy, rfa_out[[2]]$accuracy, rfa_out[[3]]$accuracy, rfa_out[[4]]$accuracy))
  
  # take mean of DIV results per parameter and sort highest to lowest
  factor_rankings <- sort(rowMeans(DIV_results), decreasing=T)
  
  # report summary
  print("Summary across DIVs:")
  print(paste("Mean accuracy = ", round(mean_accuracy, 3)))
  print(factor_rankings) 
  print("____________________________________________________________________________________")
}













# Apply normalization
norm_data <- normalize_by_plate(all_data_plus)

# Perform RFA on each DIV
runRFA(norm_data)









library("rpart")
library("rattle")
data_sub <- norm_data[norm_data[,"DIV"]==12,]
rfa_input <- data_sub[,c(ncol(norm_data), 8:(ncol(norm_data)-5))]
fit2 <- rpart(as.factor(effect) ~ ., data=rfa_input, method="class", control=rpart.control(minsplit=10, cp=0.01))
fancyRpartPlot(fit2)


