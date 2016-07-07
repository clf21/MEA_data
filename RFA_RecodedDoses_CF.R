# MEA neuronal ontogeny random forest analysis (RFA) on plate-normalized values
# Recoding doses and classifying all data together based on dose code
# Chris Frank - March 2016


# Set parameters of analysis:
set_directory <- "~/MEA_Data/PIP3_AIM1/TEST4/" # Set directory 
set_data_file <- "Final_Data_Set_SA1_DNT_Paper1 (2)(updated)_CF.csv" # Name of csv file with burst parameters per sample
set_output_file <- "RFA_results14.txt"


# Load required libraries, initialize random number generation
library(randomForest)
library(car)
set.seed(3)


# Read in data and clean it up
setwd(set_directory)
all_data <- read.delim(set_data_file, sep=",", stringsAsFactors=FALSE) #Reading in full data set ; assuming it is a csv table
all_data[is.na(all_data)] <- 0 # Replace all NAs with zeros - This may be undesirable for MEA parameters that are derived from other parameters.
all_data <- subset(all_data, trt!="Glyphosate")
all_data[all_data[,"dose"]==0,"trt"] <- "Control" # convert all zero dose rows to say "Control" for treatment
all_data <- all_data[, !(names(all_data) %in% c("cv.time","cv.network"))] #Remove "cv.time" and "cv.network" parameters from analysis





## Normalize data by % of control wells at DIV12 on same plate for a particular ontogeny parameter
normalize_by_plate <- function(all_data) {
  
  # Split data set by plate (date as well because plateIDs get reused)
  splt.by <- c('date','Plate.SN')
  per_plate_split <- split(all_data, all_data[,splt.by])
  per_plate_split2 <- per_plate_split[sapply(per_plate_split, function(x) dim(x)[1]) > 0] # remove empty plate + date combinations
  rm(per_plate_split) 
  
  # Normalize each plate by percent of control mean at DIV12 for that parameter and combine into single data frame
  norm_plates <- list()
  
  # Loop through each plate
  for (i in per_plate_split2) {
    
    # Loop through each ontogeny parameter
    for (j in names(all_data[,8:(ncol(all_data)-1)])) {
      i[,j] <- ((i[,j]) / mean(subset(i, trt=="Control" & DIV==12)[,j]))*100
    }
    
    norm_plates[[length(norm_plates)+1]] <- i
  }
  
  rm(i,j)
  norm_plates <- do.call(rbind, norm_plates) #Re-form one table of values
  
  norm_plates  
}




## Function to calculate accuracy and important features from trees
#Input is dataframe with first column of treated vs. untreated classes (data.df), output is prediction accuracy and average importance of each feature, based on its Mean Decrease in Gini
boost.tree.all <- function(data.df, num.iter=500) {
  ret<-NULL
  # data.df[,1] =  treatment factor
  classes<-data.df[,1]
  # make new data frame
  class.df<-cbind(classes, data.df[,-1])
  n.row<-nrow(class.df)
  #Set size of training set to 2/3 of total data
  Ntrain<-round(n.row*2/3)
  correct<-0 # initiate % correct to 0
  imp<-rep(0, (ncol(class.df)-1)) # initiate importance of each input variable
  #Perform 500 trials with different training sets
  
  for (i in 1:num.iter) {
    #print(i) #to see what loop we're on
    
    # randomly sample Ntrain observations from data
    traindat<-sample(1:nrow(class.df), Ntrain)
    
    # if there are not as many classes present in sampled data as full input set, discard and repeat sampling
    while(length(unique(class.df[traindat,]$classes)) < length(unique(classes))) {
      traindat <- sample(1:nrow(class.df), Ntrain)
    }
    
    # the rest in test data set
    testdat<-class.df[-traindat ,]
    # Calculate boosted classification tree, using 500 trees
    # classes are DIV and trt, '.' indicates all other variabels, 
    tree.out<-randomForest(classes~. , data=class.df[traindat,], importance=TRUE, ntree=500 )
    # using this decision tree see how well the classes were predicted.
    # print(tree.out) # Uncomment to look at confusion matrix
    # predict held-out data classes
    tree.pred<-predict(tree.out,testdat,type="class")
    # print(tree.pred)
    # compare to actual predictions ; these functions would use ROCR package to plot ROC curves
    # tree.pred.prob <- predict(tree.out, testdat, type="prob")
    # pred_eval <- prediction(tree.pred.prob, testdat$classes)
    # roc <- performance(pred_eval, "tpr", "fpr")
    # plot(roc)
    # correct is the %correctly predicted + correct before
    correct<- sum(tree.pred==testdat$classes)/(n.row-Ntrain) + correct
    # for each variable the sum of the decrease in gini is computed across all nodes that used that variable in decision
    # the larger the number the more important that variable is in classification.
    imp<-importance(tree.out)[,"MeanDecreaseGini"]+imp
  }
  
  avg.imp<-imp/num.iter # get average importance over 500 times
  ret$accuracy<-correct/num.iter # get average % correct
  ret$factors<-avg.imp
  ret
}



## Function to prepare data for RFA, run it, and display results
runRFA <- function(all_norm_data) {
  
  # loop through DIV times in recoded data and run RFA
  rfa_out <- list()
  j <- 1
  for (i in c(5,7,9,12)) {
    data_sub <- norm_data_recoded[norm_data_recoded[,"DIV"]==i,] # select data from one DIV
    rfa_input <- data.frame(dose=as.factor(data_sub$dose), data_sub[,8:(ncol(data_sub)-1)]) # format RFA input df  
    rfa_out[[j]] <- boost.tree.all(rfa_input, num.iter=500) # perform RFA, iterate 500 times
    print(paste("DIV", i, sep=""))
    print(rfa_out[[j]])
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









### Generate outputs on Aim1 data

# normalize by plate DIV12 controls
all_norm_data <- normalize_by_plate(all_data)

# Code the different doses the same between compounds
# split normalized data by compound
norm_split <- split(all_norm_data, all_norm_data[,"trt"])

# recode doses by counting down from max dose
j <- 1
for (i in norm_split) {
  ln <- length(unique(i$dose)) # get number of doses tested for that compound
  
  # Controls will stay marked as zero dose
  
  if (ln == 6) {
    norm_split[[j]]$dose <- recode(i$dose, "0 = 0 ; unique(i$dose)[ln] = 6 ; unique(i$dose)[ln-1] = 5 ; unique(i$dose)[ln-2] = 4 ; unique(i$dose)[ln-3] = 3 ; unique(i$dose)[ln-4] = 2 ; unique(i$dose)[ln-5] = 1")
  }
  
  if (ln == 7) {
    norm_split[[j]]$dose <- recode(i$dose, "0 = 0 ; unique(i$dose)[ln] = 6 ; unique(i$dose)[ln-1] = 5 ; unique(i$dose)[ln-2] = 4 ; unique(i$dose)[ln-3] = 3 ; unique(i$dose)[ln-4] = 2 ; unique(i$dose)[ln-5] = 1 ; unique(i$dose)[ln-6] = 1")
  }
  
  if (ln == 8) {
    norm_split[[j]]$dose <- recode(i$dose, "0 = 0 ; unique(i$dose)[ln] = 6 ; unique(i$dose)[ln-1] = 5 ; unique(i$dose)[ln-2] = 4 ; unique(i$dose)[ln-3] = 3 ; unique(i$dose)[ln-4] = 2 ; unique(i$dose)[ln-5] = 1 ; unique(i$dose)[ln-6] = 1 ; unique(i$dose)[ln-7] = 1")
  }
  
  j <- j+1
}

rm(j)
norm_data_recoded <- do.call(rbind, norm_split)


# Run RFA, write results to file
sink(file=set_output_file, split=TRUE) # specify output file
runRFA(norm_data_recoded)
sink() # close output file

