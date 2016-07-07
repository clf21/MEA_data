## Performing intitial steps of tipping point analysis on normalized data
## Chris Frank - February-May 2016


set_directory <- "~/MEA_Data/PIP3_AIM2-3-DNT/" # Set directory 
set_data_file1 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/SA2_allData_11-17-2014.csv" # Name of csv file with burst parameters per sample
set_data_file2 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/AllData_SA3_NoDNTRef.csv"
set_data_file3 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/DNTRef_all.csv"
set_data_file4 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/AB_CytotoxData.csv"
set_data_file5 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/LDH_CytotoxData.csv"
set_output <- "replicates"  # Can be set to "zscore_plots" (plot individual parameters),  "total_perturb" (plot sum of individual parameter z-scores), "total_perturb_slope" (plot sum of z-scores and slope changes over dose), or "replicates" to look at scalar perturb between plates run
set_normalization <- "plate" # Can be set to "global" or "plate"
set_parameters <- c("MFR","#AE","BR","#ABE","%SiB","#NS","%SiNS","R") # A vector of parameters used for the total perturbation calculations ; can include any of names(all_data[8:(ncol(all_data)-3)])              


set.seed(1)
# Load packages
library(ggplot2)


# Source and clean up dataset
setwd(set_directory)
all_data1 <- read.delim(set_data_file1, sep=",", stringsAsFactors=FALSE) # Reading in first data set ; assuming it is a csv table
all_data2 <- read.delim(set_data_file2, sep=",", stringsAsFactors=FALSE) # Reading in second data set
all_data3 <- read.delim(set_data_file3, sep=",", stringsAsFactors=FALSE) # Reading in third data set
all_data <- rbind(all_data1, all_data2, all_data3)
rm(all_data1, all_data2, all_data3)

bis_rows <- grep("12_01_", all_data$file.name, fixed=TRUE) # Index all Bicuculline-treated wells
all_data <- all_data[- bis_rows,] # Remove all bic-treated wells
all_data <- all_data[all_data[,"DIV"]!=2,] # Remove sparse DIV2 data
all_data[all_data[,"dose"]==0,"trt"] <- "Control" # Rename all zero dose wells to say Control
names(all_data) <- c("date","Plate.SN","DIV","well","trt","dose","units","MFR","BR","ISI","%SiB","BD","IBI","#AE","#ABE","#NS","NSP","NSD","%SiNS","NS-ISI","NSDsd","#SiNS","R","CVT","CVN","file.name") # Rename network parameters

#all_data[is.na(all_data)] <- 0 # Replace all NAs with zeros - This may be undesirable for MEA parameters that are derived from other parameters.
#all_data[,8:(ncol(all_data)-1)] <- log(all_data[,8:(ncol(all_data)-1)] + 1) # Uncomment to log-transform raw data values ahead of all normalization


# Source cell viability data
viability_data1 <- na.omit(read.delim(set_data_file4, sep=",", header=TRUE, stringsAsFactors=FALSE))
viability_data2 <- na.omit(read.delim(set_data_file5, sep=",", header=TRUE, stringsAsFactors=FALSE))
all_via <- merge(viability_data1, viability_data2)
all_via[,"median"] <- apply(all_via[,3:8], 1, median)
all_via[,"median"] <- all_via[,"median"]*100
names(all_via)[1] <- "trt" # changing column heading to match network data






## Normalize using global z-score (subtract median of control wells on the same DIV across ALL plates, then divide by SD of controls)
normalize_by_zscore <- function(all_data) {
  
  # Split data set by DIV
  div_split <- split(all_data, all_data[,"DIV"]) 
  
  norm_plates <- list() #initialize list
  
  # Loop through each DIV
  for (i in div_split) {
    
    # Loop through each ontogeny parameter
    for (j in names(all_data[,8:(ncol(all_data)-1)])) {
      
      cntrls <- (subset(i, trt=="Control")[,j]) # Get vector of control values for that DIV
      i[,j] <- (((i[,j]) - median(cntrls))/ sd(cntrls)) # Subtract all values on that DIV by median, divide by SD
    } 
    
    norm_plates[[length(norm_plates)+1]] <- i
  }
  
  rm(i,j)
  norm_plates <- do.call(rbind, norm_plates) #Re-form one table of values
  
  norm_plates  
}







## Normalize by same-culture controls using median-based z-score (on the same DIV)
normalize_by_date <- function(all_data) {
  
  # Split data by date and DIV
  per_date_split <- split(all_data, interaction(all_data[,"date"], all_data[,"DIV"], drop=TRUE))
  
  scaled_dates <- list() #initialize list
  
  # Loop through each date at each DIV
  for (i in per_date_split) {
    
    # Loop through each ontogeny parameter
    for (j in names(all_data[,8:(ncol(all_data)-1)])) {
      
      cntrls <- subset(i, trt=="Control")[,j] # Get vector of control values for that DIV on that date
      
      # Check if SD of controls on that date is zero, if not:
      if (sd(cntrls) != 0) {
        i[,j] <- (i[,j] - median(cntrls))/sd(cntrls) # Subtract all values on that DIV and date by median of controls and divide by SD
      
      } else {
        print(paste("The SD for culture ", i[1,"date"]," for ", j, " at DIV", i[1,"DIV"], " is zero, using global z-score instead.", sep=""))
        
        # Get SD from all control wells for that network parameter
        global_sd <- aggregate(get(j) ~ DIV, sd, data=subset(all_data, trt=="Control"))
        
        # Get median from all control wells for that network parameter
        global_med <- aggregate(get(j) ~ DIV, median, data=subset(all_data, trt=="Control"))
        
        all_cntrls_sd <- subset(global_sd, DIV==i[1,"DIV"])[1,2] # Isolate SD for same DIV
        all_cntrls_med <- subset(global_med, DIV==i[1,"DIV"])[1,2] # Isolate median for same DIV
        print(paste("Using global median of ", all_cntrls_med, sep=""))
        print(paste("Using global sd of ", all_cntrls_sd, sep=""))
        
        i[,j] <- (i[,j] - all_cntrls_med)/all_cntrls_sd # Subtract median from same culture, but now divide by global SD 
      }
    
    }
    # Add normalized data to list of dates 
    scaled_dates[[length(scaled_dates)+1]] <- i
  }
  
  rm(i,j)
  scaled_dates <- do.call(rbind, scaled_dates) #Re-form one table of values
  
  scaled_dates    
}







## Apply normalization to the plate level by dividing raw values by plate median, then scaling with z-score across plates
normalize_by_plate <- function(all_data) {
  
  # Split data by plate (date as well because plateIDs get reused) and DIV
  per_plate_split <- split(all_data, interaction(all_data[,"date"], all_data[,"Plate.SN"], all_data[,"DIV"], drop=TRUE))
  
  scaled_plates <- list() #initialize list
  
  # Loop through each plate at each DIV
  for (i in per_plate_split) {
    
    # Loop through each ontogeny parameter
    for (j in names(all_data[,8:(ncol(all_data)-1)])) {
      
      # If all plate values for a parameter are NA, set all fold-changes to 1
      if (length(na.omit(i[,j])) == 0) { 
        print(paste("All NAs for", i[,"Plate.SN"][1], i[,"DIV"][1], j, sep=" - "))
        i[,j] <- 1
      
      # If plate median is still zero, set all fold-changes to 1
      } else if(median(na.omit(i[,j])) == 0) {
        print(paste("All Zeros for", i[,"Plate.SN"][1], i[,"DIV"][1], j, sep=" - "))
        i[,j] <- 1
      
      # If plate median is not zero, divide to get fold-change
      } else {
        i[,j] <- i[,j] / median(na.omit(i[,j])) # Divide all values on that plate by median of the plate
      }      
      
    }
    
    scaled_plates[[length(scaled_plates)+1]] <- i
  }
  
  rm(i,j)
  scaled_plates <- do.call(rbind, scaled_plates) # Re-form one table of values
  
  scaled_plates[is.na(scaled_plates)] <- 1 # Change all remaining NAs to 1
  scaled_plates[scaled_plates[,"R"] < 0, "R"] <- 0 # Set correlation value minimums to zero
  
  # Log-transform fold-change values (adding 1 to set minimum log2() value to 0)
  scaled_plates[,8:(ncol(scaled_plates)-1)] <- log2(scaled_plates[,8:(ncol(scaled_plates)-1)] + 1)
  
   
  ## Z-score scale each parameter over each DIV timepoint
  # First split by DIV
  div_split <- split(scaled_plates, scaled_plates[,"DIV"]) 
  
  zscore_plates <- list() #initialize list
  
  # Loop through DIV
  for (i in div_split) {
    
    # Loop through network params
    for (j in names(scaled_plates[,8:(ncol(scaled_plates)-1)])) {
      
      i[,j] <- (i[,j] - median(i[,j])) / sd(i[,j]) # Calculate z-score
    }
  
  zscore_plates[[length(zscore_plates)+1]] <- i # Add to list
  }

  rm(i,j)
  zscore_plates <- do.call(rbind, zscore_plates) # Re-form one table of values

  zscore_plates
}







## Apply normalization to the plate level by dividing raw values by plate median of CONTROLS ONLY, then scaling with z-score across plates
normalize_by_plate_v2 <- function(all_data) {
  
  # Split data by plate (date as well because plateIDs get reused) and DIV
  per_plate_split <- split(all_data, interaction(all_data[,"date"], all_data[,"Plate.SN"], all_data[,"DIV"], drop=TRUE))
  
  scaled_plates <- list() #initialize list
  
  # Loop through each plate at each DIV
  for (i in per_plate_split) {
    
    # Loop through each ontogeny parameter
    for (j in names(all_data[,8:(ncol(all_data)-1)])) {
      
      cntrls <- (subset(i, trt=="Control")[,j]) # Get vector of control values for that DIV on that plate
      
      # If all plate values for a parameter are NA, set all fold-changes to 1
      if (length(na.omit(cntrls)) == 0) {
        print(paste("All NAs for", i[,"Plate.SN"][1], i[,"DIV"][1], j, sep=" - "))
        i[,j] <- 1
        
        # If plate median is still zero, set all fold-changes to 1
      } else if (median(na.omit(cntrls)) == 0) {
        print(paste("All Zeros for", i[,"Plate.SN"][1], i[,"DIV"][1], j, sep=" - "))
        i[,j] <- 1
        
        # If plate median is not zero, divide to get fold-change
      } else {
        i[,j] <- i[,j] / median(na.omit(cntrls)) # Divide all values on that plate by median of the plate
      }      
      
    }
    
    scaled_plates[[length(scaled_plates)+1]] <- i
  }
  
  rm(i,j)
  scaled_plates <- do.call(rbind, scaled_plates) # Re-form one table of values
  
  scaled_plates[is.na(scaled_plates)] <- 1 # Change all remaining NAs to 1
  scaled_plates[scaled_plates[,"R"] < 0, "R"] <- 0 # Set correlation value minimums to zero
  
  # Log-transform fold-change values
  scaled_plates[,8:(ncol(scaled_plates)-1)] <- log2(scaled_plates[,8:(ncol(scaled_plates)-1)] + 1)
  
  
  ## Z-score scale each parameter over each DIV timepoint
  # First split by DIV
  div_split <- split(scaled_plates, scaled_plates[,"DIV"]) 
  
  zscore_plates <- list() #initialize list
  
  # Loop through DIV
  for (i in div_split) {
    
    # Loop through network params
    for (j in names(scaled_plates[,8:(ncol(scaled_plates)-1)])) {
      
      i[,j] <- (i[,j] - median(i[,j])) / sd(i[,j]) # Calculate z-score
    }
    
    zscore_plates[[length(zscore_plates)+1]] <- i # Add to list
  }
  
  rm(i,j)
  zscore_plates <- do.call(rbind, zscore_plates) # Re-form one table of values
  
  zscore_plates
}




## Apply normalization to the plate level by dividing raw values by plate median of controls only, then scaling with z-score across plates
## The change from v2 above is now computing z-score based on control well distribution only (i.e. the z-score reflects number of SD away from control median)
normalize_by_plate_v3 <- function(all_data) {
  
  # Split data by plate (date as well because plateIDs get reused) and DIV
  per_plate_split <- split(all_data, interaction(all_data[,"date"], all_data[,"Plate.SN"], all_data[,"DIV"], drop=TRUE))
  
  scaled_plates <- list() #initialize list
  
  # Loop through each plate at each DIV
  for (i in per_plate_split) {
    
    # Loop through each ontogeny parameter
    for (j in names(all_data[,8:(ncol(all_data)-1)])) {
      
      cntrls <- (subset(i, trt=="Control")[,j]) # Get vector of control values for that DIV on that plate
      
      # If all plate values for a parameter are NA, set all fold-changes to 1
      if (length(na.omit(cntrls)) == 0) {
        print(paste("All NAs for", i[,"Plate.SN"][1], i[,"DIV"][1], j, sep=" - "))
        i[,j] <- 1
        
        # If plate median is still zero, set all fold-changes to 1
      } else if (median(na.omit(cntrls)) == 0) {
        print(paste("All Zeros for", i[,"Plate.SN"][1], i[,"DIV"][1], j, sep=" - "))
        i[,j] <- 1
        
        # If plate median is not zero, divide to get fold-change
      } else {
        i[,j] <- i[,j] / median(na.omit(cntrls)) # Divide all values on that plate by median of the plate
      }      
      
    }
    
    scaled_plates[[length(scaled_plates)+1]] <- i
  }
  
  rm(i,j)
  scaled_plates <- do.call(rbind, scaled_plates) # Re-form one table of values
  
  scaled_plates[is.na(scaled_plates)] <- 1 # Change all remaining NAs to 1
  scaled_plates[scaled_plates[,"R"] < 0, "R"] <- 0 # Set correlation value minimums to zero
  
  # Log-transform fold-change values
  scaled_plates[,8:(ncol(scaled_plates)-1)] <- log2(scaled_plates[,8:(ncol(scaled_plates)-1)] + 1)
  
  
  ## Z-score scale each parameter over each DIV timepoint
  # First split by DIV
  div_split <- split(scaled_plates, scaled_plates[,"DIV"]) 
  
  zscore_plates <- list() #initialize list
  
  # Loop through DIV
  for (i in div_split) {
    
    # Loop through network params
    for (j in names(scaled_plates[,8:(ncol(scaled_plates)-1)])) {
      
      cntrl_vals <- subset(i, trt=="Control")[,j] # Isolate control values for that parameter
      i[,j] <- (i[,j] - mean(cntrl_vals)) / sd(cntrl_vals) # Calculate z-score
    }
    
    zscore_plates[[length(zscore_plates)+1]] <- i # Add to list
  }
  
  rm(i,j)
  zscore_plates <- do.call(rbind, zscore_plates) # Re-form one table of values
  
  zscore_plates
}







## Separate out data for compound of interest and plot all network parameter median z-scores
## With cone of uncertainty (3 times the median absolute deviation)
plot_dose_perturbations <- function(norm_data, parameter, compound) {
  
  # Find median of replicate response on each DIV
  data_meds <- as.data.frame(do.call(cbind, aggregate(get(parameter) ~ trt + dose + DIV, FUN=function(x) c(med=median(x), mad=mad(x), n=length(x), low=(median(x) - (sd(x)/sqrt(length(x)))), high=(median(x) + (sd(x)/sqrt(length(x))))), data=norm_data)), stringsAsFactors=FALSE)
  
  # Subset to compound of interest and control values
  trt_meds <- subset(data_meds, trt %in% c(compound, "Control"))
  
  # Add in DIV0 placeholder data to anchor plots at zero
  trt_meds <- rbind(data.frame(trt="Compound", dose=as.numeric(as.character(unique(trt_meds[,"dose"]))), DIV=0, med=0, mad=0, n=1, low=0, high=0), trt_meds)
  trt_meds[,2:ncol(trt_meds)] <- lapply(trt_meds[,2:ncol(trt_meds)], as.numeric)
  
  # Isolate control standard deviations
  trt_cntrls <- subset(trt_meds, trt=="Control")
  trt_cntrls <- rbind(data.frame(trt="Control", dose=0, DIV=0, med=0, mad=0, n=1, low=0, high=0), trt_cntrls)
  
  # Repeat control values same number of times as total number of samples
  trt_cntrls2 <- trt_cntrls[rep(seq_len(nrow(trt_cntrls)), nrow(trt_meds)/nrow(trt_cntrls)),]
  
  # Change low and high values to three times the median absolute deviation
  trt_cntrls2[,"low"] <- trt_cntrls2[,"med"] - 3*trt_cntrls2[,"mad"]
  trt_cntrls2[,"high"] <- trt_cntrls2[,"med"] + 3*trt_cntrls2[,"mad"]
  
  # Prepare SEM intervals
  limits <- aes(ymin=trt_meds[,"low"], ymax=trt_meds[,"high"])
  
  # Plot with ggplot
  ggplot(trt_meds, aes(x=DIV, y=med, group=dose, color=log(dose+1))) +
    geom_ribbon(ymin=trt_cntrls2$low, ymax=trt_cntrls2$high, x=trt_cntrls2$DIV, fill="grey70", alpha=0.8, inherit.aes=FALSE) +  
    geom_line(size=0.4, position=position_jitter(w=0, h=0.1)) +
    # geom_errorbar(limits, size=0.2, width=0.25) +  # Uncomment to add SEM error bars
    scale_color_gradient(low="blue", high="red") + 
    coord_cartesian(ylim=c(-5,5)) + 
    theme(axis.text=element_text(size=4), axis.title.x=element_text(size=6), axis.title.y=element_text(size=0), plot.title=element_text(size=6), legend.text=element_text(size=2), legend.position = "none") +
    ggtitle(paste(compound, parameter, sep=" - "))
}







## Multiple plot function (sourced from http://peterhaschke.com/Code/multiplot.R)
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}






## Plot all parameters together for a given compound
make_multi <- function(compound) {
  
  p <- list() #initialize plot list
  i <- 1 # intialize counter
  
  # loop through all network parameters
  for (j in names(all_data[8:(ncol(all_data)-3)])){
    
    p[[i]] <- plot_dose_perturbations(z_norm_data, j, compound) # generate plot for each parameter
    i <- i + 1
  }
  
  # write multiplot to pdf file
  pdf(paste(compound, "_endpoints.pdf", sep=""))
  multiplot(plotlist=p, cols=4)
  dev.off() 
  
}






## Calculate scalar measure of total network parameter perturbation
scal_perturb <- function(norm_data, compound, slopes=FALSE) {
  
  # Loop through each network parameter to get median values across replicates
  # Choosing parameters among the 16 based on "set_parameters" character vector
  data_meds <- list()
  for (j in set_parameters) {
    # Find median and MAD of replicate responses on each DIV
    data_meds[[j]] <- as.data.frame(do.call(cbind, aggregate(get(j) ~ trt + dose + DIV, FUN=function(x) c(med=median(x), mad=3*mad(x), n=length(x)), data=norm_data)), stringsAsFactors=FALSE)
  }
  
  data_meds <- do.call(rbind, data_meds) # Reform dataframe
  data_meds[,2:ncol(data_meds)] <- lapply(data_meds[,2:ncol(data_meds)], as.numeric)
  
  # Aggregate across all network parameters to obtain sum (Euclidean norm) of absolute valued Zscore medians
  data_sums <- as.data.frame(do.call(cbind, aggregate(med ~ trt + dose + DIV, FUN=function(x) c(sum=sqrt(sum(x^2)), n=length(x)), data=data_meds)), stringsAsFactors=FALSE)
  #data_sums <- as.data.frame(do.call(cbind, aggregate(med ~ trt + dose + DIV, FUN=function(x) c(sum=sum(abs(x)), n=length(x)), data=data_meds)), stringsAsFactors=FALSE)
  
  # Subset to compound of interest and control values
  trt_sums <- subset(data_sums, trt %in% c(compound, "Control"))
  
  # Add in DIV0 placeholder data to anchor plots at zero
  trt_sums <- rbind(data.frame(trt="Compound", dose=as.numeric(as.character(unique(trt_sums[,"dose"]))), DIV=0, sum=0, n=1), trt_sums)
  trt_sums[,2:ncol(trt_sums)] <- lapply(trt_sums[,2:ncol(trt_sums)], as.numeric)
  #print(trt_sums)
  
  # Calculate sum (Euclidean norm) of control well 3*MAD variability
  cntrl_var <- as.data.frame(do.call(cbind, aggregate(mad ~ DIV, FUN=function(x) c(dose=0, sum=sqrt(sum(x^2)), n=length(x)), data=subset(data_meds, trt=="Control"))), stringsAsFactors=FALSE)
  cntrl_var <- as.data.frame(lapply(cntrl_var[,1:ncol(cntrl_var)], as.numeric))
  cntrl_var[,"low"] <- 0
  cntrl_var[,"high"] <- cntrl_var[,"sum"]
  
  # Add in DIV0 value and Repeat control values same number of times as total number of samples (minus control rows)
  cntrl_var <- rbind(data.frame(DIV=0, dose=0, sum=0, n=1, low=0, high=0), cntrl_var)
  cntrl_var <- cntrl_var[rep(seq_len(nrow(cntrl_var)), nrow(subset(trt_sums, dose!=0))/nrow(cntrl_var)),]
  
  # Find cytotox viability data and get mean value (of medians) for each dose
  via_meds <- aggregate(median ~ dose, mean, data=subset(all_via, trt==compound))
  via_meds[,"DIV"] <- 12.5
  via_meds[,"sum"] <- (100 - via_meds[,"median"])/(100/8) # Change % alive to % dead ; Divide to rescale to y-axis
  via_meds$sum[via_meds$sum < 0] <- 0 # Set lower limit of 0%
  via_meds <- data.frame(trt=compound, dose=via_meds[,"dose"], DIV=via_meds[,"DIV"], sum=via_meds[,"sum"])
  #print(via_meds)

  # Plot with ggplot2
  vis <- ggplot(subset(trt_sums, dose!=0), aes(x=DIV, y=sum, group=dose, color=log(dose))) +
        geom_ribbon(ymin=cntrl_var$low, ymax=cntrl_var$high, x=cntrl_var$DIV, fill="grey70", alpha=0.8, inherit.aes=FALSE) +
        geom_smooth(size=1) +
        scale_color_gradient(low="blue", high="red") +
        coord_cartesian(ylim=c(0,8)) +
        ggtitle(compound) +
        theme(axis.text=element_text(size=14), axis.title.x=element_text(size=16), axis.title.y=element_text(size=0), plot.title=element_text(size=18), legend.text=element_text(size=14)) 
  
  vis <- vis + geom_point()
  vis <- vis + geom_point(data=subset(via_meds, dose!=0), position=position_jitter(w=0.25, h=0)) # Adds viability points
  
  # Print plot
  print(vis)
  
  
  if (slopes==TRUE) {
    ts_split <- split(trt_sums, trt_sums[,"dose"]) # Split on different doses
    ts_slopes <- list() # initialize list of slopes
    
    # Loop through each dose, calculate slope between adjacent timepoints
    for (i in 1:length(ts_split)) {
      ts_slopes[[i]] <- diff(ts_split[[i]]$sum) / diff(ts_split[[i]]$DIV) # Diff takes difference between each element of the vector and the previous element
      print(ts_split[[i]]) # Uncomment to look at slope calculations
    }
    
    ts_slopes <- do.call(rbind, ts_slopes) # Put all concentrations together
    ts_slopes <- data.frame(cbind(unique(trt_sums[,"dose"]), ts_slopes), stringsAsFactors=FALSE) # Form dataframe with dose information in first column
    names(ts_slopes) <- c("dose","div5","div7","div9","div12")
    
    # Add column of mean slope from DIV7 to DIV12
    ts_slopes[,"slp"] <- rowMeans(ts_slopes[,4:5])
    
    # Calculate change in slope over change in dose (derivative)
    ts_slopes[1,"slp"] <- 0  # Setting zero dose (control wells) slope to zero for derivative calculations
    ts_slopes[,"deriv"] <- c(0, diff(ts_slopes[,"slp"]) / diff(ts_slopes[,"dose"]))
    
    print(ts_slopes)
    print(c(0, diff(ts_slopes[,"slp"])))
    print(diff(ts_slopes[,"dose"]))
    
    # Remove zero dose row and log-transform doses
    ts_slopes <- subset(ts_slopes, dose!=0)
    ts_slopes[,"dose"] <- log(ts_slopes[,"dose"])
    
    # Multiplying slope and derivative by constant for plotting on same scale as scalar perturbation
    ts_slopes[,"slp"] <- ts_slopes[,"slp"]*5
    ts_slopes[,"deriv"] <- ts_slopes[,"deriv"]*5
    
    # Setting first point of derivative to zero to reflect no confidence in lower dose than tested slope
    ts_slopes[1,"deriv"] <- 0
    
    # Gather endpoint scalar perturbation values at DIV12
    trt_sums_ep <- subset(trt_sums, DIV==12 & dose!=0)
    trt_sums_ep[,"dose"] <- log(trt_sums_ep[,"dose"])
    
    # Print tables used for slope plot
    print(trt_sums_ep)
    print(ts_slopes)
    
    # Create loess-smoothened lines
    scal_loess <- loess(sum ~ dose, data=trt_sums_ep, span=1) 
    vel_loess <- loess(slp ~ dose, data=ts_slopes, span=1)
    dv_loess <- loess(deriv ~ dose, data=ts_slopes, span=1)
    
    # Plot slopes at DIV9 to DIV12 over different doses of compound
    plot(ts_slopes[,"dose"], ts_slopes[,"slp"], ylim=c(-10,10), xlab="log(concentration)", ylab="", pch=19, col="darkblue", main=compound, cex.lab=1.5, cex.main=1.6, cex.axis=1.2)
    abline(h=0, col="darkgray", lty=2, lwd=1) # add y=0 line for reference
    abline(h=(subset(cntrl_var, DIV==12)[1,"high"]), lty=2, lwd=2, col="darkgray") # add y=DIV12 3*MAD
    
    lines(ts_slopes[,"dose"], vel_loess$fitted, lwd=2, col="darkblue") # plot mean slope between DIV7 and DIV12 (change in perturbation over time)
    
    points(sum ~ dose, data=trt_sums_ep, pch=19, col="darkgreen")
    lines(trt_sums_ep[,"dose"], scal_loess$fitted, lwd=2, col="darkgreen") # plot scalar perturbation over dose at DIV12
    
    points(ts_slopes[,"dose"], ts_slopes[,"deriv"], pch=19, col="red")
    lines(ts_slopes[,"dose"], dv_loess$fitted, lwd=2, col="red") # plot derivative of DIV7 to DIV12 slope (slope over change in dose)
    
    legend("bottomright", legend=c("|X|", "V=dX/dt", "dV/dC"), col=c("darkgreen", "darkblue", "red"), lwd=2)
    
  }
}






## Calculate scalar measure of total network parameter perturbation for EACH REPLICATE
scal_perturb_reps <- function(norm_data, compound) {
  
  # Subset normalized data to compound of interest
  comp_norm <- subset(norm_data, trt==compound)
  
  # Split compound data by date and plate ID to separate replicates
  cnd_reps <- split(comp_norm, interaction(comp_norm[,"date"], comp_norm[,"Plate.SN"], comp_norm[,"well"], drop=TRUE))
  
  # Loop through replicates to compute Euclidean norm
  trt_sums <- list()
  trt_derivs <- list()
  n <- 1
  for (i in cnd_reps) {
    
    # Subset to set_parameters
    cnd_reps_sp <- cbind(i[,c("date","Plate.SN","DIV","well","trt","dose","units")], i[,set_parameters])
    
    # Calculate Euclidean norm across rows
    cnd_reps_sp[,"sum"] <- apply(cnd_reps_sp[,8:ncol(cnd_reps_sp)], 1, function(x) sqrt(sum(x^2)))
    
    # Add data frame to list
    trt_sums[[n]] <- cnd_reps_sp 
    
    # Isolate sums and calculate slopes
    trt_slopes <- cnd_reps_sp[,c("date","Plate.SN","DIV","well","trt","dose","sum")]
    trt_slopes[,"div5to12"] <- (subset(trt_slopes, DIV==12)[,"sum"][1] - subset(trt_slopes, DIV==5)[,"sum"][1])/7
    trt_slopes[,"div7to12"] <- (subset(trt_slopes, DIV==12)[,"sum"][1] - subset(trt_slopes, DIV==7)[,"sum"][1])/5
    trt_slopes[,"div9to12"] <- (subset(trt_slopes, DIV==12)[,"sum"][1] - subset(trt_slopes, DIV==9)[,"sum"][1])/3
    trt_slopes <- trt_slopes[nrow(trt_slopes),]
    
    # Add data frame of slopes to second list
    trt_derivs[[n]] <- trt_slopes   
    
    # Increment n
    n <- n + 1
  }  
  
  # Convert list of replicates to one dataframe
  trt_sums <- do.call(rbind, trt_sums)
  trt_derivs <- do.call(rbind, trt_derivs)
  
  # Subset to needed data only
  trt_sums <- trt_sums[,c("date","Plate.SN","DIV","well","trt","dose","sum")]
  
  # Add in DIV0 placeholder data to anchor plots at zero
  trt_sums <- rbind(data.frame(date=0, Plate.SN=0, DIV=0, well="AO", trt="Compound", dose=as.numeric(as.character(unique(trt_sums[,"dose"]))), sum=0), trt_sums)
  #print(trt_sums)
  
  
  
  ## Now find 3*MAD variability of control wells
  # Get control well data across experiment
  cntrl_norm <- subset(norm_data, trt=="Control")
  
  # Split compound data by date, plate ID, and well to separate replicates
  cntrl_reps <- split(cntrl_norm, interaction(cntrl_norm[,"date"], cntrl_norm[,"Plate.SN"], cntrl_norm[,"well"], drop=TRUE))
  
  # Loop through control well replicates to compute Euclidean norm
  cntrl_sums <- list()
  n <- 1
  for (i in cntrl_reps) {
    
    # Subset to set_parameters
    cntrl_reps_sp <- cbind(i[,c("date","Plate.SN","DIV","well","trt","dose","units")], i[,set_parameters])
    
    # Calculate Euclidean norm across rows
    cntrl_reps_sp[,"sum"] <- apply(cntrl_reps_sp[,8:ncol(cntrl_reps_sp)], 1, function(x) sqrt(sum(x^2)))
    
    # Add data frame to list
    cntrl_sums[[n]] <- cntrl_reps_sp 
    n <- n + 1
  } 
  
  # Convert list of replicates to one dataframe
  cntrl_sums <- do.call(rbind, cntrl_sums)
  
  # Aggregate to find 3*MAD across DIV timepoints
  cntrl_var <- as.data.frame(do.call(cbind, aggregate(sum ~ DIV, FUN=function(x) c(med=median(x), mad=3*mad(x), n=length(x)), data=cntrl_sums)), stringsAsFactors=FALSE)
  
  # Add min and max values for plotting
  cntrl_var[,"min"] <- 0
  cntrl_var[,"max"] <- cntrl_var[,"mad"] + cntrl_var[,"med"] # adding the 3*MAD to the median perturbation of controls
  
  # Add in DIV0 placeholder value 
  cntrl_var <- rbind(data.frame(DIV=0, med=0, mad=0, n=1, min=0, max=0), cntrl_var)
  
  # Get DIV12 3*MAD for derivative plotting
  final_var <- subset(cntrl_var, DIV==12)[,"max"][1]
  
  
  ## Get cytotoxicity data associated with compound
  # Get mean value (of medians) for each dose
  via_meds <- aggregate(median ~ dose, mean, data=subset(all_via, trt==compound))
  via_meds[,"DIV"] <- 14
  via_meds[,"sum"] <- (100 - via_meds[,"median"])/(100/20) # Change % alive to % dead ; Divide to rescale to y-axis
  via_meds$sum[via_meds$sum < 0] <- 0 # Set lower limit of 0%
  via_meds <- data.frame(trt=compound, dose=via_meds[,"dose"], DIV=via_meds[,"DIV"], sum=via_meds[,"sum"])
  
  
  ## Plot scalar perturbs with ggplot2
  vis <- ggplot(trt_sums, aes(x=DIV, y=sum, group=dose, color=log(dose))) +
    geom_ribbon(data=cntrl_var, ymin=cntrl_var[,"min"], ymax=cntrl_var[,"max"], x=cntrl_var[,"DIV"], fill="grey70", alpha=0.8, linetype=2, inherit.aes=FALSE) +
    geom_smooth(size=1.4, se=FALSE, span=1) +
    geom_point(position=position_jitter(w=0.25, h=0.2)) +
    scale_color_gradient(low="blue", high="red") +
    coord_cartesian(ylim=c(0,20)) +
    ggtitle(compound) +
    scale_x_continuous(breaks = c(0,2,5,7,9,12)) +
    theme(axis.text=element_text(size=14), axis.title.x=element_text(size=16), axis.title.y=element_text(size=0), plot.title=element_text(size=18), legend.text=element_text(size=14)) +
    geom_point(data=subset(via_meds, dose!=0), shape=17, position=position_jitter(w=0.25, h=0)) # Adds viability points
  
  print(vis)
  
  
  # Plot DIV9 to DIV12 change
  vis2 <- ggplot(trt_derivs, aes(x=log(dose))) +
    geom_hline(yintercept=final_var, linetype=2) +
    geom_hline(yintercept=0, linetype=2) +
    geom_smooth(aes(y=sum), size=1.4, se=TRUE, span=1, col="blue", fill="blue") +
    geom_smooth(aes(y=div9to12*5), size=1.4, se=TRUE, span=1, col="forestgreen", fill="forestgreen") +
    geom_point(aes(y=sum), col="blue", position=position_jitter(w=0.2, h=0.25)) +
    geom_point(aes(y=div9to12*5), col="forestgreen", position=position_jitter(w=0.2, h=0.25)) +
    coord_cartesian(ylim=c(-20,20)) +
    ggtitle(paste(compound, " DIV9-12", sep="")) +
    theme(axis.text=element_text(size=14), axis.title.x=element_text(size=16), axis.title.y=element_text(size=0), plot.title=element_text(size=18), legend.text=element_text(size=14)) 
  
  
  # Plot DIV7 to DIV12 change
  vis3 <- ggplot(trt_derivs, aes(x=log(dose))) +
    geom_hline(yintercept=final_var, linetype=2) +
    geom_hline(yintercept=0, linetype=2) +
    geom_smooth(aes(y=sum), size=1.4, se=TRUE, span=1, col="blue", fill="blue") +
    geom_smooth(aes(y=div7to12*5), size=1.4, se=TRUE, span=1, col="forestgreen", fill="forestgreen") +
    geom_point(aes(y=sum), col="blue", position=position_jitter(w=0.2, h=0.25)) +
    geom_point(aes(y=div7to12*5), col="forestgreen", position=position_jitter(w=0.2, h=0.25)) +
    coord_cartesian(ylim=c(-20,20)) +
    ggtitle(paste(compound, " DIV7-12", sep="")) +
    theme(axis.text=element_text(size=14), axis.title.x=element_text(size=16), axis.title.y=element_text(size=0), plot.title=element_text(size=18), legend.text=element_text(size=14)) 

  
  # Plot DIV5 to DIV12 change
  vis4 <- ggplot(trt_derivs, aes(x=log(dose))) +
    geom_hline(yintercept=final_var, linetype=2) +
    geom_hline(yintercept=0, linetype=2) +
    geom_smooth(aes(y=sum), size=1.4, se=TRUE, span=1, col="blue", fill="blue") +
    geom_smooth(aes(y=div5to12*5), size=1.4, se=TRUE, span=1, col="forestgreen", fill="forestgreen") +
    geom_point(aes(y=sum), col="blue", position=position_jitter(w=0.2, h=0.25)) +
    geom_point(aes(y=div5to12*5), col="forestgreen", position=position_jitter(w=0.2, h=0.25)) +
    coord_cartesian(ylim=c(-20,20)) +
    ggtitle(paste(compound, " DIV5-12", sep="")) +
    theme(axis.text=element_text(size=14), axis.title.x=element_text(size=16), axis.title.y=element_text(size=0), plot.title=element_text(size=18), legend.text=element_text(size=14)) 
            
  
  # Find derivative of slope function
  add_deriv <- function(dv_plot) {
    # Get change in scalar over time loess fitted values
    # The fourth function in the vis2-4 ggplot calls is the slope
    x_vals <- ggplot_build(dv_plot)$data[[4]]$x
    y_vals <- ggplot_build(dv_plot)$data[[4]]$y
    y_mins <- ggplot_build(dv_plot)$data[[4]]$ymin
    y_maxs <- ggplot_build(dv_plot)$data[[4]]$ymax
    
    # Approximate derivative with change in y over change in x
    xy_prime <- diff(y_vals)/diff(x_vals)
    xy_prime_min <- diff(y_mins)/diff(x_vals)
    xy_prime_max <- diff(y_maxs)/diff(x_vals)
    deriv <- data.frame(x=x_vals[2:length(x_vals)], y=xy_prime, ymin=xy_prime_min, ymax=xy_prime_max)
    
    # Add derivative to ggplot object
    dv_plot <- dv_plot + geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), data=deriv, inherit.aes=FALSE, fill="red", alpha=0.5) +
      geom_line(aes(x=x, y=y), data=deriv, inherit.aes=FALSE, col="red", size=1.4)
    
    # If median scalar perturb at highest dose exceeds 3*MAD of controls,
    # And at least one y-value of the velocity curve is negative (recovering trajectory)
    # And at least one y-value of the derivative is distinctly positive,
    # Find first dose (x-value) where derivative crosses y=0 line
    if (median(subset(trt_derivs, dose==max(trt_derivs[,"dose"]))[,"sum"]) > final_var && min(y_vals[1:40]) < 0 && max(deriv[,"y"]) > 0.001) {
      
      # If the highest dose derivative values are negative, remove them 
      # This represents topping out of scalar perturbation, which is not of interest 
      # Use cutoff slightly greater than zero to also catch practically zero slopes
      if (deriv[nrow(deriv),"y"] < 0.001) {
        last_pos <- max(which(deriv[,"y"] > 0.001)) # index of last positive value
        deriv <- deriv[1:last_pos,] # subset to remove last string of negative values
      }
      
      # If there is still a negative derivative value, find index of last negative and estimate tipping point
      if (min(deriv[,"y"]) < 0) {
        last_neg <- max(which(deriv[,"y"] < 0)) # index of last negative value
        root1 <- round(exp(deriv[last_neg+1,"x"]), 3) # Get dose of next (positive) value
        cat(paste(compound, root1, sep=','), "\n") 
        
      # Otherwise, the first dose produced a positive slope so set the root to lowest dose  
      } else {
        root1 <- min(trt_derivs[,"dose"])
        cat(paste(compound, root1, sep=','), "\n")
      }
      
      # Add tipping point estimate to ggplot object
      dv_plot <- dv_plot + geom_vline(xintercept=log(root1), linetype=2)
    
    # Otherwise, we could not define a tipping point concentration
    } else {
      root1 <- NA
      cat(paste(compound, root1, sep=','), "\n")
    }
    
    # Print plot      
    print(dv_plot)
  }
  
  # write results to file
  sink(file="tipping_point_list.csv", append=TRUE, split=TRUE)
  add_deriv(vis2)
  add_deriv(vis3)
  add_deriv(vis4)
  sink()
}
  
  
  


  
  


## Calculate scalar measure of total network parameter perturbation based on iterative random sampling of network parameters
rand_scal_perturb <- function(norm_data, compound, samples=100) {
  
  # Subset to compound of interest and control values
  comp_norm_data <- subset(norm_data, trt %in% c(compound, "Control"))
  
  # Loop through each network parameter and get median values across replicates
  param_names <- c("MFR","BR","ISI","%SiB","BD","IBI","#AE","#ABE","#NS","NSP","NSD","%SiNS","NS-ISI","NSDsd","#SiNS","R")
  param_weights <- c(0.8907,0.8869,0.9830,2.2564,0.8014,0.7904,3.5272,2.6909,0.6294,1.7568,1.2278,0.9404,0.5686,0.8481,0.8288,1.0627) # These are based on CV of DIV9-12 control wells
  #param_weights <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1) # Uncomment to set all parameters to equal probability of being selected
  
  data_meds <- list()
  for (j in param_names) {
    # Find median and MAD of replicate responses on each DIV
    data_meds[[j]] <- as.data.frame(do.call(cbind, aggregate(get(j) ~ trt + dose + DIV, FUN=function(x) c(med=median(x), mad=3*mad(x), n=length(x)), data=comp_norm_data)), stringsAsFactors=FALSE)
  }
  
  # Iteratively select subset of network parameters to sum up
  all_trt_sums <- list()
  for (i in 1:samples) {
    
    set_parameters <- sample(param_names, 8, replace=FALSE, prob=param_weights) # randomly sample 8 parameters
    data_meds_sub <- data_meds[set_parameters] # Select list elements matching 8 parameter names
    data_meds_sub <- do.call(rbind, data_meds_sub) # Convert list to dataframe
    data_meds_sub[,2:ncol(data_meds_sub)] <- lapply(data_meds_sub[,2:ncol(data_meds_sub)], as.numeric) # Make sure its numeric

    # Aggregate across all network parameters to obtain sum (Euclidean norm) of absolute valued Zscore medians
    trt_sums <- as.data.frame(do.call(cbind, aggregate(med ~ trt + dose + DIV, FUN=function(x) c(sum=sqrt(sum(x^2)), n=i), data=data_meds_sub)), stringsAsFactors=FALSE)
  
    # Add in DIV0 placeholder data to anchor plots at zero
    trt_sums <- rbind(data.frame(trt="Compound", dose=as.numeric(as.character(unique(trt_sums[,"dose"]))), DIV=0, sum=0, n=i), trt_sums)
    trt_sums[,2:ncol(trt_sums)] <- lapply(trt_sums[,2:ncol(trt_sums)], as.numeric)
    
    # Add this sample of parameters to growing list
    all_trt_sums[[i]] <- trt_sums
  }  
    
  all_trt_sums <- do.call(rbind, all_trt_sums) # Convert list to dataframe
  
  # Split up data by DIV and dose to calculate CI limits of each
  ats_split <- split(all_trt_sums, interaction(all_trt_sums[,"dose"], all_trt_sums[,"DIV"]))
  
  ats_sum <- data.frame() # initialize new data frame
  for (tmt in ats_split){
    quant_vals <- quantile(tmt[,"sum"], probs=c(0.05,0.5,0.95)) # calculate CI
    ats_sum <- rbind(ats_sum, data.frame(trt=tmt[,"trt"][1], dose=tmt[,"dose"][1], DIV=tmt[,"DIV"][1], sum=quant_vals, n=c(1,2,3)))
  }
  
  ats_low <- subset(ats_sum, n==1) # X% of data below this value
  ats_med <- subset(ats_sum, n==2) # median value
  ats_high <- subset(ats_sum, n==3) # X% of data above this value
  ats_ci <- cbind(ats_med, ats_low[,"sum"], ats_high[,"sum"]) # Combining low and high for geom_ribbon call below
  names(ats_ci) <- c("trt","dose","DIV","sum","n","low","high")
  
  # Find cytotox viability data and get mean value (of medians) for each dose
  via_meds <- aggregate(median ~ dose, mean, data=subset(all_via, trt==compound))
  via_meds[,"DIV"] <- 12.5
  via_meds[,"sum"] <- (100 - via_meds[,"median"])/(100/8) # Change % alive to % dead ; Divide by value to rescale to plot y-axis
  via_meds$sum[via_meds$sum < 0] <- 0 # Set lower limit of 0%
  via_meds <- data.frame(trt=compound, dose=via_meds[,"dose"], DIV=via_meds[,"DIV"], sum=via_meds[,"sum"], n=1)
  
  # Plot with ggplot2
  vis <- ggplot(ats_med, aes(x=DIV, y=sum, group=interaction(dose,n), color=log(dose+1))) +
    geom_ribbon(ymin=ats_ci$low, ymax=ats_ci$high, x=ats_ci$DIV, fill="grey70", alpha=0.3) +
    geom_line(size=1) +
    geom_line(data=ats_low, linetype="dotted", size=1, alpha=0.5) +
    geom_line(data=ats_high, linetype="dotted", size=1, alpha=0.5) +
    geom_point() +
    scale_color_gradient(low="blue", high="red") +
    coord_cartesian(ylim=c(0,8)) +
    ggtitle(compound) +
    theme(axis.text=element_text(size=14), axis.title.x=element_text(size=16), axis.title.y=element_text(size=0), plot.title=element_text(size=18), legend.text=element_text(size=14)) 
  
  vis <- vis + geom_point(data=via_meds, position=position_jitter(w=0.25, h=0)) # Adds viability points
  
  # Print plot
  print(vis)
}










## Look at variability between control wells (one line per plate of controls)
scal_perturb_cntrls <- function(norm_data) {
  
  data_meds <- list()
  
  # Loop through each network parameter to get median of controls for each plate
  # Choosing parameters among the 16 based on "set_parameters" character vector
  for (j in set_parameters) {
    # Find median of replicate responses on each DIV
    data_meds[[j]] <- as.data.frame(do.call(cbind, aggregate(get(j) ~ Plate.SN + DIV, FUN=function(x) c(med=median(x), sd=sd(x), n=length(x)), data=subset(norm_data, trt=="Control"))), stringsAsFactors=FALSE)
  }
  
  data_meds <- do.call(rbind, data_meds) # Reform dataframe
  data_meds[,2:ncol(data_meds)] <- lapply(data_meds[,2:ncol(data_meds)], as.numeric)
  
  # Aggregate across all network parameters to obtain sum of absolute valued Zscore medians
  data_sums <- as.data.frame(do.call(cbind, aggregate(med ~ Plate.SN + DIV, FUN=function(x) c(sum=sum(abs(x)), n=length(x)), data=data_meds)), stringsAsFactors=FALSE)
  
  # Add in DIV0 placeholder data to anchor plots at zero
  trt_sums <- rbind(data.frame(Plate.SN=(unique(data_sums[,"Plate.SN"])), DIV=0, sum=0, n=1), data_sums)
  trt_sums[,2:ncol(trt_sums)] <- lapply(trt_sums[,2:ncol(trt_sums)], as.numeric)
  #print(trt_sums)
  
  
  # Plot with ggplot2
  # Note this includes jitter to keep lines from covering each other up
  vis <- ggplot(trt_sums, aes(x=DIV, y=sum, group=Plate.SN)) + 
    geom_line(size=1) +
    geom_point() +
    coord_cartesian(ylim=c(0,50)) +
    ggtitle("Control Medians") +
    theme(axis.text=element_text(size=14), axis.title.x=element_text(size=16), axis.title.y=element_text(size=0), plot.title=element_text(size=18), legend.text=element_text(size=14))
  
  # Print plot
  print(vis)
  
}











## function to examine individual replicate data (across different plateIDs)
replicate_perturb <- function(norm_data, compound) {
  
  trt_vals <- subset(norm_data, trt==compound) # subset normalized data to compound of interest
  trt_vals <- split(trt_vals, trt_vals[,"Plate.SN"]) # split by plate
  
  reps <- list() # initialize lists and counter
  rep_plots <- list()
  p <- 1
  
  # Loop through each plate to sum up 
  for (i in trt_vals) {
    
    # Calculate network perturb sums for set parameters
    reps[[p]] <- cbind(rowSums(abs(subset(i, select=set_parameters))), i)
    names(reps[[p]])[1] <- "sum"
    
    # Generate plot for this replicate
    rep_plots[[p]] <- ggplot(reps[[p]], aes(x=DIV, y=sum, group=dose, color=log(dose+1))) + geom_line(size=1,  position=position_jitter(w=0, h=0.1)) + scale_color_gradient(low="blue", high="red") + coord_cartesian(ylim=c(0,50)) + ggtitle(paste(compound, reps[[p]]$Plate.SN[1], sep=" ")) + theme(axis.text=element_text(size=12), axis.title.x=element_text(size=14), axis.title.y=element_text(size=0), plot.title=element_text(size=12), legend.text=element_text(size=12), legend.position = "none")
    
    p <- p+1 #increment p
  }
  
  # Plot replicates together
  multiplot(plotlist=rep_plots, cols=2)
}










## Generate outputs

if (set_output == "zscore_plots") {
  
  # normalize data by median subtracted z-score
  if (set_normalization == "plate") {
    z_norm_data <- normalize_by_plate_v3(all_data)
  } else {
    z_norm_data <- normalize_by_zscore(all_data)
  }
  
  # loop through each compound in dataset
  for (c in unique(all_data$trt)) {
    make_multi(c)
  } 
}




if (set_output == "total_perturb") {
  
  # normalize data by median subtracted z-score
  if (set_normalization == "plate") {
    z_norm_data <- normalize_by_plate_v3(all_data)
  } else {
    z_norm_data <- normalize_by_zscore(all_data)
  }
  
  # loop through each compound in dataset
  for (c in unique(all_data$trt)) {
    
    if (c != "Control") {
      # make plot of total perturbation (sum of parameters)
      pdf(paste(c, "_scalarPerturb.pdf", sep=""))
      scal_perturb(z_norm_data, c)
      dev.off()
    }
  }  
}



if (set_output == "rand_total_perturb") {
  
  # normalize data by median subtracted z-score
  if (set_normalization == "plate") {
    z_norm_data <- normalize_by_plate_v3(all_data)
  } else {
    z_norm_data <- normalize_by_zscore(all_data)
  }

  # loop through each compound in dataset
  for (c in unique(all_data$trt)) {
    
    if (c != "Control") {
      # make plot of total perturbation (sum of parameters)
      pdf(paste(c, "_scalarPerturbRS.pdf", sep=""))
      rand_scal_perturb(z_norm_data, c, samples=100)
      dev.off()
    }
  }  
}




if (set_output == "total_perturb_slope") {
  
  # normalize data by median subtracted z-score
  if (set_normalization == "plate") {
    z_norm_data <- normalize_by_plate_v3(all_data)
  } else {
    z_norm_data <- normalize_by_zscore(all_data)
  }
  
  # loop through each compound in dataset
  for (c in unique(all_data$trt)) {
    
    if (c != "Control") {
      # make plot of total perturbation (sum of parameters)
      pdf(paste(c, "_scalarPerturb_slope.pdf", sep=""))
      scal_perturb(z_norm_data, c, slopes=TRUE)
      dev.off()
    }
  }  
}





if (set_output == "replicates") {
  
  # normalize data by median subtracted z-score
  if (set_normalization == "plate") {
    z_norm_data <- normalize_by_plate_v3(all_data)
  } else {
    z_norm_data <- normalize_by_zscore(all_data)
  } 
  
  # loop through each compound in dataset
  for (c in unique(all_data$trt)) {
    
    if (c != "Control") {
      # make plot of total perturbation (sum of parameters) for each replicate plateID
      pdf(paste(c, "_scalarPerturb_reps.pdf", sep=""))
      scal_perturb_reps(z_norm_data, c)
      dev.off()
    }
  }  

}













##Cutouts

## Normalize by same-plate controls using median-based z-score (on the same DIV)
#alt_plate_norm <- function(all_data) {
#  
#  # Split data by plate (date as well because plateIDs get reused) and DIV
#  per_plate_split <- split(all_data, interaction(all_data[,"date"], all_data[,"Plate.SN"], all_data[,"DIV"], drop=TRUE))
#  
#  scaled_plates <- list() #initialize list
#  
#  # Loop through each plate at each DIV
#  for (i in per_plate_split) {
#    
#    # Loop through each ontogeny parameter
#    for (j in names(all_data[,8:(ncol(all_data)-1)])) {
#     
#      cntrls <- (subset(i, trt=="Control")[,j]) # Get vector of control values for that DIV on that plate
#      i[,j] <- ((i[,j]) - median(cntrls)) # Subtract all values on that DIV and plate by median of controls
#      
#    }
#    
#    scaled_plates[[length(scaled_plates)+1]] <- i
#  }
#  
#  rm(i,j)
#  scaled_plates <- do.call(rbind, scaled_plates) #Re-form one table of values
#  
#  # Now find and divide by global standard deviations - Split scaled data set by DIV
#  # We do this because all control values on the same plate may be zero or the same, leading to division by zero otherwise
#  div_split <- split(scaled_plates, scaled_plates[,"DIV"]) 
#  
#  norm_plates <- list() #initialize list
#  
#  # Loop through each DIV
#  for (i in div_split) {
#    
#    # Loop through each ontogeny parameter
#    for (j in names(all_data[,8:(ncol(all_data)-1)])) {
#     
#      cntrls <- (subset(i, trt=="Control")[,j]) # Get vector of control values for that DIV
#      i[,j] <- ((i[,j]) / sd(cntrls)) # Divide by SD of all controls on that DIV
#      
#    }  
#    norm_plates[[length(norm_plates)+1]] <- i
#  }
#  
#  rm(i,j)
#  norm_plates <- do.call(rbind, norm_plates) #Re-form one table of values
#  
#  norm_plates   
#}

# Return one value per well that sums median Z-score for each network parameter
# norm_data_sums <- cbind(norm_data[,1:7], rowSums(norm_data[,8:(ncol(norm_data)-3)])) # Fourth from last column is last parameter used because we are excluding cv.time and cv.network for now
# names(norm_data_sums)[8] <- "Zsum"




## Separate out data for compound of interest and plot all network parameter median z-scores
## With cone of uncertainty (twice the standard deviation of control values)
#plot_dose_perturbations <- function(norm_data, parameter, compound) {
  
  # Find median of replicate response on each DIV
#  data_meds <- as.data.frame(do.call(cbind, aggregate(get(parameter) ~ trt + dose + DIV, FUN=function(x) c(med=median(x), sd=sd(x), n=length(x), low=(median(x) - (sd(x)/sqrt(length(x)))), high=(median(x) + (sd(x)/sqrt(length(x))))), data=norm_data)), stringsAsFactors=FALSE)
  
  # Subset to compound of interest and control values
#  trt_meds <- subset(data_meds, trt %in% c(compound, "Control"))
  
  # Add in DIV0 placeholder data to anchor plots at zero
#  trt_meds <- rbind(data.frame(trt="Compound", dose=as.numeric(as.character(unique(trt_meds[,"dose"]))), DIV=0, med=0, sd=0, n=1, low=0, high=0), trt_meds)
#  trt_meds[,2:ncol(trt_meds)] <- lapply(trt_meds[,2:ncol(trt_meds)], as.numeric)
  
  # Isolate control standard deviations
#  trt_cntrls <- subset(trt_meds, trt=="Control")
#  trt_cntrls <- rbind(data.frame(trt="Control", dose=0, DIV=0, med=0, sd=0, n=1, low=0, high=0), trt_cntrls)
  
  # Repeat control values same number of times as total number of samples
#  trt_cntrls2 <- trt_cntrls[rep(seq_len(nrow(trt_cntrls)), nrow(trt_meds)/nrow(trt_cntrls)),]
  
  # Change low and high values to twice the standard deviation
#  trt_cntrls2[,"low"] <- trt_cntrls2[,"med"] - 2*trt_cntrls2[,"sd"]
#  trt_cntrls2[,"high"] <- trt_cntrls2[,"med"] + 2*trt_cntrls2[,"sd"]
  
  # Prepare SEM intervals
#  limits <- aes(ymin=trt_meds[,"low"], ymax=trt_meds[,"high"])
  
  # Plot with ggplot
#  ggplot(trt_meds, aes(x=DIV, y=med, group=dose, color=log(dose+1))) +
#    geom_ribbon(ymin=trt_cntrls2$low, ymax=trt_cntrls2$high, x=trt_cntrls2$DIV, fill="grey70", alpha=0.8, inherit.aes=FALSE) +  
#    geom_line(size=0.4, position=position_jitter(w=0, h=0.1)) +
#    # geom_errorbar(limits, size=0.2, width=0.25) +  # Uncomment to add SEM error bars
#    scale_color_gradient(low="blue", high="red") + 
#    coord_cartesian(ylim=c(-5,5)) + 
#    theme(axis.text=element_text(size=4), axis.title.x=element_text(size=6), axis.title.y=element_text(size=0), plot.title=element_text(size=6), legend.text=element_text(size=2), legend.position = "none") +
#    ggtitle(paste(compound, parameter, sep=" - "))
#}




