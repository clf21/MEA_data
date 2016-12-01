## Performing intitial steps of tipping point analysis on normalized data
## Chris Frank - February-May 2016


set_directory <- "~/MEA_Data/PIP3_AIM2-3-DNT/" # Set directory 
set_data_file1 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/SA2_allData_11-17-2014.csv" # Name of csv file with burst parameters per sample
set_data_file2 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/AllData_SA3_NoDNTRef.csv"
set_data_file3 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/DNTRef_all.csv"
set_data_file4 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/AB_CytotoxData.csv"
set_data_file5 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/LDH_CytotoxData.csv"
set_output <- "recovery"  # Can be set to "zscore_plots" (plot individual parameters), "total_perturb" (look at scalar perturb and calculate velocity and critical concentration), or "rand_total_perturb"
set_normalization <- "plate" # Can be set to "global" or "plate"
set_parameters <- c("MFR","#AE","BR","#ABE","%SiB","#NS","%SiNS","R") # A vector of parameters used for the total perturbation calculations ; can include any of names(all_data[8:(ncol(all_data)-3)])              


set.seed(1)
# Load packages
library(ggplot2)


## Source and clean up dataset
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



## Source cell viability data
viability_data1 <- na.omit(read.delim(set_data_file4, sep=",", header=TRUE, stringsAsFactors=FALSE))
viability_data2 <- na.omit(read.delim(set_data_file5, sep=",", header=TRUE, stringsAsFactors=FALSE))
all_via <- merge(viability_data1, viability_data2)
all_via[,"median"] <- apply(all_via[,3:8], 1, median)
all_via[,"median"] <- all_via[,"median"]*100
names(all_via)[1] <- "trt" # changing column heading to match network data




## Apply normalization to the plate level by dividing raw values by plate median of controls only, then scaling with z-score across plates
## The change from v2 above is now computing z-score based on control well distribution only (i.e. the z-score reflects number of SD away from control median)
normalize_by_plate <- function(all_data) {
  
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









## Calculate scalar measure of total network parameter perturbation for EACH REPLICATE
scal_perturb_reps <- function(norm_data, compound) {
  
  # Subset normalized data to compound of interest
  comp_norm <- subset(norm_data, trt==compound)
  
  # Split compound data by date, plate ID, and well to separate replicates
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
    geom_smooth(aes(y=sum), size=1.4, se=TRUE, span=0.8, col="blue", fill="blue") +
    geom_smooth(aes(y=div9to12*5), size=1.4, se=TRUE, span=0.8, col="forestgreen", fill="forestgreen") +
    geom_point(aes(y=sum), col="blue", position=position_jitter(w=0.2, h=0.25)) +
    geom_point(aes(y=div9to12*5), col="forestgreen", position=position_jitter(w=0.2, h=0.25)) +
    coord_cartesian(ylim=c(-20,20)) +
    ggtitle(paste(compound, " DIV9-12", sep="")) +
    theme(axis.text=element_text(size=14), axis.title.x=element_text(size=16), axis.title.y=element_text(size=0), plot.title=element_text(size=18), legend.text=element_text(size=14)) 
  
  
  # Plot DIV7 to DIV12 change
  vis3 <- ggplot(trt_derivs, aes(x=log(dose))) +
    geom_hline(yintercept=final_var, linetype=2) +
    geom_hline(yintercept=0, linetype=2) +
    geom_smooth(aes(y=sum), size=1.4, se=TRUE, span=0.8, col="blue", fill="blue") +
    geom_smooth(aes(y=div7to12*5), size=1.4, se=TRUE, span=0.8, col="forestgreen", fill="forestgreen") +
    geom_point(aes(y=sum), col="blue", position=position_jitter(w=0.2, h=0.25)) +
    geom_point(aes(y=div7to12*5), col="forestgreen", position=position_jitter(w=0.2, h=0.25)) +
    coord_cartesian(ylim=c(-20,20)) +
    ggtitle(paste(compound, " DIV7-12", sep="")) +
    theme(axis.text=element_text(size=14), axis.title.x=element_text(size=16), axis.title.y=element_text(size=0), plot.title=element_text(size=18), legend.text=element_text(size=14)) 

  
  # Plot DIV5 to DIV12 change
  vis4 <- ggplot(trt_derivs, aes(x=log(dose))) +
    geom_hline(yintercept=final_var, linetype=2) +
    geom_hline(yintercept=0, linetype=2) +
    geom_smooth(aes(y=sum), size=1.4, se=TRUE, span=0.8, col="blue", fill="blue") +
    geom_smooth(aes(y=div5to12*5), size=1.4, se=TRUE, span=0.8, col="forestgreen", fill="forestgreen") +
    geom_point(aes(y=sum), col="blue", position=position_jitter(w=0.2, h=0.25)) +
    geom_point(aes(y=div5to12*5), col="forestgreen", position=position_jitter(w=0.2, h=0.25)) +
    coord_cartesian(ylim=c(-20,20)) +
    ggtitle(paste(compound, " DIV5-12", sep="")) +
    theme(axis.text=element_text(size=14), axis.title.x=element_text(size=16), axis.title.y=element_text(size=0), plot.title=element_text(size=18), legend.text=element_text(size=14)) 
            
  
  # Find derivative of slope function
  add_deriv <- function(dv_plot, timeframe) {
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
    dv_plot <- dv_plot + geom_smooth(aes(x=x, y=y), data=deriv, inherit.aes=FALSE, size=1.4, se=FALSE, span=0.8, col="red", fill="red") +
      #geom_ribbon(aes(x=x, ymax=ymax, ymin=ymin), data=deriv, inherit.aes=FALSE, alpha=0.5, fill="red")
      geom_smooth(aes(x=x, y=ymin), data=deriv, inherit.aes=FALSE, size=0.8, se=FALSE, span=0.8, col="red", fill="red") +
      geom_smooth(aes(x=x, y=ymax), data=deriv, inherit.aes=FALSE, size=0.8, se=FALSE, span=0.8, col="red", fill="red")
    
    # Get derivative loess fitted values to remove noise in derivative function
    # The 7th function in the vis2-4 ggplot calls is now the derivative (slope change over concentration change)
    x_vals <- ggplot_build(dv_plot)$data[[7]]$x
    y_vals <- ggplot_build(dv_plot)$data[[7]]$y
    y_mins <- ggplot_build(dv_plot)$data[[8]]$y
    y_maxs <- ggplot_build(dv_plot)$data[[9]]$y
    
    # Replace derivative with smoothened values
    deriv <- data.frame(x=x_vals, y=y_vals, ymin=y_mins, ymax=y_maxs)
    
    # Get cytotoxicity data for the compound
    via_meds <- aggregate(median ~ dose, mean, data=subset(all_via, trt==compound))
    # Add a vertical line to the plot marking >50% cytotoxicity
    toxlvl <- min(c(via_meds[via_meds[,"median"] < 50,"dose"])) # find all doses with median cytotox value less than 50% and use min. dose
    dv_plot <- dv_plot + geom_vline(xintercept=log(toxlvl), col="goldenrod1")
    toxlvl <- min(c(via_meds[via_meds[,"median"] < 50,"dose"], 1000)) # find all doses with median cytotox value less than 50% and use min. dose ; if no dose exceeds 50%, set to arbitrarily high 1000
    
    # If low end of confidence interval for fitted scalar perturb at highest dose exceeds 3*MAD of controls,
    # And at least 5 y-values of the velocity curve fit are negative (recovering trajectory)
    # And at least one y-value of the derivative is distinctly positive and one y-value is negative,
    # Find first dose (x-value) where derivative crosses y=0 line
    if (max(ggplot_build(dv_plot)$data[[3]]$ymin[75:80]) > final_var && sort(ggplot_build(dv_plot)$data[[4]]$ymin[1:40])[5] < 0 && max(deriv[,"y"]) > 0.001 && min(deriv[1:40,"y"]) < -0.001) {
      
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
        dv_plot <- dv_plot + geom_vline(xintercept=log(root1), linetype=2) # Add tipping point estimate to ggplot object
        
        # If root value is at lower concentration than 50% cytotox dose, report it
        if (root1 < toxlvl) {
          cat(paste(compound, root1, sep=','))
        } else { 
          cat(paste(compound, "NA", "NA", "NA", sep=','), "\n")
        }
         
      # Otherwise, the first dose produced a positive slope so set the root to lowest dose  
      #} else {
      #  root1 <- min(trt_derivs[,"dose"])
      #  cat(paste(compound, root1, sep=','), "\n")
      
      # Otherwise, the first dose produced a positive slope, so set the tipping point to undefined
      } else {
        root1 <- NA
        cat(paste(compound, root1, "NA", "NA", sep=','), "\n")
      }
    
    # Else, we could not define a tipping point concentration
    } else {
      root1 <- NA
      cat(paste(compound, root1, "NA", "NA", sep=','), "\n")
    }
    
    
    # special fix for Dieldrin, which is lacking DIV5 data
    if (dv_plot$labels$title == "Dieldrin DIV5-12") {
      root1 <- NA
    }
      
    
    if (!is.na(root1)) {
      ## Estimating critical concentration error by fitting loess models to randomly sampled data subsets
      trt_deriv_split <- split(trt_derivs, trt_derivs[,"dose"]) # split up slope values by dose tested
      
      samps <- c() # intialize vector of roots
      while(length(na.omit(samps)) < 1000) { #loop through sampling 1000 times
        
        smpld_trt_deriv <- list() # intialize sample list
        n <- 1 # initialize counter
        for (i in trt_deriv_split){ #loop through each concentration
          smpld_trt_deriv[[n]] <- i[sample(nrow(i), nrow(i)-1),] # sample one less than number of rows (replicates) from each concentration tested
          n <- n+1
        }
    
        smpld_trt_deriv <- do.call(rbind, smpld_trt_deriv) # form sample data frame
        
        timeframe_loess <- loess(get(timeframe) ~ log(dose), smpld_trt_deriv, span=0.8) # fit localized regression
        pred_range <- seq(from=min(log(smpld_trt_deriv[,"dose"])), to=max(log(smpld_trt_deriv[,"dose"])), length.out=80) # use 80 x-values to predict y-values for the model (matching geom_smooth)
        vel_curve <- 5*predict(timeframe_loess, pred_range) # output model y value predictions
        
        xy_prime <- diff(vel_curve)/diff(pred_range) # approximate derivative by taking change in velocity curve over change in concentration
        deriv <- data.frame(x=pred_range[2:length(pred_range)], y=xy_prime)
        deriv_loess <- loess(y ~ x, deriv, span=0.8) # model derivative with loess smoothing
        deriv_curve <- predict(deriv_loess, pred_range[2:length(pred_range)]) # get smoothened derivative
        deriv <- data.frame(x=pred_range[2:length(pred_range)], y=deriv_curve) # replace df with fitted values
        
        # make sure the derivative fit crosses y=0
        if (max(deriv[,"y"]) > 0.001 && min(deriv[,"y"]) < -0.001) {
          
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
            samp_root <- round(exp(deriv[last_neg+1,"x"]), 3) # Get dose of next (positive) value
          } else {
            samp_root <- NA
          } 
          
          # else, the derivative fit did not cross y=0, set to NA  
        } else { 
          samp_root <- NA
        }
        
        samps <- append(samps, samp_root) # add sampled data root to growing vector 
        }
      
      #print(sort(na.omit(samps)))
      #print(quantile(na.omit(samps), probs=c(0.025, 0.975)))
      ci_95 <- quantile(na.omit(samps), probs=c(0.025, 0.975)) # Get 95% CI from 1000 sampled roots
      ci_low <- round(ci_95[[1]], 3)
      ci_hi <- round(ci_95[[2]], 3)
      cat(",", paste(ci_low, ci_hi, sep=","), "\n", sep="") # Report 95% CI limits
    }  
  
      # Print plot      
      print(dv_plot)
    }
  
  # write results to .csv file
  sink(file="tipping_point_list.csv", append=TRUE, split=TRUE)
  add_deriv(vis2, "div9to12")
  add_deriv(vis3, "div7to12")
  add_deriv(vis4, "div5to12")
  sink()
}
  
  
  


  
  


## Calculate scalar measure of total network parameter perturbation based on iterative random sampling of network parameters
rand_scal_perturb <- function(norm_data, compound, samples=100) {
  
  # Subset to compound of interest values
  comp_norm_data <- subset(norm_data, trt==compound)
  
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
    quant_vals <- quantile(tmt[,"sum"], probs=c(0.025,0.5,0.975)) # calculate CI
    ats_sum <- rbind(ats_sum, data.frame(trt=tmt[,"trt"][1], dose=tmt[,"dose"][1], DIV=tmt[,"DIV"][1], sum=quant_vals, n=c(1,2,3)))
  }
  
  ats_low <- subset(ats_sum, n==1) # X% of data below this value
  ats_med <- subset(ats_sum, n==2) # median value
  ats_high <- subset(ats_sum, n==3) # X% of data above this value
  ats_ci <- cbind(ats_med, ats_low[,"sum"], ats_high[,"sum"]) # Combining low and high for geom_ribbon call below
  names(ats_ci) <- c("trt","dose","DIV","sum","n","low","high")
  
  # Find cytotox viability data and get mean value (of medians) for each dose
  #via_meds <- aggregate(median ~ dose, mean, data=subset(all_via, trt==compound))
  #via_meds[,"DIV"] <- 12.5
  #via_meds[,"sum"] <- (100 - via_meds[,"median"])/(100/20) # Change % alive to % dead ; Divide by value to rescale to plot y-axis
  #via_meds$sum[via_meds$sum < 0] <- 0 # Set lower limit of 0%
  #via_meds <- data.frame(trt=compound, dose=via_meds[,"dose"], DIV=via_meds[,"DIV"], sum=via_meds[,"sum"], n=1)
  
  # Plot with ggplot2
  vis <- ggplot(ats_med, aes(x=DIV, y=sum, group=interaction(dose,n), color=log(dose))) +
    geom_ribbon(ymin=ats_ci$low, ymax=ats_ci$high, x=ats_ci$DIV, fill="grey70", alpha=0.3) +
    geom_line(size=1) +
    geom_line(data=ats_low, linetype="dotted", size=1, alpha=0.5) +
    geom_line(data=ats_high, linetype="dotted", size=1, alpha=0.5) +
    geom_point() +
    scale_color_gradient(low="blue", high="red") +
    coord_cartesian(ylim=c(0,20)) +
    scale_x_continuous(breaks = c(0,2,5,7,9,12)) +
    ggtitle(compound) +
    theme(axis.text=element_text(size=14), axis.title.x=element_text(size=16), axis.title.y=element_text(size=0), plot.title=element_text(size=18), legend.text=element_text(size=14)) 
  
  #vis <- vis + geom_point(data=via_meds, position=position_jitter(w=0.25, h=0)) # Adds viability points
  
  # Print plot
  print(vis)
}







## Calculate scalar measure of total network parameter perturbation; find first dose that causes departure from normal range and recovery.
## These trajectories might be indicative of system adaptation that makes the tissue susceptible to further challenge later
scal_perturb_recovery <- function(norm_data, compound) {
  
  # Subset normalized data to compound of interest
  comp_norm <- subset(norm_data, trt==compound)
  
  # Split compound data by date, plate ID, and well to separate replicates
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
  
  # Get median of each DIV and dose 
  trt_meds <- aggregate(sum ~ DIV + dose, data=trt_sums, FUN=median)
  
  
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
  
  # Add min and max values
  cntrl_var[,"min"] <- 0
  cntrl_var[,"max"] <- cntrl_var[,"mad"] + cntrl_var[,"med"] # adding the 3*MAD to the median perturbation of controls

  
  # Split up the different doses
  trt_meds_split <- split(trt_meds, trt_meds[,"dose"])
  
  # Loop through each dose ; keep every dose that shows recovery trend
  rcvr_dose <- c()
  for (i in trt_meds_split) {
    
    i[,"madlim"] <- cntrl_var[,"max"] # add 3xMAD cutoffs to dataframe
    i[,"diff"] <- i[,"sum"] - i[,"madlim"] # take difference for each timepoint
    
    # If one of DIV5, 7, or 9 medians is greater than 3*MAD and DIV12 shows recovery (i.e. has a lower z-score)
    max_index <- which.max(i[1:3,"diff"]) # get index of max departure from normal among DIV5,7,9
    if (i[max_index,"diff"] > 0 && i[4,"sum"] < i[max_index,"sum"] && i[4,"diff"] < i[max_index,"diff"]) {
      #print(i)
      rcvr_dose <- append(rcvr_dose, i[,"dose"][1]) # Add dose that shows recovery to vector
    }
  }

  if (length(rcvr_dose) < 1) {rcvr_dose <- NA} # If no dose showed recovery, set to NA
  cat(paste(compound, rcvr_dose, "\n")) # Print doses that produced recovery trends
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











## Normalize data by z-score of fold-change from same-plate controls
z_norm_data <- normalize_by_plate(all_data)



## Generate outputs

if (set_output == "zscore_plots") {
  
  # loop through each compound in dataset
  for (c in unique(all_data$trt)) {
    make_multi(c)
  } 
}



if (set_output == "rand_total_perturb") {

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




if (set_output == "total_perturb") {
  
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




if (set_output == "recovery") {
  
  sink("recovery_doses.txt")
  # loop through each compound in dataset
  for (c in unique(all_data$trt)) {
    
    if (c != "Control") {
      scal_perturb_recovery(z_norm_data, c)
    }
  }
  sink()
}




