## MEA neuronal ontogeny AUC calculations, plotting, and summary statistics
## Chris Frank - February-March 2016


## Set parameters of analysis:
set_compound <- "Acetaminophen" # This has to match compound in the data set.
set_ontogeny_AUC <- "meanfiringrate_AUC" # This can be set to one of meanfiringrate_AUC, burst.per.min_AUC, mean.isis_AUC, per.spikes.in.burst_AUC, mean.dur_AUC, mean.IBIs_AUC, nAE_AUC, nABE_AUC, ns.n_AUC, ns.peak.m_auc, ns.durn.m_AUC, ns.percent.of.spikes.in.ns_AUC, ns.mean.insis_AUC, ns.durn.sd_AUC, ns.mean.spikes.in.ns_AUC, r_AUC, cv.time_AUC, cv.network_AUC, ns.n_AUC, ns.durn.m_AUC, ns.mean.spikes.in.ns_AUC, r_AUC, cv.time_AUC, cv.network_AUC
set_ontogeny_name <- "MFR" # Any name to mark y-axis and file names
set_ylimit <- 200 # Plot max y-value (percent of control)
set_color <- "blue" # Color of the plot lines
set_cutoff <- 0 # Threshold for AUC dose-response curve (percent of control); Set to zero to produce no cutoff data
set_directory <- "~/MEA_Data/PIP3_AIM2-3-DNT/TEST5/" # Set directory 
set_data_file1 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/SA2_allData_11-17-2014.csv" # Name of csv file with burst parameters per sample
set_data_file2 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/AllData_SA3_NoDNTRef.csv"
set_data_file3 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/DNTRef_all.csv"
set_mi_data_file <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/MutualInformation/update_10-31-16/all_mi.csv"
set_output <- "ec50_table" # This can be set to "replicates", "dr_plot", "all_dr_plots", "hill_plot", "all_hill_plots", "ec50_table","toxcast_table", or "all_toxcast_tables" ; if none of the above are specified, nothing will be produced 
set_ec_value <- 50 # set EC value for threshold to be reported
set_direction <- "down"  # choose to look at "up" (increasing parameter values) or "down" (decreasing parameter values)



## Read in data and clean it up
setwd(set_directory)
all_data1 <- read.delim(set_data_file1, sep=",", stringsAsFactors=FALSE) #Reading in full data set ; assuming it is a csv table
all_data2 <- read.delim(set_data_file2, sep=",", stringsAsFactors=FALSE) # Reading in second data set
all_data3 <- read.delim(set_data_file3, sep=",", stringsAsFactors=FALSE) # Reading in third data set
all_data <- rbind(all_data1, all_data2, all_data3)
rm(all_data1, all_data2, all_data3)

bis_rows <- grep("12_01_", all_data$file.name, fixed=TRUE) #index all Bicuculline-treated wells
all_data <- all_data[- bis_rows,] #Remove all bic-treated wells
all_data[is.na(all_data)] <- 0 #Replace all NAs with zeros for AUC calculations - This may be undesirable for MEA parameters that are derived from other parameters.
all_data <- subset(all_data, DIV %in% c(5,7,9,12)) # Remove sparse DIV2 data
all_data <- unique(all_data) # eliminate duplicate data



## Read in mutual information parameter and attach to full dataset
mi_data <- read.delim(set_mi_data_file, sep=",", stringsAsFactors=FALSE)
mi_data <- subset(mi_data, DIV %in% c(5,7,9,12)) # Remove sparse DIV2 data
mi_data <- unique(mi_data) # eliminate duplicate data
all_data <- merge(all_data, mi_data, all.x=TRUE) # merge two data frames on common columns (plateID, date, well, treatment, etc.)
all_data <- all_data[, c(1:23, 27, 24:26)] # rearrange data frame to accomodate additional parameter
names(all_data)[[24]] <- "mi"



## Split data frame by individual wells over time (interaction term speeds this up greatly for larger datasets)
all_data_split <- split(all_data, interaction(all_data[,"date"], all_data[,"Plate.SN"], all_data[,"well"], all_data[,"trt"], all_data[,"dose"], drop=TRUE)) # Split data into bins of single well across time (DIV)




## Function to calculate area under the curve (AUC) for each ontogeny parameter for each bin (each experiment)
calc_auc <- function(all_data_split, sqrt=FALSE) {
  require(pracma) #pracma package has trapz function that computes AUC based on trapezoidal geometry (no curve fitting)
  
  out <- lapply(1:length(all_data_split), function(i) {
    
    all_data_split[[i]] <- all_data_split[[i]][order(all_data_split[[i]][,"DIV"]),]  # Make sure order of rows follows DIV time
    
    date <- all_data_split[[i]]$date[1]
    plate.SN <- all_data_split[[i]]$Plate.SN[1]
    well <- all_data_split[[i]]$well[1]
    trt <- all_data_split[[i]]$trt[1]
    dose <- all_data_split[[i]]$dose[1]
    units <- all_data_split[[i]]$units[1]
    
    # calculating AUC for each parameter after adding div=2, param=0 data point to each AUC calculation with append() to represent no activity at first timepoint
    for (j in names(all_data_split[[i]][,8:(ncol(all_data_split[[i]])-1)])) {
      param_name <- paste(j, "_auc", sep="") # create variable name
      assign(param_name, round(trapz(append(all_data_split[[i]][,"DIV"], 2, after=0), append(all_data_split[[i]][,j], 0, after=0)),3), inherits=TRUE) # calculate auc, assign to variable name
    }
    
    # put vector of AUC values together
    c(date, as.character(plate.SN), as.character(well), as.character(trt), dose, as.character(units), meanfiringrate_auc, burst.per.min_auc, mean.isis_auc, per.spikes.in.burst_auc, mean.dur_auc, mean.IBIs_auc, nAE_auc, nABE_auc, ns.n_auc, ns.peak.m_auc, ns.durn.m_auc, ns.percent.of.spikes.in.ns_auc, ns.mean.insis_auc, ns.durn.sd_auc, ns.mean.spikes.in.ns_auc, r_auc, mi_auc, cv.time_auc, cv.network_auc)
  })
  
  sum_table <- as.data.frame(do.call(rbind, out), stringsAsFactors=FALSE) # Re-form data frame
  names(sum_table) <- c("date","plate.SN","well","treatment","dose","units","meanfiringrate_AUC","burst.per.min_AUC","mean.isis_AUC","per.spikes.in.burst_AUC","mean.dur_AUC","mean.IBIs_AUC", "nAE_AUC","nABE_AUC","ns.n_AUC","ns.peak.m_AUC","ns.durn.m_AUC","ns.percent.of.spikes.in.ns_AUC","ns.mean.insis_AUC","ns.durn.sd_AUC","ns.mean.spikes.in.ns_AUC","r_AUC","mi_AUC","cv.time_AUC","cv.network_AUC")
  sum_table[,7:ncol(sum_table)] <- lapply(sum_table[,7:ncol(sum_table)], as.numeric) # Change output columns class to numeric
  sum_table[sum_table[,"dose"]==0,"treatment"] <- "Control" # Convert all zero dose rows to say "Control" for treatment
  
  if (sqrt==TRUE){
    sum_table <- cbind(sum_table[,1:6], sqrt(sum_table[,7:25]))
  }
  
  sum_table
}




## function to normalize by plate for each ontogeny parameter
auc_summary <- function(ontogeny, summary_table, direction=set_direction) {
  
  # Split data set by plate (date as well because plateIDs get reused)
  per_plate_split <- split(summary_table, interaction(summary_table[,"date"], summary_table[,"plate.SN"], drop=TRUE))
  
  # Normalize each plate by percent of control median and combine into single data frame
  norm_plates <- list() # initialize list
  
  # Loop through each plate
  for (i in per_plate_split) {
    
    if (i[,"plate.SN"][[1]]!="MW1008-41") { # Exclude any individual plates where untreated wells are off 
    
      # Loop through each ontogeny parameter
      for (j in names(summary_table[,7:ncol(summary_table)])) {
    
        cntrls <- (subset(i, treatment=="Control")[,j]) # Get vector of control values for that DIV on that plate
        i[,j] <- (i[,j] / median(na.omit(cntrls)))*100 # Divide all values on that plate by controls median
      }
    
    } else { # For plates singled-out above, use culture median instead of plate median
      
      for (j in names(summary_table[,7:ncol(summary_table)])) {
        
        cntrls <- (subset(summary_table, treatment=="Control" & date==i[,"date"][[1]])[,j]) # identify same-date control wells
        i[,j] <- (i[,j] / median(na.omit(cntrls)))*100  # Divide by same-culture control wells median     
      }
    }
    
    norm_plates[[length(norm_plates)+1]] <- i  # Add to growing list of normalized plates
  }
  
  rm(i,j) 
  norm_plates <- do.call(rbind, norm_plates) #Re-form one table of values
  
  # if set_direction is up, take positive difference as negative
  if (direction=="up") {
    norm_plates[,7:ncol(norm_plates)] <- 100-(norm_plates[,7:ncol(norm_plates)]-100)
  }
  
  #print(norm_plates) #checkpoint for normalization
  norm_plates
}





## function to calculate CV of control wells for a given (normalized) ontogeny parameter
calc_cntrl_mad <- function(ontogeny, norm_data) {
 
  cntrl_vals <- subset(norm_data, treatment=="Control")[,ontogeny] # isolate control values for that parameter
  
  cntrl_thresh <- 100 - (3*mad(cntrl_vals)) # Calculate 3 times the median absolute deviation 
  
  # If 3 times the MAD takes us below zero, set to zero
  if (cntrl_thresh < 0) {
    cntrl_thresh <- 0
  }
  
  cntrl_thresh
}



 
## function to plot each replicate of AUC calculations for the given ontogeny parameter 
plotOverTime <- function(split_data, ontogeny, compound, ontogeny_label, ylim) {
  require(pracma)
  
  ontogeny <- substr(ontogeny, 1, nchar(ontogeny)-4) #Removing "_AUC" from ontogeny name
  
  pdf(paste(paste(ontogeny, "AUC-reps", sep="-"),".pdf",sep=""))
  par(mar=c(5.1,5.1,4.1,2.1))
  
  for (i in split_data) {
    auc=round(trapz(append(i[,"DIV"], 2, after=0), append(i[,ontogeny], 0, after=0)),3)
    plot(append(i[,"DIV"], 2, after=0), append(i[,ontogeny], 0, after=0), ylim=c(0, ylim), pch=19, xlab="DIV", ylab=ontogeny_label, cex=1.2, cex.lab=1.5, cex.axis=1.25)
    lines(append(i[,"DIV"], 2, after=0), append(i[,ontogeny], 0, after=0),lwd=2)
    text(x=4, y=0.8*ylim, labels=paste("AUC=",auc, sep=" "), cex=1.5)
    title(main=paste(i$trt[1], i$dose[1], i$units[1], i$date[1], i$Plate.SN[1], sep=" "))
  }
  
  dev.off()
}







## function to model AUC-dose response with log-logistic curves (Hill curves) and report EC50 values
plotHillCurve <- function(compound, ontogeny, ontogeny_label) {
  
  require(drc)
  
  comp_data <- subset(norm_plates, treatment==compound) # isolate data for compound of interest
  
  x <- as.numeric(comp_data[,"dose"])
  y <- comp_data[,ontogeny]
  
  # Get summary statistics for plotting mean +/- SD
  auc_stats <- as.data.frame(do.call(cbind, aggregate(get(ontogeny) ~ dose, FUN=function(x) c(mean=mean(x), sd=sd(x), sem=(sd(x)/length(x)), n=length(x), min=min(x), max=max(x)), data=comp_data)), stringsAsFactors=FALSE)
  auc_stats[,1:ncol(auc_stats)] <- lapply(auc_stats[,1:ncol(auc_stats)], as.numeric)
  b <- auc_stats[order(auc_stats[,"dose"]),] # Make sure it's sorted by ascending dose
  
  print(b)
  
  d <- c()
  for (i in 1:(nrow(b)-1)) { d[i] <- (b[i+1,"mean"] - b[i,"mean"]) } # get difference between consecutive dose mean values
  
  # Calculate 3x median absolute deviation (MAD) threshold for given parameter
  lim <- calc_cntrl_mad(ontogeny, norm_plates)
  print(paste("3 x MAD limit is ", lim, sep=""))
  
  # Alternative is to use fixed 50% reduction from control value as threshold
  # lim <- 50
  
  # Set up estimation of threshold dose without curve fitting function for when curve fitting is very poor
  est_wo_curve <- function() {
    # Find dose at which threshold is crossed for conservative EC50 estimate
    conc <- b[,"dose"]
    mn_auc <- b[,"mean"]
    cutoff_dose <- conc[min(which(mn_auc < set_ec_value))] # determines vector position of first value below cutoff and grabs corresponding dose
    prev_dose <- conc[(min(which(mn_auc < set_ec_value)))-1] # get previous dose
    cutoff_mn <- mn_auc[min(which(mn_auc < set_ec_value))] # get mean AUC that first crosses cutoff
    prev_mn <- mn_auc[(min(which(mn_auc < set_ec_value)))-1] # get previous mean AUC
    if (cutoff_dose == min(b[,"dose"])) { prev_dose <- 0 ; prev_mn <- 100 } # if lowest dose crosses cutoff, set previous dose to zero
    yval <- c(prev_dose, cutoff_dose)
    xval <- c(prev_mn, cutoff_mn)
    f <- approxfun(xval, yval)
    approx_dose <- round(f(lim), digits=2)
    approx_dose
  }  
  
  # Set up smooth spline curve visual in case log-logistic fails
  est_curve <- smooth.spline(b[,"dose"], y=b[,"mean"], spar=0.6)

  
  ## Check if the mean AUC falls 3xMAD below control mean (100%) for greatest two doses tested; If not the model may never converge and likely we cannot reliably estimate EC50 value.
  if (min(b[(nrow(b) - 1):nrow(b),"mean"]) < min(c(lim, set_ec_value))) {
    
    # Try to fit a 4-parameter log-logistic model with lower limit at 0 and upper limit at 100 (hill model)
    m1 <- try(drm(y ~ x, fct=LL.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), control=drmc(maxIt=3000, method="Nelder-Mead")), silent=TRUE)    
    
    # Generate plot
    pdf(paste(paste(compound, ontogeny_label, "EC50", sep="-"),".pdf",sep=""))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    # If model converged, plot the model
    if (!is(m1, "try-error")) {
      
      print(summary(m1))
      
      # Check if model is an OK fit of data; if it is, find values and plot fit
      # Currently changed to zero because we are accepting all model fits
      fit_test <- modelFit(m1)
      print(fit_test)
      if (fit_test[,"p value"][2] > 0) {
      
        # Calculate EC30 and EC50 values and associated standard errors
        cutoff_doses <- ED(m1, c(30,50), interval="delta")
        ec30_dose <- round(cutoff_doses[1,1], digits=2)
        ec30_se <- round(cutoff_doses[1,2], digits=2)
        ec30_low <- ec30_dose - ec30_se
        ec50_dose <- round(cutoff_doses[2,1], digits=2)
        ec50_se <- round(cutoff_doses[2,2], digits=2)
        ec50_low <- ec50_dose - ec50_se
        upper_limit <- 100
        
        # If EC30 is above maximum tested dose, ignore it
        if (ec30_dose > max(x)) {
          ec30_dose <- paste(">", max(x), sep="")
        }
        
        # If EC50 is above maximum tested dose, ignore it
        if (ec50_dose > max(x)) {
          ec50_dose <- paste(">", max(x), sep="")
        }
        
        plot(m1, lwd=2, pch="",xlab=expression(paste("Concentration (", mu,"M)")), ylab="Percent of control", ylim=c(0,200), cex=1.2, cex.lab=1.5, cex.axis=1.25, col="blue")
        #points(m1$data, col="black", pch=19) # Uncomment to plot individual replicate datapoints
        #points(b[,"dose"], b[,"mean"], pch=3, col="black", cex=2)
        points(b[,"dose"], b[,"mean"], pch=19, cex=1.2, col="black")
        #points(0, 100, pch=3, col="black", cex=2)
        #points(0, 100, pch=19, col="black", cex=1.2)
        suppressWarnings(arrows(b[,"dose"], b[,"mean"] - b[,"sd"], b[,"dose"], b[,"mean"] + b[,"sd"], length=0.10, angle=90, code=3, lwd=1.5, col="black"))
        #abline(h=.70*upper_limit, col="grey70", lwd=2, lty=2)
        abline(h=.50*upper_limit, col="grey70", lwd=2, lty=2)
        abline(h=lim, col="darkorange1", lwd=2, lty=2)
        title(main=paste(compound, ontogeny_label, sep=" - "), cex.main=1.5)
        legend("topright", paste("EC50 = ", ec50_dose, "然"), col="red", cex=1.2)
        
        # If standard error of EC estimate is undetermined,
        # Or if minimum of confidence interval is less than lowest tested dose and slope is close to 0, 
        # Fit and EC30/50 is poorly defined ; add warning message and estimate threshold dose
        if ( is.na(ec50_se) || (ec50_low < min(x[x > 0]) && m1$fit$par[1] < 0.2) || (ec30_low < min(x[x > 0]) && m1$fit$par[1] < 0.2)) {
          
          legend("topleft", "Fit is poor!", cex=1.2) 
          
          EC50 <- est_wo_curve()
          
          # Plot with conservative EC50 estimate
          plot(b[,"dose"], b[,"mean"], pch=3, cex=0, xlab=expression(paste("Concentration (", mu,"M)")), ylab="Percent of control", ylim=c(0,200), log="x", cex.lab=1.5, cex.axis=1.25, col="black", las=1)
          lines(est_curve, col="blue", lwd=2)
          suppressWarnings(arrows(b[,"dose"], b[,"mean"] - b[,"sd"], b[,"dose"], b[,"mean"] + b[,"sd"], length=0.10, angle=90, code=3, lwd=1.5, col="black"))
          points(b[,"dose"], b[,"mean"], pch=19, cex=1.2, col="black")
          #points(x, y, col="black", pch=19) # Uncomment to plot individual replicate datapoints
          #abline(h=70, col="grey70", lwd=2, lty=2)
          abline(h=50, col="grey70", lwd=2, lty=2)
          abline(h=lim, col="darkorange1", lwd=2, lty=2)
          title(main=paste(compound, ontogeny_label, sep=" - "), cex.main=1.5)
          legend("topright", paste("EC50 = ", EC50, "然"), col="red", cex=1.2)
        }
          
        
      # If model is poor fit by ANOVA test, just plot data points and estimate threshold-crossing dose  
      } else {
        
        print("Fit is poor! Manually estimating threshold dose.")
        
        EC50 <- est_wo_curve()
        
        # Plot with conservative EC50 estimate
        plot(b[,"dose"], b[,"mean"], pch=3, cex=0, xlab=expression(paste("Concentration (", mu,"M)")), ylab="Percent of control", ylim=c(0,200), log="x", cex.lab=1.5, cex.axis=1.25, col="black", las=1)
        lines(est_curve, col="blue", lwd=2)
        suppressWarnings(arrows(b[,"dose"], b[,"mean"] - b[,"sd"], b[,"dose"], b[,"mean"] + b[,"sd"], length=0.10, angle=90, code=3, lwd=1.5, col="black"))
        points(b[,"dose"], b[,"mean"], pch=19, cex=1.2, col="black")
        #points(x, y, col="black", pch=19) # Uncomment to plot individual replicate datapoints
        #abline(h=70, col="grey70", lwd=2, lty=2)
        abline(h=50, col="grey70", lwd=2, lty=2)
        abline(h=lim, col="darkorange1", lwd=2, lty=2)
        title(main=paste(compound, ontogeny_label, sep=" - "), cex.main=1.5)
        legend("topright", paste("EC50 = ", EC50, "然"), col="red", cex=1.2) 
      }   
        
        
    # If model failed to converge, plot data points and estimate threshold-crossing dose  
    } else {
      
      EC50 <- est_wo_curve()
      
      # Plot with conservative EC50 estimate
      plot(b[,"dose"], b[,"mean"], pch=3, cex=0, xlab=expression(paste("Concentration (", mu,"M)")), ylab="Percent of control", ylim=c(0,200), log="x", cex.lab=1.5, cex.axis=1.25, col="black", las=1)
      lines(est_curve, col="blue", lwd=2)
      suppressWarnings(arrows(b[,"dose"], b[,"mean"] - b[,"sd"], b[,"dose"], b[,"mean"] + b[,"sd"], length=0.10, angle=90, code=3, lwd=1.5, col="black"))
      points(b[,"dose"], b[,"mean"], pch=19, cex=1.2, col="black")
      #points(x, y, col="black", pch=19) # Uncomment to plot individual replicate datapoints
      #abline(h=70, col="grey70", lwd=2, lty=2)
      abline(h=50, col="grey70", lwd=2, lty=2)
      abline(h=lim, col="darkorange1", lwd=2, lty=2)
      title(main=paste(compound, ontogeny_label, sep=" - "), cex.main=1.5)
      legend("topright", paste("EC50 = ", EC50, "然"), col="red", cex=1.2)
    }  
      
    dev.off()
    
    # If model converged, print EC50 value
    if (!is(m1, "try-error") && exists("ec50_dose")) {
      print(paste(compound, " EC50 value is ", ec50_dose, " microMolar", sep=""))
    # If not, set to greater than max dose
    } else { 
      print(paste(compound, " EC50 value is ", EC50, " microMolar", sep=""))
    }
    print("______________________________________________________________")
    
  
  # If hill model was not attempted, plot data points and set EC50 to above max dose  
  } else {
      pdf(paste(paste(compound, ontogeny_label, "EC50", sep="-"),".pdf",sep=""))
    
      par(mar=c(5.1,5.1,4.1,2.1))
    
      EC50 <- max(b[,"dose"])
      plot(b[,"dose"], b[,"mean"], pch=3, cex=0, xlab=expression(paste("Concentration (", mu,"M)")), ylab="Percent of control", ylim=c(0,200), log="x", cex.lab=1.5, cex.axis=1.25, col="black", las=1)
      lines(est_curve, col="blue", lwd=2)
      suppressWarnings(arrows(b[,"dose"], b[,"mean"] - b[,"sd"], b[,"dose"], b[,"mean"] + b[,"sd"], length=0.10, angle=90, code=3, lwd=1.5, col="black"))
      points(b[,"dose"], b[,"mean"], pch=19, cex=1.2, col="black")
      #points(x, y, col="black", pch=19) # Uncomment to plot individual replicate datapoints
      abline(h=lim, col="darkorange1", lwd=2, lty=2)
      title(main=paste(compound, ontogeny_label, sep=" - "), cex.main=1.5)
      legend("topright", paste("EC50 > ", EC50, "然"), col="red", cex=1.2)
    
      dev.off()
    
      # Print EC50 value not found
      print(paste(compound, " EC50 value is greater than max dose of ", EC50, " microMolar", sep=""))
      print("_____________________________________________________________")
  }
}







## function to make summary table of all AUC EC50 values for all ontogeny parameters
generate_EC50_table <- function(compound, ontogeny, ec_level=50) {
  
  require(drc)
  
  comp_data <- subset(norm_plates, treatment==compound) # isolate data for compound of interest
  
  x <- as.numeric(comp_data[,"dose"])
  y <- comp_data[,ontogeny]
  
  # Get summary statistics for plotting mean +/- SD
  auc_stats <- as.data.frame(do.call(cbind, aggregate(get(ontogeny) ~ dose, FUN=function(x) c(mean=mean(x), sd=sd(x), sem=(sd(x)/length(x)), n=length(x), min=min(x), max=max(x)), data=comp_data)), stringsAsFactors=FALSE)
  auc_stats[,1:ncol(auc_stats)] <- lapply(auc_stats[,1:ncol(auc_stats)], as.numeric)
  b <- auc_stats[order(auc_stats[,"dose"]),] # Make sure it's sorted by ascending dose
  
  print(b)
  
  d <- c()
  for (i in 1:(nrow(b)-1)) { d[i] <- (b[i+1,"mean"] - b[i,"mean"]) } # get difference between consecutive dose mean values
  
  # Calculate 3x median absolute deviation (MAD) threshold for given parameter
  lim <- calc_cntrl_mad(ontogeny, norm_plates)
  print(paste("3 x MAD limit is ", lim, sep=""))
  
  # Alternatively, we could use fixed 50% reduction from control value as threshold beyond which we fit a curve
  # lim <- 50
  
  # Set up estimation of threshold dose without curve fitting function
  est_wo_curve <- function() {
    # Find dose at which threshold is crossed for conservative EC50 estimate
    conc <- b[,"dose"]
    mn_auc <- b[,"mean"]
    cutoff_dose <- conc[min(which(mn_auc < set_ec_value))] # determines vector position of first value below cutoff and grabs corresponding dose
    prev_dose <- conc[(min(which(mn_auc < set_ec_value)))-1] # get previous dose
    cutoff_mn <- mn_auc[min(which(mn_auc < set_ec_value))] # get mean AUC that first crosses cutoff
    prev_mn <- mn_auc[(min(which(mn_auc < set_ec_value)))-1] # get previous mean AUC
    if (cutoff_dose == min(b[,"dose"])) { prev_dose <- 0 ; prev_mn <- 100 } # if lowest dose crosses cutoff, set previous dose to zero
    yval <- c(prev_dose, cutoff_dose)
    xval <- c(prev_mn, cutoff_mn)
    f <- approxfun(xval, yval) # create linear function between two points
    approx_dose <- round(f(lim), digits=2) # estimate threshold dose along the line
    approx_dose_se <- 0.5*(cutoff_dose - prev_dose) # estimate stderr by 1/2 times difference of flanking doses
    approx_dose_low <- approx_dose - approx_dose_se
    approx_dose_hi <- approx_dose + approx_dose_se
    return(c(approx_dose, approx_dose_se, approx_dose_low, approx_dose_hi))
  } 
  
  ## Check if the mean AUC falls >3*MAD lower for greatest two doses tested, if so attempt hill model fit
  if (min(b[(nrow(b) - 1):nrow(b),"mean"]) < min(c(lim, set_ec_value))) {
    
    # Try to fit a 4-parameter log-logistic model with lower limit at 0 and upper limit at 100 (hill model)
    m1 <- try(drm(y ~ x, fct=LL.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), control=drmc(maxIt=3000, method="Nelder-Mead")), silent=TRUE)
    
    # If hill model converged, report EC value and standard error
    if (!is(m1, "try-error")) {
      
      # Check if model is an OK fit of data
      # Currently changed to zero because we are accepting all model fits
      fit_test <- modelFit(m1)
      if (fit_test[,"p value"][2] > 0) {
        
        # Calculate EC50 values and associated standard errors
        cutoff_dose <- ED(m1, ec_level, interval="delta")
        ec_dose <- round(cutoff_dose[1,1], digits=2)
        ec_se <- round(cutoff_dose[1,2], digits=2)
        ci_low <- round(cutoff_dose[1,1] - cutoff_dose[1,2], digits=2)
        ci_high <- round(cutoff_dose[1,1] + cutoff_dose[1,2], digits=2)
        
        #rmse <- round(sqrt(mean((fitted(m1) - y)^2)), digits=3)  # This gives the root mean square error (RMSE)
        #aic <- round(AIC(m1), digits=3) # This gives the Akaike information criterion where a lower score is better 
        
        # If standard error of EC estimate is undetermined (usually due to very poor fit, lack of numerical accuracy),
        # Or if minimum of confidence interval is less than lowest tested dose and slope is close to 0, EC50 is poorly defined by curve; determining value by straight lines
        if (is.na(ec_se) || (ci_low < min(x[x > 0]) && m1$fit$par[1] < 0.2)) {
          ec_est <- est_wo_curve()
          ec_dose <- ec_est[1]
          ec_se <- ec_est[2]
          ci_low <- ec_est[3]
          ci_high <- ec_est[4]
        }
        
      # When model fit is poor by ANOVA test, ignore model and determine by straight line estimate 
      } else {
          ec_est <- est_wo_curve()
          ec_dose <- ec_est[1]
          ec_se <- ec_est[2]
          ci_low <- ec_est[3]
          ci_high <- ec_est[4]
      }
      
    # If hill model failed to converge, determine by straight line estimate  
    } else {
      ec_est <- est_wo_curve()
      ec_dose <- ec_est[1]
      ec_se <- ec_est[2]
      ci_low <- ec_est[3]
      ci_high <- ec_est[4]
    }

  # If declined to try hill model at all, set to NA        
  } else {
    ec_dose <- NA
    ec_se <- NA
    ci_low <- NA
    ci_high <- NA
  }
  
  
  # If model EC50 estimate is still greater than tested dose range, set to NA
  if (!is.na(ec_dose)){
    if (ec_dose > max(x)) {
      ec_dose <- NA
      ec_se <- NA
      ci_low <- NA
      ci_high <- NA
    }
  }
  
  
  ec50_out <- data.frame("compound"=compound, "ontogeny"=ontogeny, "ec"=ec_dose, "stderr"=ec_se, "ci_low"=ci_low, "ci_high"=ci_high)
  data.frame(ec50_out)
  
} 






## Generate potential outputs:

sum_table <- calc_auc(all_data_split) # Compute AUC values (needed for everything below)
norm_plates <- auc_summary(ontogeny, sum_table) # normalize by same-plate control medians



# Making dose-response curve based on AUC values
if (set_output=="dr_plot") {  
  auc_vs_conc_plot(sum_table, set_ontogeny_AUC, set_compound, set_ontogeny_name, set_ylimit, set_color, set_cutoff)
}


# Making AUC dose-response curves for all compounds in data set
if (set_output=="all_dr_plots") { 
  for (i in unique(all_data[,"trt"])) {
    auc_vs_conc_plot(sum_table, set_ontogeny_AUC, i, set_ontogeny_name, set_ylimit, set_color, set_cutoff)
  }
}


# Making output table of AUC values per well
if (set_output=="toxcast_table") {
  require(car)
  
  sum_table[,2] <- paste(sum_table[,2], sum_table[,1], sep="_") #Make plateID unique by combining date and plateSN
  sum_table <- sum_table[,2:ncol(sum_table)]
  
  auc_table <- cbind(sum_table[,1], substr(sum_table[,2], 1, 1), substr(sum_table[,2], 2, 2), sum_table[,3:4], sum_table[,set_ontogeny_AUC]) #Select ontogeny parameter for AUC numbers ; Separate out well info into two columns
  names(auc_table) <- c("experimentID","row","column","treatment","dose",set_ontogeny_AUC)
  
  auc_table <- apply(auc_table, 2, function(x) {x <- recode(x, "c('A','a')=1 ; c('B','b')=2 ; c('C','c')=3 ; c('D','d')=4 ; c('E','e')=5 ; c('F','f')=6") ; x}) #Replace well letters with numbers
  
  write.table(auc_table, file=paste(set_ontogeny_AUC,"_tc.csv",sep=""), quote=FALSE, sep=",", row.name=FALSE) #Create csv text file of plateID, well, compound, dose, AUC
}


# Making output table of AUC values per well for all ontogeny paramters
if (set_output=="all_toxcast_tables") {
  require(car)
  
  sum_table[,2] <- paste(sum_table[,2], sum_table[,1], sep="_") #Make plateID unique by combining date and plateSN
  sum_table <- sum_table[,2:ncol(sum_table)]
  
  for (i in c("meanfiringrate_AUC","burst.per.min_AUC","mean.isis_AUC","per.spikes.in.burst_AUC","mean.dur_AUC","mean.IBIs_AUC", "nAE_AUC","nABE_AUC","ns.n_AUC","ns.peak.m_auc","ns.durn.m_AUC","ns.percent.of.spikes.in.ns_AUC","ns.mean.insis_AUC","ns.durn.sd_AUC","ns.mean.spikes.in.ns_AUC","r_AUC", "cv.time_AUC", "cv.network_AUC")){

    auc_table <- cbind(sum_table[,1], substr(sum_table[,2], 1, 1), substr(sum_table[,2], 2, 2), sum_table[,3:4], sum_table[,i]) #Select ontogeny parameter for AUC numbers ; Separate out well info into two columns
    names(auc_table) <- c("experimentID","row","column","treatment","dose",i)
    
    auc_table <- apply(auc_table, 2, function(x) {x <- recode(x, "c('A','a')=1 ; c('B','b')=2 ; c('C','c')=3 ; c('D','d')=4 ; c('E','e')=5 ; c('F','f')=6") ; x}) #Replace well letters with numbers
    
    write.table(auc_table, file=paste(i,"_tc.csv",sep=""), quote=FALSE, sep=",", row.name=FALSE) #Create csv text file of plateID, well, compound, dose, AUC    
  } 
}



# Making plots of raw data per replicate upon which AUC calculations are made
if (set_output=="replicates") {
  plotOverTime(all_data_split, set_ontogeny_AUC, set_compound, set_ontogeny_name, 16)
}


# Making AUC dose-response plot with Hill curve and EC50 value
if (set_output=="hill_plot") {
  plotHillCurve(set_compound, set_ontogeny_AUC, set_ontogeny_name)
}


# Making AUC dose-response plots with Hill curve and EC50 values for all compounds in data set
if (set_output=="all_hill_plots") { 
  for (i in unique(all_data[,"trt"])) {
    plotHillCurve(i, set_ontogeny_AUC, set_ontogeny_name)
  }
}


# Making table with all EC50 values for each ontogeny parameter
if (set_output=="ec50_table") {
  ont_names <- names(sum_table[,7:ncol(sum_table)])
  results=list()
  n <- 1
  for (i in unique(all_data[,"trt"])) {
    for (j in ont_names) {
      results[[n]] <- generate_EC50_table(i, j, ec_level=set_ec_value)
      n <- n + 1
    }  
  }
  results2 <- as.data.frame(do.call(rbind, results), stringsAsFactors=FALSE)
  write.table(results2, paste("ec", set_ec_value, "_allOntogeny.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
}







# Cut-outs:

#rmse <- sqrt(mean((fitted(m1) - y)^2))  # This gives the root mean square error (RMSE)
#aic <- AIC(m1) # This gives the Akaike information criterion where a lower score is better 


