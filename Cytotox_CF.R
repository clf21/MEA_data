# MEA neuronal ontogeny cytotoxicity EC50 values
# Chris Frank - March 2016


# Set parameters of analysis:
set_directory <- "~/MEA_Data/PIP3_AIM2-3-DNT/TEST4/Cytotox/" # Set directory
set_data_file <- "LDH_CytotoxData.csv" # Name of csv file with cytotox data (compound name, dose, then three columns of AB cytotox values)
set_output <- "LDHcytotox_EC50_summary.txt" # Name for EC50 summary table file
set_assay <- "LDHcytotox" # Name of assay to include in .pdf plots



# Load required packages
require(drc)
require(reshape2)


# Read in data and clean it up
setwd(set_directory)
all_data <- na.omit(read.delim(set_data_file, sep=",", header=TRUE, stringsAsFactors=FALSE)) # read in data, remove empty rows (marked with NA)

# Calculate median of triplicate values
# We are using median because there are occasional outliers
#all_data[,"median"] <- apply(all_data[,3:5], 1, median)
#all_data[,"median"] <- all_data[,"median"]*100

# Melt data frame
all_data_mlt <- melt(all_data, id.vars=c("compound","dose"))

# Multiply all cytotox values by 100 to get percent of control viability
all_data_mlt[,"value"] <- all_data_mlt[,"value"]*100 

# Split data frame by compound
all_data_split <- split(all_data_mlt, all_data[,"compound"])



## Attempt hill curve fitting for each compound, if it fails, plot data points without curve
# Report EC50 estimate; if no curve fit, report as above max dose
sink(file=set_output, split=TRUE)
for (i in all_data_split) {
  
  # Set compound name
  compound <- i[,"compound"][1]
  
  x <- i[,"dose"]
  y <- i[,"value"]
  
  cmp_means <- aggregate(value ~ dose, data=i, FUN=mean) # get mean percent for each dose
  
  if (min(cmp_means[,"value"]) < 50) {
    
    m1 <- try(drm(y ~ x, fct=LL.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), control=drmc(maxIt=3000, method="Nelder-Mead")), silent=TRUE) 
    
    # Generate plot
    pdf(paste(paste(i[,"compound"][1], set_assay, "EC50", sep="-"),".pdf",sep=""))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    # If model converged, plot the model
    if (!is(m1, "try-error")) {
      
      #print(summary(m1))
      
      # Check if model is an OK fit of data; if it is, find values and plot fit
      # Currently changed to zero because we are accepting all model fits
      fit_test <- modelFit(m1)
      #print(fit_test)
      if (fit_test[,"p value"][2] > 0) {
        
        # Calculate EC30 and EC50 values and associated standard errors
        cutoff_doses <- ED(m1, c(30,50), interval="delta", display=FALSE)
        ec30_dose <- round(cutoff_doses[1,1], digits=2)
        ec30_se <- round(cutoff_doses[1,2], digits=2)
        ec50_dose <- round(cutoff_doses[2,1], digits=2)
        ec50_se <- round(cutoff_doses[2,2], digits=2)
        ec50_low <- ec50_dose - ec50_se
        ec50_hi <- ec50_dose + ec50_se
        upper_limit <- 100
        
        # If EC30 is above maximum tested dose, ignore it
        if (ec30_dose > max(x)) {
          ec30_dose <- paste(">", max(x), sep="")
        }
        
        # If EC50 is above maximum tested dose, ignore it
        if (ec50_dose > max(x)) {
          ec50_dose <- paste(">", max(x), sep="")
        }     
 
        # If EC50 is below minimum tested dose, set to minimum dose
        if (ec50_dose < min(x[x>0])) {
          ec50_dose <- min(x[x>0])
        }
        
        
      plot(m1, lwd=2, pch="", xlab=expression(paste("Concentration (", mu,"M)")), ylab="% Viability", ylim=c(0,200), cex.lab=1.5, cex.axis=1.25, col="blue")
      points(cmp_means[,"dose"], cmp_means[,"value"], pch=3, cex=2, col="black")
      points(x, y, pch=19, cex=1.2, col="blue")
      title(main=paste(compound, set_assay, sep=" - "), cex.main=1.5)
      abline(h=.50*upper_limit, col="grey70", lwd=2, lty=2)
      legend("topright", paste("EC50 = ", ec50_dose, "µM"), col="red", cex=1.2)
      }
      
      # If model failed to converge, plot data points and set EC50 to above max dose  
    } else {
      EC50 <- max(x)
      plot(cmp_means[,"dose"], cmp_means[,"value"], pch=3, cex=2, xlab=expression(paste("Concentration (", mu,"M)")), ylab="% viability", ylim=c(0,200), log="x", cex.lab=1.5, cex.axis=1.25, col="black")
      points(x, y, pch=19, cex=1.2, col="blue")
      title(main=paste(compound, set_assay, sep=" - "), cex.main=1.5)
      legend("topright", paste("EC50 > ", EC50, "µM"), col="red", cex=1.2)
    }  
    
    dev.off()
    
    # If model converged, print EC50 value and associated std error
    if (!is(m1, "try-error")) {
      cat(paste(compound, set_assay, ec50_dose, ec50_se, ec50_low, ec50_hi, sep="\t"))
      
      # If not, set to greater than max dose
    } else { 
      cat(paste(compound, set_assay, "NA", "NA", "NA", "NA", sep="\t"))
    }
    
    
    # If hill model was not attempted, plot data points and set EC50 to above max dose  
  } else {
    
    pdf(paste(paste(i[,"compound"][1], set_assay, "EC50", sep="-"),".pdf",sep=""))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    EC50 <- max(x)
    plot(cmp_means[,"dose"], cmp_means[,"value"], pch=3, cex=2, xlab=expression(paste("Concentration (", mu,"M)")), ylab="% viability", ylim=c(0,200), log="x", cex.lab=1.5, cex.axis=1.25, col="black")
    points(x, y, pch=19, cex=1.2, col="blue")
    title(main=paste(compound, set_assay, sep=" - "), cex.main=1.5)
    legend("topright", paste("EC50 > ", EC50, "µM"), col="red", cex=1.2)
    
    dev.off()
    
    # Print EC50 value not found
    cat(paste(compound, set_assay, "NA", "NA", "NA", "NA", sep="\t"))
  }
cat("\n") # add blank line after each entry
}
sink()





