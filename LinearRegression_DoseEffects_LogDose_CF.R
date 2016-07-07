# MEA neuronal ontogeny, nested linear regression
# This determines whether dose term has a significant improvement on linear model
# Chris Frank - February-March 2016


# Set parameters of analysis:
set_directory <- "~/MEA_Data/PIP3_AIM2-3-DNT/" # Set directory 
set_data_file1 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/SA2_allData_11-17-2014.csv" # Name of csv file with burst parameters per sample
set_data_file2 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/AllData_SA3_NoDNTRef.csv"
set_data_file3 <- "~/MEA_Data/PIP3_AIM2-3-DNT/sourceData/DNTRef_all.csv"
set_compound <- "Cadmium" # This has to match compound in the data set. If output set to "all" below, this is ignored.
set_parameter <- "meanfiringrate"
set_output_file <- "LM_results.txt"
set_output <- "single_plot"  # Can be set to "single_plot" or "all" which does all combinations of compound and ontogeny parameter


# Load required packages
require(rgl)
require(car)


# Read in data and clean it up
setwd(set_directory)
all_data1 <- read.delim(set_data_file1, sep=",", stringsAsFactors=FALSE) # Reading in full data set ; assuming it is a csv table
all_data2 <- read.delim(set_data_file2, sep=",", stringsAsFactors=FALSE) # Reading in second data set
all_data3 <- read.delim(set_data_file3, sep=",", stringsAsFactors=FALSE) # Reading in third data set

extra_acet <- (all_data1$trt=="Acetaminophen" & all_data1$dose!=0) # Index extra acetaminophen treated wells
all_data1 <- all_data1[-which(extra_acet),] # Remove extra acetaminophen wells from Aim2 data
extra_acet <- (all_data3$trt=="Acetaminophen" & all_data3$dose!=0) # Index extra acetaminophen treated wells
all_data3 <- all_data3[-which(extra_acet),] # Remove extra acetaminophen wells from DNT Ref data


all_data <- rbind(all_data1, all_data2, all_data3)
rm(all_data1, all_data2, all_data3)

bis_rows <- grep("12_01_", all_data$file.name, fixed=TRUE) # Index all Bicuculline-treated wells
all_data <- all_data[- bis_rows,] # Remove all bic-treated wells
all_data[is.na(all_data)] <- 0 #Replace all NAs with zeros - This may be undesirable for MEA parameters that are derived from other parameters.
all_data[all_data[,"dose"]==0,"trt"] <- "Control" # Convert all zero dose rows to say "Control" for treatment
all_data <- all_data[all_data[,"DIV"]!=2,] # Remove all sparse DIV2 data



## Normalize data by % of control wells at DIV12 on same plate for a particular ontogeny parameter
normalize_by_plate <- function(ontogeny_param) {
  
  # first calculate median of control wells or any other duplicate treatment on each plate
  per_plate <- as.data.frame(do.call(cbind, aggregate(get(ontogeny_param) ~ date + Plate.SN + DIV + trt + dose, FUN=function(x) c(med=median(x), n=length(x)), data=all_data)), stringsAsFactors=FALSE)

  # Split data set by plate (date as well because plateIDs get reused)
  splt.by <- c('date','Plate.SN')
  per_plate_split <- split(per_plate, per_plate[,splt.by])
  per_plate_split2 <- per_plate_split[sapply(per_plate_split, function(x) dim(x)[1]) > 0]
  rm(per_plate_split)

  # Normalize each plate by percent of control median and combine into single data frame
  norm_plates <- list()
  for (i in per_plate_split2) {
    i[,6] <- (as.numeric(i[,6]) / as.numeric(subset(i, trt=="Control" & DIV==12)[1,6]))*100
    norm_plates[[length(norm_plates)+1]] <- i
  }

  norm_plates <- do.call(rbind, norm_plates) #Re-form one table of values
  class(norm_plates$DIV) <- "numeric" #make sure DIV and dose columns are considered numeric
  class(norm_plates$dose) <- "numeric"
  
  norm_plates  
}



## Generate linear models on normalized data and test whether fit significantly improves with dose term
make_lm <- function(compound, ontogeny_param, plot=TRUE){
  
  # Get normalized data
  norm_data <- normalize_by_plate(ontogeny_param)
  
  # Subset data to controls and compound of interest
  data_sub <- norm_data[which(norm_data$trt %in% c("Control", compound)),]
  
  # Add column of date plus plateID
  data_sub[,"id"] <- paste(data_sub[,1], data_sub[,2], sep="_")
  
  # Keep only control wells run on the same plate as treated
  data_sub <- data_sub[which(data_sub[,"id"] %in% subset(data_sub, trt==compound)[,"id"]),]
  
  # Log-transform dose
  data_sub[,"logdose"] <- log(data_sub[,"dose"] + 1)
  
  #print(data_sub)
  
  # First linear model without dose term ; using normalized data
  m1 <- lm(med ~ DIV, data=data_sub)
  
  # Second linear model with dose ; using normalized data
  m2 <- lm(med ~ DIV + logdose, data=data_sub)
  
  if (plot==TRUE) {
    # Aggregate to get mean of plate-medians for plotting
    data_avg <- aggregate(med ~ trt + DIV + dose, FUN=mean, data=data_sub)
    
    # Generate 3d scatterplot of mean values
    r3dDefaults$windowRect <- c(136,44,1218,953)
    plot3d(data_avg[,"DIV"], log(data_avg[,"dose"]+1), data_avg[,"med"], col="blue", size=8, xlab="DIV", ylab="log(dose)", zlab=ontogeny_param, main=compound, box=TRUE)
    #rgl.postscript("plot.pdf", "pdf") # To save the image, may not work on all systems
    #movie3d(spin3d(axis = c(0,0,1), rpm = 10), duration=6,  type = "png") # To create a movie of the image spinning
    #browseURL(paste("file://", writeWebGL(dir=file.path(tempdir(), "webGL"), width=500), sep=""))
    
    # Another 3d plotting option (car package)
    #scatter3d(med ~ DIV + log(dose+1), data=data_avg, axis.scales=TRUE, model.summary=TRUE)
    
    # Add regression plane to 3d plot 
    coefs <- coef(m2)
    a <- coefs["DIV"]
    b <- coefs["logdose"]
    c <- -1
    d <- coefs["(Intercept)"]
    rgl.planes(a, b, c, d, alpha=0.5, color = "steelblue")
    
    #rgl.viewpoint(zoom=0.8)
    #rgl.snapshot(filename=paste(compound, "_", ontogeny_param, "_1.png", sep=""))
    
  }

  
  # Test if model is significantly improved by dose term inclusion in second model
  results <- anova(m1, m2)
  
  exact_pval <- results$'Pr(>F)'[2] # get exact p-value from ANOVA
  
  print(paste(compound, ontogeny_param, sep=" - "))
  print(results)
  print(paste("exact p-value = ", exact_pval, sep=""))
  
  if(exists("num_tests")) {
    cor_pval <- p.adjust(exact_pval, method="bonferroni", n=num_tests) # correct p-values by Bonferroni method (multiply by number of comparisons made)
    # Note - in order to compute accurate FDR or other such corrected p-values, the complete vector of p-values must be supplied to p.adjust and let n be specified automatically
    print(paste("Bonferroni corrected p-value = ", cor_pval, sep=""))
  }
  
  print("-----------------------------------------------------")
  print("")
  
  
  return(exact_pval)
  
}






## Generate outputs

if (set_output=="single_plot"){
  make_lm(set_compound, set_parameter, plot=TRUE)
}




if (set_output=="all"){
  ont_names <- names(all_data[,8:(ncol(all_data)-3)]) # all parameter names, excluding cv.time and cv.network
  comp_names <- unique(all_data[,"trt"]) # get all compound names
  comp_names <- comp_names[!comp_names %in% ("Control")] # remove control from testing
  num_tests <- length(ont_names)*length(comp_names) # total number of statistical tests run
  pvals <- list() #initialize list of p-values
  q <- 1
  
  sink(file=set_output_file, split=TRUE) # write output to file
  
  # Loop through all compound and ontogeny parameter combinations, store p-values
  for (i in comp_names) {
    for (j in ont_names) {
      pvals[q] <- make_lm(i, j, plot=FALSE)
      q <- q+1
    }
  } 
  
  sink() # close main output file
  
  # Adjust p-values and organize tables
  adj_pvals <- p.adjust(pvals, method="BH")
  bf_pvals <- p.adjust(pvals, method="bonferroni")
  
  pval_list <- split(unlist(pvals), ceiling(seq_along(pvals)/length(ont_names)))
  pvals_df <- as.data.frame(lapply(pval_list, cbind), stringsAsFactors=FALSE)
  row.names(pvals_df) <- ont_names
  names(pvals_df) <- comp_names
  write.table(pvals_df, "lm_raw_pvalues.csv", sep=",", quote=FALSE)
  
  adj_pvals_list <- split(adj_pvals, ceiling(seq_along(adj_pvals)/length(ont_names)))
  adj_pvals_df <- as.data.frame(lapply(adj_pvals_list, cbind), stringsAsFactors=FALSE)
  row.names(adj_pvals_df) <- ont_names
  names(adj_pvals_df) <- comp_names
  write.table(adj_pvals_df, "lm_BHadjusted_pvalues.csv", sep=",", quote=FALSE)
  
  bf_pvals_list <- split(bf_pvals, ceiling(seq_along(bf_pvals)/length(ont_names)))
  bf_pvals_df <- as.data.frame(lapply(bf_pvals_list, cbind), stringsAsFactors=FALSE)
  row.names(bf_pvals_df) <- ont_names
  names(bf_pvals_df) <- comp_names
  write.table(bf_pvals_df, "lm_bonferroni_pvalues.csv", sep=",", quote=FALSE)
  
}

