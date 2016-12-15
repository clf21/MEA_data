## Batch converting spike counts .csv files in a selected directory to format compatible with DARB.exe program.
##
## Chris Frank - May 2016



## Prompt user to select directory with spike counts .csv files
working_dir <- choose.dir(caption="Please select directory with spike counts .csv files")
setwd(working_dir)


## Loop through each spike counts .csv file
for (i in list.files(pattern="spike_counts.csv")) {
  
  # Ignore the file name if it already has "DARB" or "results" in the name
  if (grepl("DARB", i) || grepl("results", i)) {
    print(paste(i, " looks like it's already been processed, skipping!", sep=""))
    print("_________________________________________________________________________")
  
  } else {
    print(paste("Processing file ", i, sep=""))
    current_file <- read.delim(i, sep=",", header=F, colClasses="character")
    current_file <- current_file[,c(3:4,55:ncol(current_file))] # Subset retained columns
    current_file <- current_file[1:(nrow(current_file)-11),] # Subset retained rows
    current_file[1,] <- gsub("_", " ", current_file[1,]) # Replace all underscores with spaces in first row
    
    outfile <- paste(strsplit(i, "\\.")[[1]][1], "_DARBr.", strsplit(i, "\\.")[[1]][2], sep="")
    print(paste("Writing to ", outfile, sep=""))
    write.table(current_file, file=outfile, sep=",", quote=F, row.names=F, col.names=F) # Output reformatted table, remains comma-delim
    print("_________________________________________________________________________")
  }
  
}


## Open up the DARB.exe program for the user
print("Done!  Opening DARB...")
print("Use the ...DARBr.csv files as input to DARB.exe")
system2("L:/Lab//NHEERL_MEA/MAESTRO SYSTEM/Axions Data Processing/D.A.R.B.exe", invisible=FALSE)
