#merge all X2K databases

# clearAndLoadPackages ----------------------------------------------------
rm(list=ls())

#Check which packages are installed and load all needed packages
packages <- c("xlsx", "magrittr", "ggplot2", "reshape2", "readr")

check <- as.data.frame(installed.packages())
toBeInstalled <- is.na(match(packages, check$Package))

for (i in 1:length(toBeInstalled)) {
  if (toBeInstalled[i]) {
    install.packages(packages[i])
  }
}
lapply(packages, require, character.only = TRUE)
#end ----

# load files, define refNetworks, join tables ---------------------------------------
#Make a vector of the reaction files 

#Change working directory 
setwd("Z:/Matt Van de Graaf/Hypertrophy/RyallNetworkComparison/downloadedDatabases/X2K/")

toMerge <- list.files(path = "individualFiles/")

output = data.frame()

names <-  c("participantA", "NA.", "NA..1", "NA..2", "NA..3",
            "participantB", "NA..4", "NA..5", "NA..6", "NA..7",
            "reaction", "type", "pubMedID")
            
for (i in 1:length(toMerge)) {
  thisDatabase <- read.table(file = paste0("individualFiles/",toMerge[i]),
                             sep = "", col.names = names)
  
  thisDatabase$database <- strsplit(x = toMerge[i], split = ".sig")
  
  output <- rbind(output, thisDatabase)
}

output <- output[,c(1,6,11,12,13,14)]

write.table(x = as.matrix(output), file = "X2K_MergedDatabases.txt", 
            sep = "\t", row.names = FALSE)