# clearAndLoadPackages ----------------------------------------------------
rm(list=ls())

#Check which packages are installed and load all needed packages
packages <- c("openxlsx", "magrittr", "ggplot2", "reshape2", "plyr", "xlsx")

check <- as.data.frame(installed.packages())
toBeInstalled <- is.na(match(packages, check$Package))

for (i in 1:length(toBeInstalled)) {
  if (toBeInstalled[i]) {
    install.packages(packages[i])
  }
}
lapply(packages, require, character.only = TRUE)

source("Z:/Archived/Matt Van de Graaf/Hypertrophy/RyallNetworkComparison/downloadedDatabases/Signor/signorFunctions.R")

#end ----

# load files, define refNetworks, join tables ---------------------------------------
#Make a vector of the reaction files 
filePaths <- c(#"../../KEGG_CardiacMuscleContraction.csv", 
  #"../../KEGG_HypertrophicCardiomyopathy.csv",
  #"../../KEGG_PathwaysInCancer.csv",
  "refNetwork_ZeiglerFibroblast.csv",
  "refNetwork_TanMechanosignaling.csv",
  "refNetwork_RyallHypertrophy.csv",
  "refPath4_hypertrophy.csv")

database <- read_csv("OmnipathData.csv")
mergedWithKey <- database

mergedWithKey$inGene = mergedWithKey$source_genesymbol
mergedWithKey$outGene = mergedWithKey$target_genesymbol
mergedWithKey$Directed = mergedWithKey$is_directed

for (i in 1:length(filePaths)) {
  
  thisPathway <- networkRecall(filePaths[i], mergedWithKey)
  
  write.xlsx2(file = "Omnipath_Data.xlsx", x = thisPathway, 
              sheetName = gsub(".csv", "", gsub("../", "", filePaths[i])), 
              append = TRUE, col.names = TRUE, row.names = FALSE)
}
