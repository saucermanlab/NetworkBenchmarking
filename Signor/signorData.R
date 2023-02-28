#This script will evaluate the results of using the Signor database

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

#Change working directory 
setwd("Z:/Archived/Matt Van de Graaf/Hypertrophy/RyallNetworkComparison/downloadedDatabases/Signor/")
 
source("signorFunctions.R")

#end ----

# load files, define refNetworks, join tables ---------------------------------------
#Make a vector of the reaction files 
filePaths <- c(#"../../KEGG_CardiacMuscleContraction.csv", 
  #"../../KEGG_HypertrophicCardiomyopathy.csv",
  #"../../KEGG_PathwaysInCancer.csv",
  "../../refNetwork_ZeiglerFibroblast.csv",
  "../../refNetwork_TanMechanosignaling.csv",
  "../../refNetwork_RyallHypertrophy.csv",
  "../../refPath4_hypertrophy.csv")

database <- read_tsv(file = "SignorDatabases.txt")

interactionKey <- read.xlsx("InteractionTypeKey.xlsx",1)

endResult <- data.frame(matrix(ncol = length(filePaths), 
                               nrow = 1))
colnames(endResult)[1] <- "SignorData"

#Join the directed column of interactionKey to allDatabases 
mergedWithKey <- join(x = database, y = interactionKey, "EFFECT", type = "inner")

mergedWithKey$inGene <- as.character(NA)
mergedWithKey$outGene <- as.character(NA)

#add inGene and outGene parameters to the databass file 
#Step 1: for ones that are undirected, assign inGene and outGene and then flip 
#   the in and outs and concatenate to the bottom of the allDatabases data frame 
noFlip <- mergedWithKey$FLIP == 0
mergedWithKey$inGene[noFlip] = mergedWithKey$ENTITYA[noFlip]
mergedWithKey$outGene[noFlip] = mergedWithKey$ENTITYB[noFlip]

Flip <- mergedWithKey$FLIP == 1
mergedWithKey$inGene[Flip] = mergedWithKey$ENTITYB[Flip]
mergedWithKey$outGene[Flip] = mergedWithKey$ENTITYA[Flip]

for (i in 1:length(filePaths)) {

  thisPathway <- networkRecall(filePaths[i], mergedWithKey)
  
  write.xlsx2(file = "signorData.xlsx", x = thisPathway, 
                            sheetName = gsub(".csv", "", gsub("../", "", filePaths[i])), 
                            append = TRUE, col.names = TRUE, row.names = FALSE)
}

