#This script will produce a Venn Diagram for unique interactions  

# clearAndLoadPackages ----------------------------------------------------
rm(list=ls())

#Check which packages are installed and load all needed packages
packages <- c("magrittr", "ggplot2", "reshape2", "readr", "plyr", "xlsx")

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
filePaths <- c(
  "../refNetwork_ZeiglerFibroblast.csv",
  "../refNetwork_TanMechanosignaling.csv",
  "../refNetwork_RyallHypertrophy.csv"
  )

setwd("Y:/Archived/Matt Van de Graaf/Hypertrophy/RyallNetworkComparison/downloadedDatabases/")

#Get the numbers for each section of the Venn Diagram ---- 

Ryall <- read.csv(filePaths[3], colClasses = "character")
Ryall <- paste0(Ryall$inGene, Ryall$outGene)

Tan <- read.csv(filePaths[2], colClasses = "character")
Tan <- paste0(Tan$inGene, Tan$outGene)

Zeigler <- read.csv(filePaths[1], colClasses = "character")
Zeigler <- paste0(Zeigler$inGene, Zeigler$outGene)


RyallTan <- intersect(Ryall,Tan)

RyallTanZeigler <- intersect(RyallTan, Zeigler)
RyallTanZeiglerNum <- length(RyallTanZeigler)

RyallTanNum <- length(RyallTan) - RyallTanZeiglerNum 

RyallZeigler <- intersect(Ryall, Zeigler)
RyallZeiglerNum <- length(RyallZeigler) - RyallTanZeiglerNum

RyallNum <- length(Ryall) - RyallZeiglerNum - RyallTanNum - RyallTanZeiglerNum

TanZeigler <- intersect(Tan, Zeigler)
TanZeiglerNum <- length(TanZeigler) - RyallTanZeiglerNum

TanNum <- length(Tan) - TanZeiglerNum - RyallTanNum - RyallTanZeiglerNum

ZeiglerNum <- length(Zeigler) - TanZeiglerNum - RyallZeiglerNum - RyallTanZeiglerNum

