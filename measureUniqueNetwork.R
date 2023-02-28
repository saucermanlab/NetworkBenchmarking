#This script will measure unique nodes and interactions in each network 

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
filePaths <- c(#"../KEGG_CardiacMuscleContraction.csv", 
  #"../KEGG_HypertrophicCardiomyopathy.csv",
  #"../KEGG_PathwaysInCancer.csv",
  "../refNetwork_ZeiglerFibroblast.csv",
  "../refNetwork_TanMechanosignaling.csv",
  "../refNetwork_RyallHypertrophy.csv"
)

setwd("Z:/Matt Van de Graaf/Hypertrophy/RyallNetworkComparison/downloadedDatabases/")

allNodes <- as.character(array(dim = 0))
allInts <- as.character(array(dim = 0))


for (i in 1:length(filePaths)) {
  thisFile <- read.csv(filePaths[i], colClasses = "character")
  nodes <- unique(c(thisFile$inGene, thisFile$outGene))
  indices <- seq(length(allNodes) + 1, length(nodes) + length(allNodes))
  
  allNodes[indices] <- nodes
  
  interactions <- unique(paste(thisFile$inGene, thisFile$outGene))
  indices <- seq(length(allInts) + 1, length(interactions) + length(allInts))
  
  allInts[indices] <- interactions
}

genesForCompare <- data.frame(table(allNodes))
intsForCompare <- data.frame(table(allInts))

uniqueMeasures <- data.frame(matrix(nrow = length(filePaths), ncol = 3))

colnames(uniqueMeasures) <- c("Network", 
                              "Genes",
                              "Interactions")
                              

uniqueMeasures$Network <- gsub(".csv", "", gsub("../", "", filePaths))


for (i in 1:length(filePaths)) {
  thisFile <- read.csv(filePaths[i], colClasses = "character")
  nodes <- unique(c(thisFile$inGene, thisFile$outGene))
  indices <- match(nodes, genesForCompare$allNodes, nomatch = 0)
  uniqueGenes <- which(genesForCompare$Freq[indices] == 1)
  
  uniqueMeasures[i,2] <- (length(uniqueGenes)) /  length(nodes)
  
  interactions <- unique(paste(thisFile$inGene, thisFile$outGene))
  indices <- match(interactions, intsForCompare$allInts, nomatch = 0)
  uniqueInts <- which(intsForCompare$Freq[indices] == 1)
  
  uniqueMeasures[i,3] <- (length(uniqueInts)) /  length(interactions)
}

uniqueMeasures$Network <- c("Zeigler", "Tan", "Ryall")

uniqueGraph <- melt(uniqueMeasures, id.vars = "Network")

uniquePlot <- ggplot(uniqueGraph, 
                  aes(x = Network, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip(ylim = c(0,1)) +
  labs(y = "Fraction Unique To Network", 
       x = "Network",
       title = "Comparing Network Genes and Interactions") +
  scale_fill_grey(name="Measures:")  +
  theme(text = element_text(size=24), 
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 24))


ggsave("uniquePlot.png", uniquePlot, width = 12, height = 5, units = "in")