#This script will evaluate the results across multiple databases 

# clearAndLoadPackages ----------------------------------------------------
rm(list=ls())

#Check which packages are installed and load all needed packages
packages <- c("openxlsx", "magrittr", "ggplot2", "reshape2", "plyr",
              "xlsx", "gridExtra", "extrafont")

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
fileNames <- c(#"KEGG_CardiacMuscleContraction", 
  #"KEGG_HypertrophicCardiomyopathy",
  #"KEGG_PathwaysInCancer",
  "refNetwork_ZeiglerFibroblast",
  "refNetwork_TanMechanosignaling",
  "refNetwork_RyallHypertrophy")
  #"refPath4_hypertrophy")

resultPaths <- c("Signor/signorData.xlsx",
                 "Reactome/ReactomeData.xlsx",
                 "pathwayCommonsData.xlsx",
                 "X2K/X2KData.xlsx",
                 "Omnipath/OmnipathData.xlsx")

databaseNames <-  c("Signor", "Reactome", "Pathway Commons", "X2K", "Omnipath", "Combined")

numsForGraphUndirected <- data.frame(matrix(ncol = length(fileNames) + 1, nrow = length(databaseNames)))
colnames(numsForGraphUndirected) <- c("Source", fileNames)
numsForGraphUndirected$Source <- databaseNames

numsForGraphDirected <- data.frame(matrix(ncol = length(fileNames) + 1, nrow = length(databaseNames)))
colnames(numsForGraphDirected) <- c("Source", fileNames)
numsForGraphDirected$Source <-databaseNames

numsForGraphAll <- data.frame(matrix(ncol = length(fileNames) + 1, nrow = length(databaseNames)))
colnames(numsForGraphAll) <- c("Source", fileNames)
numsForGraphAll$Source <-databaseNames

for (iData in 1:length(fileNames)) {
  
  resultOne<- read.xlsx2(file = resultPaths[1], sheetName = fileNames[iData])
  
  output <- data.frame(matrix(ncol = 6, nrow = nrow(resultOne)))
  
  outputColNames <- c(colnames(resultOne)[c(1,2,3, ncol(resultOne))], 'undirected', 'all')
  colnames(output) <- outputColNames
  output$name <- resultOne$name
  output$inGene <- resultOne$inGene
  output$outGene <- resultOne$outGene
  output$undirected <-as.numeric(as.character(resultOne$all))
  output$directed <- as.numeric(as.character(resultOne$directed))
  
  output[paste0(databaseNames[1], "Undirected")] <- resultOne$all
  output[paste0(databaseNames[1], "Directed")] <- resultOne$directed
  output[paste0(databaseNames[1], "All")] <- as.numeric(as.character(resultOne$all)) + 
    as.numeric(as.character(resultOne$directed))
  
  rxns <- unique(resultOne$name)
  score <- 0
  for (iRxn in 1:(length(rxns) - 1)) {
    name <- rxns[iRxn]
    thisRxn <- output[paste0(databaseNames[1], "All")][as.numeric(
      as.character(which(name == output$name))),]
    if (sum(thisRxn) > 0) {
      score <- score + 1
    }
  }
  output[paste0(databaseNames[1], "All")][nrow(output),] <- 
    score/(length(rxns) - 1)
  numsForGraphUndirected[1,1 + iData] <- output$undirected[nrow(output)]
  numsForGraphDirected[1,1 + iData] <- output$directed[nrow(output)]
  numsForGraphAll[1,1 + iData] <- score/(length(rxns) - 1)
  
  for (jResult in 2:length(resultPaths)) {
    thisResult <- read.xlsx2(file = resultPaths[jResult], sheetName = fileNames[iData])
    
    output$undirected <- output$undirected + as.numeric(as.character(thisResult$all_all))
    output$directed <- output$directed + as.numeric(as.character(thisResult$all_directed))
    output$all <- output$all + as.numeric(as.character(thisResult$all_all)) + 
      as.numeric(as.character(thisResult$all_directed))
    
    output[paste0(databaseNames[jResult], "Undirected")] <- thisResult$all_all
    output[paste0(databaseNames[jResult], "Directed")] <- thisResult$all_directed
    output[paste0(databaseNames[jResult], "All")] <- as.numeric(as.character(thisResult$all_all)) + 
      as.numeric(as.character(thisResult$all_directed))
    
    rxns <- unique(thisResult$name)
    score <- 0
    for (iRxn in 1:(length(rxns) - 1)) {
      name <- rxns[iRxn]
      thisRxn <- output[paste0(databaseNames[jResult], "All")][as.numeric(
        as.character(which(name == output$name))),]
      if (sum(thisRxn) > 0) {
        score <- score + 1
      }
    }
    numsForGraphUndirected[jResult,1 + iData] <-  as.numeric(as.character(
      thisResult$all_all[nrow(thisResult)]))
    numsForGraphDirected[jResult,1 + iData] <- as.numeric(as.character(
      thisResult$all_directed[nrow(thisResult)]))
    
    output[paste0(databaseNames[jResult], "All")][nrow(output),] <- 
      score/(length(rxns) - 1)
    numsForGraphAll[jResult,1 + iData] <- score/(length(rxns) - 1)
  }
  
  output$all <- output$directed + output$undirected
  
  numPresentUndirected <- 0
  numPresentDirected <- 0
  numPresentAll <- 0
  
  rxnNames <- unique(output$name)
  
  for (iName in 1:(length(rxnNames) - 1)) {
    thisRxn <- (rxnNames[iName] == output$name)
    
    rxnNums <- output[thisRxn,]
    
    if (sum(rxnNums$all) > 0) {
      numPresentAll <- numPresentAll + 1
    }
    if (sum(rxnNums$directed) > 0) {
      numPresentDirected <- numPresentDirected + 1
    }
    if (sum(rxnNums$undirected) > 0) {
      numPresentUndirected <- numPresentUndirected + 1
    }
  }
  output$undirected[nrow(output)] <- numPresentUndirected/(length(rxnNames) - 1)
  output$directed[nrow(output)] <- numPresentDirected/(length(rxnNames) - 1)
  output$all[nrow(output)] <- numPresentAll/(length(rxnNames) - 1)
  
  numsForGraphUndirected[nrow(numsForGraphUndirected),1 + iData] <-  
    numPresentUndirected/(length(rxnNames) - 1)
  numsForGraphDirected[nrow(numsForGraphDirected),1 + iData] <-  
    numPresentDirected/(length(rxnNames) - 1)
  numsForGraphAll[nrow(numsForGraphAll),1 + iData] <-  
    numPresentAll/(length(rxnNames) - 1)
  
  write.xlsx2(x = output, file = "compareOverlap.xlsx",
              sheetName = fileNames[iData], append = TRUE, row.names = FALSE)
}
#Add excel sheets for both the directed and all versions 

write.xlsx2(x = numsForGraphUndirected, file = "compareOverlap.xlsx",
            sheetName = "summaryAll", append = TRUE, row.names = FALSE)

write.xlsx2(x = numsForGraphDirected, file = "compareOverlap.xlsx",
            sheetName = "summaryDirected", append = TRUE, row.names = FALSE)


#Making plots for the combined data comparing it to individual data sets
#All Plot-------------------
colnames(numsForGraphAll) <- c("Source", "Cardiac\nFibroblast", "Mechano-\nSignaling",
                                    "Cardiac\nHypertophy")
dataAll.m <- melt(numsForGraphAll, id.vars = 'Source')

emptyAll <- dataAll.m$value == 0
dataAll.m$value[emptyAll] <- 0.005

allPlot <- ggplot(dataAll.m, aes(variable, value)) +
  geom_bar(aes(fill = Source), position = "dodge", stat = "identity") +
  coord_flip(ylim = c(0,1)) +
  labs(y = "", 
       x = "",
       title = "All Interactions") +
  scale_fill_grey(guide = FALSE)  +
  theme(text = element_text(size=20, family = "Arial")) 

ggsave("compareOverlapAll.png", allPlot, width = 10, 
       height = 5, units = "in")

####Make directed plot ###--------------------
colnames(numsForGraphDirected) <- c("Source", "Cardiac\nFibroblast", "Mechano-\nSignaling",
                                    "Cardiac\nHypertophy")
dataDirected.m <- melt(numsForGraphDirected, id.vars = 'Source')

emptyDirected <- dataDirected.m$value == 0
dataDirected.m$value[emptyDirected] <- 0.005

directedPlot <- ggplot(dataDirected.m, aes(variable, value)) +
  geom_bar(aes(fill = Source), position = "dodge", stat = "identity") +
  coord_flip(ylim = c(0,1)) +
  #ylim(0, 1) + 
  labs(y = "", 
       x = "Curated Network Model",
       title = "Directed Interactions") +
  scale_fill_grey(name="Database")  +
  theme(text = element_text(size=20, family = "Arial"), 
        legend.justification = c(1,.5), legend.position=c(.95,.45)) +
  guides(fill = guide_legend(reverse = TRUE))

ggsave("compareOverlapDirected.png", directedPlot, width = 10, 
       height = 5, units = "in")
#Undirected-------------------------
colnames(numsForGraphUndirected) <- c("Source", "Cardiac\nFibroblast", "Mechano-\nSignaling",
                                      "Cardiac\nHypertophy")
dataUndirected.m <- melt(numsForGraphUndirected, id.vars = 'Source')

emptyUndirected <- dataUndirected.m$value == 0
dataUndirected.m$value[emptyUndirected] <- 0.005

undirectedPlot <- ggplot(dataUndirected.m, aes(variable, value)) +
  geom_bar(aes(fill = Source), position = "dodge", stat = "identity") +
  coord_flip(ylim = c(0,1)) +
  labs(y = "Fraction of Reactions Present in Database", 
       x = "",
       title = "Undirected Interactions") +
  scale_fill_grey(guide = FALSE)  +
  theme(text = element_text(size=20, family = "Arial"))

ggsave("compareOverlapUndirected.png", undirectedPlot, width = 10, 
       height = 5, units = "in")

#make one large figure -----------

combinedPlot <- grid.arrange(allPlot, directedPlot, undirectedPlot, ncol = 1)
ggsave("compareOverlapCombined.png", combinedPlot, width = 12, 
       height = 10, units = "in")
