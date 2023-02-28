#This script will generate a spreadsheet that has a tab for each interaction type in the 
#X2K downloaded spreadsheet

# clearAndLoadPackages ----------------------------------------------------
rm(list=ls())

#Check which packages are installed and load all needed packages
packages <- c("openxlsx", "magrittr", "ggplot2", "reshape2", "readr",
              "plyr")

check <- as.data.frame(installed.packages())
toBeInstalled <- is.na(match(packages, check$Package))

for (i in 1:length(toBeInstalled)) {
  if (toBeInstalled[i]) {
    install.packages(packages[i])
  }
}
lapply(packages, require, character.only = TRUE)
i#end ----

# load files, define refNetworks, join tables ---------------------------------------
#Make a vector of the reaction files 
filePaths <- c(#"../../KEGG_CardiacMuscleContraction.csv", 
               #"../../KEGG_HypertrophicCardiomyopathy.csv",
               #"../../KEGG_PathwaysInCancer.csv",
                "../../refNetwork_ZeiglerFibroblast.csv",
                "../../refNetwork_TanMechanosignaling.csv",
               "../../refNetwork_RyallHypertrophy_v2.csv",
               "../../refPath4_hypertrophy.csv")

#Change working directory 
setwd("Y:/Archived/Matt Van de Graaf/Hypertrophy/RyallNetworkComparison/downloadedDatabases/X2K/")

allDatabases <- read_tsv("X2K_MergedDatabases.txt", col_types = "ccccic")
interactionKey <- read.xlsx("InteractionTypeKey.xlsx", sheet = 1) #changed from sheetIndex to sheet

#Join the directed column of interactionKey to allDatabases 
allDatabases <- join(x = allDatabases, y = interactionKey, "reaction", type = "inner")

#Loop through each database and compare 
individualDatabases <- sort(unique(allDatabases$database))

#add inGene and outGene parameters to the databass file 
#Step 1: for ones that are undirected, assign inGene and outGene and then flip 
#   the in and outs and concatenate to the bottom of the allDatabases data frame 
noFlip <- allDatabases$FLIP == 0
allDatabases$inGene[noFlip] = allDatabases$participantA[noFlip]
allDatabases$outGene[noFlip] = allDatabases$participantB[noFlip]

Flip <- allDatabases$FLIP == 1
allDatabases$inGene[Flip] = allDatabases$participantB[Flip]
allDatabases$outGene[Flip] = allDatabases$participantA[Flip]
#end ----

# load function and call for processing ----
source("X2K_Functions.R")

endResult <- data.frame(matrix(ncol = length(filePaths), 
                               nrow = length(individualDatabases) *2 + 2))
colnames(endResult)[1] <- "databases"

for (iDatabase in 1:length(filePaths)) {
  thisPathway <- networkRecall(filePaths[iDatabase], 
                               allDatabases, individualDatabases) 
  write.xlsx2(file = "X2KData_v2.xlsx", x = thisPathway, 
              sheetName = gsub(".csv", "", gsub("../", "", filePaths[iDatabase])), 
              append = TRUE, col.names = TRUE, row.names = FALSE)
  
  endResult$databases <- colnames(thisPathway)[4:ncol(thisPathway)]
  endResult[,iDatabase + 1] <- t(thisPathway[nrow(thisPathway),
                                             4:ncol(thisPathway)])
  
  colnames(endResult)[1 + iDatabase] <- 
    gsub(".csv", "", gsub("../", "", filePaths[iDatabase]))
}

write.xlsx2(file = "X2KData_v2.xlsx", x = endResult, 
            sheetName = "allPaths", 
            append = TRUE, col.names = TRUE, row.names = FALSE)

# Split the results df into direccted and undirected for visualization ----

#Split the databases 
databaseNames <- c(individualDatabases, "All")

endResultsDir <- endResult[seq(2, nrow(endResult), by = 2),]
endResultsDir$databases <- databaseNames
endResultsAll <- endResult[seq(1, nrow(endResult) - 1, by = 2),]
endResultsAll$databases <- databaseNames

#Start visualization 
#directed 
endResultsDir <- melt(data = endResultsDir, id = c("databases"))

dirPlot <- ggplot(endResultsDir, 
          aes(x = reorder(databases, value), y = value, fill = variable)) +
          geom_bar(stat = "identity", position = "dodge") +
          coord_flip(ylim = c(0,1)) +
          labs(y = "Fraction of Reactions Present", 
               x = "Database",
               title = "Comparing X2K Databases (directed interactions only)")+
          scale_fill_discrete(name="Reference Network") +
          theme(text = element_text(size=20)) +
          guides(fill = guide_legend(reverse = TRUE))
  
ggsave("directedInteractionsPlot.png", dirPlot, width = 400,
       height = 400, units = "mm")
#undirected 
endResultsAll <- melt(data = endResultsAll, id = c("databases"))

allPlot <- ggplot(endResultsAll, 
          aes(x = reorder(databases, value), y = value, fill = variable)) +
          geom_bar(stat = "identity", position = "dodge") +
          coord_flip(ylim = c(0,1)) +
          labs(y = "Fraction of Reactions Present", 
                x = "Database",
                title = "Comparing X2K Databases") +
          scale_fill_discrete(name="Reference Network")  +
          theme(text = element_text(size=20))+
          guides(fill = guide_legend(reverse = TRUE))

ggsave("allInteractionsPlot.png", allPlot, width = 400, 
       height = 400, units = "mm")

#### end ####
#### summarize the interaction database counts ####
countList <- as.data.frame(table(sort(unique(allDatabases$database))), 
                                 stringsAsFactors=FALSE)

colnames(countList) <- c("Database", "Occurrences") 
countList <- rbind(countList, c("Total Num of ints", nrow(allDatabases)))

directedallDatabases <- allDatabases[allDatabases$DIRECTED == 1,]
directedCountList <- as.data.frame(table(sort(unique(directedallDatabases$database))),
                                                     stringsAsFactors=FALSE)

colnames(directedCountList) <- c("Database", "DirectedOccurrences") 
directedCountList <- rbind(directedCountList, c("Total Num of ints", nrow(directedallDatabases)))

countList <- join(countList, directedCountList, "Database", "left")

countList$percentDirected <-  as.numeric(countList$DirectedOccurrences) / 
  as.numeric(countList$Occurrences)

countList$percentOfTotalIns <- as.numeric(countList$Occurrences) /
  as.numeric(countList$Occurrences[nrow(countList)])

write.xlsx2(file = "X2KData_v2.xlsx", x = countList, 
            sheetName = "databaseNums", 
            append = TRUE, col.names = TRUE, row.names = FALSE)
