#Make a figure for benchmarking performance for the interactions 
#that are common to all networks

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
#end ----

# load files, define refNetworks, join tables ---------------------------------------
#Make a vector of the reaction files 
fileNames <- c(
  "refNetwork_ZeiglerFibroblast",
  "refNetwork_TanMechanosignaling",
  "refNetwork_RyallHypertrophy")

resultPaths <- c("Signor/signorData.xlsx",
                 "Reactome/ReactomeData.xlsx",
                 "pathwayCommonsData.xlsx",
                 "X2K/X2KData.xlsx")

databaseNames <-  c("Signor", "Reactome", "Pathway Commons", "X2K")

numsForGraph <- data.frame(matrix(ncol = 3, nrow = length(databaseNames)))
colnames(numsForGraph) <- c("Source", "Undirected", "Directed")
numsForGraph$Source <- databaseNames


#Read in the results
getCommonInts <- function(network1, network2, network3){
  network1[-nrow(network1),]
  network2[-nrow(network2),]
  network3[-nrow(network3),]
  
  ints1 <- paste0(network1$inGene,network1$outGene)
  ints2 <- paste0(network2$inGene,network2$outGene)
  overlapNetworkIndex <- intersect(ints1,ints2) %>% 
    match(ints1) 
  overlapNetwork <- network1[overlapNetworkIndex,]
  
  intsOverlap <- paste0(overlapNetwork$inGene, overlapNetwork$outGene)
  ints3 <- paste0(network3$inGene,network3$outGene)
  overlapNetworkIndex <- intersect(intsOverlap,ints3) %>% 
    match(intsOverlap) 
  overlapNetwork <- overlapNetwork[overlapNetworkIndex,]  
  
  return(overlapNetwork)
}
#need to go to a single database and get all 3 resuls files

for (i in 1:length(resultPaths)) {
  zeigler = read.xlsx2(file = resultPaths[i], sheetName = fileNames[1])
  tan = read.xlsx2(file = resultPaths[i], sheetName = fileNames[2])
  ryall = read.xlsx2(file = resultPaths[i], sheetName = fileNames[3])
  
  if (i == 1) {
    colnames(zeigler) <- c("name", "inGene", "outGene",
                           "all_all", "all_directed")
    colnames(tan) <- c("name", "inGene", "outGene",
                           "all_all", "all_directed")
    colnames(ryall) <- c("name", "inGene", "outGene",
                           "all_all", "all_directed")
  }
  
  overlapNetwork <- getCommonInts(ryall, tan, zeigler)
  numsForGraph$Undirected[i] <- 
    sum(as.numeric(as.character(
      overlapNetwork$all_all)))/(nrow(overlapNetwork) - 1)
  numsForGraph$Directed[i] <- 
    sum(as.numeric(as.character(
      overlapNetwork$all_directed)))/(nrow(overlapNetwork) - 1)
}

#Plot results 
data <- melt(numsForGraph)

commonPlot <- ggplot(data, aes(fill=variable, y=value, x=Source)) + 
  geom_bar(position="dodge", stat="identity")  +
  ylim(0, 1) + 
  labs(y = "Fraction of Reactions Present in Database", 
     x = "Database",
     title = "Benchmarking Results for Interactions Common to 3 Networks") +
  scale_fill_grey(name="Interaction Type")  +
  theme(text = element_text(size=20))+
  guides(fill = guide_legend(reverse = FALSE)) 

ggsave("benchmarkCommonInts.png", commonPlot, width = 12, 
       height = 10, units = "in")
