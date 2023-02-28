#This script has functions to be used in pathwayCommonsData.R

networkRecall <- function(filePath, allDatabases){
  #Adjust pathway data frame ----
  thisPathway <- read.csv(filePath)
  addColNames <- c( "all_all", "all_directed")
  #Add zero columns 
  pathwayNames <- colnames(thisPathway)
  thisPathway[,(ncol(thisPathway) + 1):(ncol(thisPathway) + length(addColNames))] <- 
    matrix(data = 0, nrow = nrow(thisPathway), ncol = length(addColNames))
  colnames(thisPathway) <- c(pathwayNames, t(addColNames)) 
  #### end #####
  #### Break data structure into directed and undirected ####
  directedInts <- allDatabases[allDatabases$DIRECTED == 1,]
  undirectedInts <- allDatabases[allDatabases$DIRECTED == 0,]
  for (iPathway in 1:nrow(thisPathway)) {
    
    ##### Directed matches first ####
    directedMatches <- tolower(directedInts$inGene) == tolower(thisPathway$inGene[iPathway]) & 
      tolower(directedInts$outGene) == tolower(thisPathway$outGene[iPathway])
    
    directedMatchRxns <- directedInts[directedMatches,]

    if (nrow(directedMatchRxns) > 0) {
      thisPathway[iPathway,c(5)] <- 1;
    }

    #### undirectedMatches ####
    undirectedMatches <- 
      (tolower(undirectedInts$inGene) == tolower(thisPathway$inGene[iPathway]) & 
        tolower(undirectedInts$outGene) == tolower(thisPathway$outGene[iPathway])) |
      (tolower(undirectedInts$inGene) == tolower(thisPathway$outGene[iPathway]) & 
        tolower(undirectedInts$outGene) == tolower(thisPathway$inGene[iPathway])) 
    
    undirectedMatchRxns <- undirectedInts[undirectedMatches,]    
    
    if (nrow(undirectedMatchRxns) > 0) {
      thisPathway[iPathway,4] <- 1;
    }
    #### end ####
    }
  #### end ####
  #### Calculate recall for each column ####
  
  rxns <- unique(thisPathway$name)
  
  thisPathway[nrow(thisPathway) + 1, ] <- matrix(0, nrow = 1, ncol = ncol(thisPathway))
  
    for (kRxn in 1:length(rxns)) {
      index <-  rxns[kRxn] == thisPathway$name
      rxnSum <- colSums(
        thisPathway[index, (length(pathwayNames) + 1):ncol(thisPathway)],
        na.rm = TRUE)
      if (sum(rxnSum > 0)) {
        thisPathway[nrow(thisPathway),c(FALSE, FALSE, FALSE,rxnSum != 0 )] <- 
          thisPathway[nrow(thisPathway),c(FALSE, FALSE, FALSE,rxnSum != 0 )] + 1
      } 

    }
  
  thisPathway[nrow(thisPathway),] <- 
    as.numeric(thisPathway[nrow(thisPathway), ]) / length(rxns)
  
  
  
  return(thisPathway)
}