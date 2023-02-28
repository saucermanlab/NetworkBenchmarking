#This script has functions to be used in pathwayCommonsData.R

networkRecall <- function(filePath, allDatabases, individualDatabases){
  #Adjust pathway data frame ----
  thisPathway <- read.csv(filePath)
  addColNames <- paste(individualDatabases, c("directedOnly"), sep = "_") %>%
    cbind(paste(individualDatabases, c("all"), sep = "_")) %>%
    sort() %>%
    c( "all_all", "all_directed")
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
    databaseIndexes <-  unique(directedMatchRxns$database) %>%
      sort() %>%
      match(individualDatabases) * 2 - 1
    databaseIndexes <- sort(c(databaseIndexes + 1)) + 
      length(pathwayNames)
    
    thisPathway[iPathway,(databaseIndexes)] <- 1;

    #### undirectedMatches ####
    undirectedMatches <- 
      (tolower(undirectedInts$inGene) == tolower(thisPathway$inGene[iPathway]) & 
        tolower(undirectedInts$outGene) == tolower(thisPathway$outGene[iPathway])) |
      (tolower(undirectedInts$inGene) == tolower(thisPathway$outGene[iPathway]) & 
        tolower(undirectedInts$outGene) == tolower(thisPathway$inGene[iPathway])) 
    
    undirectedMatchRxns <- undirectedInts[undirectedMatches,]    
    databaseIndexes <-  unique(undirectedMatchRxns$database) %>%
      sort() %>%
      match(individualDatabases) * 2 - 1
    databaseIndexes <- databaseIndexes + length(pathwayNames)
    
    thisPathway[iPathway,databaseIndexes] <- 1;
    #### end ####
    
    #### calculate totals for each row ####
    if (sum(thisPathway[iPathway, seq(length(pathwayNames) + 2, 
                                      ncol(thisPathway) - 2, 2)]) > 0) {
      thisPathway$all_directed[iPathway] <- 1
    }
    if (sum(thisPathway[iPathway, seq(length(pathwayNames) + 1, 
                                      ncol(thisPathway) - 2, 2)]) > 0) {
      thisPathway$all_all[iPathway] <- 1
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