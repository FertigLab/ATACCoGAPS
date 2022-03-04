#' Match cells to patterns
#'
#' Use the patternMarker statistic to determine which cells belong to each
#' pattern in the data
#' 
#' @param cgapsResult a CoGAPSResult object
#' @return list containing a prediction matrix and vector classifying cells to patterns
#' @examples repmis::source_data("https://github.com/FertigLab/ATACCoGAPS/blob/PaperVersion/data/schepCogapsResult.rda?raw=true")
#' pClass <- patternMarkerClassifier(schepCogapsResult)
#' @export
patternMarkerCellClassifier <- function(cgapsResult) {
  
  #get PatternMarker cells
  PM <- CoGAPS::patternMarkers(cgapsResult, axis = 2)
  patMarkers <- PM$PatternMarkers
  #get number of patternMarker cells found for each pattern
  lengths <- lapply(patMarkers, length)
  #get PM ranks
  patRanks <- PM$PatternMarkerRanks
  
  #create matrix of cells predicted to be in each pattern
  predictionMatrix <- matrix(NA, nrow(patRanks), ncol(patRanks)) 
  for(i in seq(ncol(patRanks))){
    prediction <- rep(0, nrow(patRanks))
    preds <- which(patRanks[,i] %in% seq(lengths[[i]]))
    prediction[preds] <- 1
    predictionMatrix[,i] <- prediction
  }
  
  #create list of which pattern each cell belongs to, with unclassified cells labelled 0
  classList <- apply(predictionMatrix, 2, function(x){
    cellselect <- which(x == 1)
  })
  cellClassifier <- vector("numeric", nrow(predictionMatrix))
  for(i in seq(ncol(predictionMatrix))){
    cellClassifier[classList[[i]]] <- i
  }
  #return both prediction Matrix (for benchmarking against ground truth) and
  #classifier (for downstream analysis)
  return(list(predictionMatrix = predictionMatrix, cellClassifier = cellClassifier))
}
