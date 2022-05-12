#' Find Enrichment of GO Terms in PatternMarker Peaks using GREAT
#'
#' Use the rGREAT package to find enrichment of GO terms or genes for
#' the peaks found to be most pattern differentiating using the
#' PatternMarker statistic.
#'
#' @param cogapsResult result object from CoGAPS
#' @param granges GRanges object corresponding to the peaks of the scATAC-seq
#'   data CoGAPS was applied to
#' @param genome UCSC genome designation for input to the sumbitGreatJob
#'   function from the rGREAT package (e.g. "hg19")
#' @param scoreThreshold threshold of PatternMarker score to take peaks for
#'   analysis, higher values return more peaks. Defaults to use all
#'   PatternMarker genes with value NULL
#' @param GREATCategory input to the category argument of the rGREAT 
#' getEnrichmentTables function. Usually "GO" or "Genes"
#' @return list containing enrichment results for each pattern
#' @examples data("schepCogapsResult")
#' data(schepGranges)
#'
#' GOenrichment <- applyGREAT(cogapsResult = schepCogapsResult,
#'  granges = schepGranges, genome = "hg19")
#' @export


applyGREAT <- function(cogapsResult, granges, genome,
                       scoreThreshold = NULL, GREATCategory = "GO") {
  
  #get PatternMarker peak indices
  patMarkers = CoGAPS::patternMarkers(cogapsResult)
  if(is.null(scoreThreshold)){
    peaks = patMarkers$PatternMarkers
    nPeaks = lapply(peaks, length)
    PMRanks = patMarkers$PatternMarkerRanks
    regionPatList = vector(mode = "list", length = ncol(PMRanks))
    for(i in seq(ncol(PMRanks))) {
      topPeaksPat = which(PMRanks[,i] %in% seq(nPeaks[[i]]))
      regionPatList[[i]] = topPeaksPat
    }
  }
  else{
    patScores = as.data.frame(patMarkers$PatternMarkerScores)
    chr_regions = rownames(patScores)
    regionPatList = vector(mode=  "list", length = ncol(patScores))
    for(i in seq(ncol(patScores))) {
      topPeaksPat = which(patScores[,i] < scoreThreshold)
      regionPatList[[i]] = topPeaksPat
    }
  }
  
  
  #print number of peaks used based on patternMarker score threshold
  numPeaks = unlist(lapply(regionPatList, length))
  names(numPeaks) = lapply(seq(length(regionPatList)),
                           function(x) {paste("Pattern", x)})
  print("Number of peaks used for each pattern:", quote = FALSE)
  print(numPeaks)
  
  #call GREAT and get enrichment for GO terms
  GREATResults = vector("list", length(regionPatList))
  for(i in seq(length(regionPatList))) {
    patRanges = granges[regionPatList[[i]]]
    
    greatJob = rGREAT::submitGreatJob(gr = patRanges, species = genome,
                                      request_interval = 1)
    tbl = rGREAT::getEnrichmentTables(greatJob, category = GREATCategory)
    
    GREATResults[[i]] = tbl
  }
  
  return(GREATResults)
}
