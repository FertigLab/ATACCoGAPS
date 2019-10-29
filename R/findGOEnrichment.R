#' Find Enrichment of GO Terms in PatternMarker Peaks
#'
#' Use the rGREAT package to find enrichment of GO terms for the peaks found to
#' be most pattern differentiating using the PatternMarker statistic.
#'
#' @param cogapsResult result object from CoGAPS
#' @param granges GRanges object corresponding to the peaks of the scATAC-seq
#'   data CoGAPS was applied to
#' @param genome UCSC genome designation for input to the sumbitGreatJob
#'   function from the rGREAT package (e.g. "hg19")
#' @param scoreThreshold threshold of PatternMarker score to take peaks for
#'   analysis, higher values return more peaks, default is 0.03
#' @return list containing GO enrichment result for each pattern
#' @examples data(schepCogapsResult)
#' data(schepGranges)
#'
#' GOenrichment <- findGOEnrichment(cogapsResult = schepCogapsResult, granges = schepGranges, genome = "hg19")
#' @export

findGOEnrichment <- function(cogapsResult, granges, genome, scoreThreshold = 0.05) {
  
  #get PatternMarker regions
  patMarkers = CoGAPS::patternMarkers(cogapsResult)
  patScores = as.data.frame(patMarkers$PatternMarkerScores)
  chr_regions = rownames(patScores)
  regionPatList = vector("list", length = ncol(patScores))
  for(i in seq(ncol(patScores))) {
    topPeaksPat = which(patScores[,i] < scoreThreshold)
    regionPatList[[i]] = topPeaksPat
  }
  
  #print number of peaks used based on patternMarker score threshold
  numPeaks = unlist(lapply(regionPatList, length))
  names(numPeaks) = lapply(seq(length(regionPatList)),
                           function(x) {paste("Pattern", x)})
  print("Number of peaks used for each pattern:", quote = FALSE)
  print(numPeaks)
  
  #call GREAT and get enrichment for GO terms
  GOList = vector("list", length(regionPatList))
  for(i in seq(length(regionPatList))) {
    patRanges = granges[regionPatList[[i]]]
  
    greatJob = rGREAT::submitGreatJob(gr = patRanges, species = genome,
                                      request_interval = 1)
    tbl = rGREAT::getEnrichmentTables(greatJob)
    GOBPaths = tbl$`GO Biological Process`
    GOBPaths = GOBPaths[order(GOBPaths$Binom_Adjp_BH),
                        ][, c("ID", "name", "Binom_Fold_Enrichment", "Binom_Adjp_BH")]
    GOList[[i]] = GOBPaths
  }
  
  return(GOList)
}
