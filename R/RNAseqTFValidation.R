#' Validate TF Findings with RNA-seq CoGAPS
#'
#' Use results from CoGAPS run on matched RNA-seq data to verify TF activity
#' suggested by motif matching analysis of ATAC CoGAPS output. Uses the fgsea
#' package to find enrichment of PatternMarker genes among genes regulated by
#' identified candidate TFs
#'
#' @param TFGenes genes regulated by the TFs as returned by simpleMotifTFMatch()
#'   or findRegulatoryNetworks()
#' @param RNACoGAPSResult CoGAPSResult object from matched RNA-seq data, or, if
#'   matrix = TRUE, a matrix containing patternMarker gene ranks. Must contain
#'   gene names
#' @param ATACPatternSet vector of patterns found by CoGAPS in the ATAC data to
#'   match against patterns found in RNA
#' @param RNAPatternSet vector of patterns found by CoGAPS in RNA to match
#'   against those found in ATAC
#' @param matrix TRUE if inputting matrix of PatternMarker genes, FALSE if
#'   inputting CoGAPS result object. FALSE by default
#'
#' @return Result matrices from the fgsea function for each pattern comparison
#' @examples \dontrun{
#' gseaList = RNAseqTFValidation(TFMatchResult$RegulatoryNetworks, RNACoGAPS,
#'  c(1,3), c(2,7), matrix = FALSE)
#' }
#' @export

RNAseqTFValidation = function(TFGenes, RNACoGAPSResult, ATACPatternSet, 
                              RNAPatternSet, matrix = FALSE) {
  
  if(matrix == FALSE) {
  patMarkers = CoGAPS::patternMarkers(RNACoGAPSResult)
  genesRanks = patMarkers$PatternMarkerRanks
  }
  else{
    genesRanks = RNACoGAPSResult
  }
  
  gseaResults = vector(mode = "list",
                       length = length(ATACPatternSet)*length(RNAPatternSet))
  k=1
  for(i in ATACPatternSet){
    for(j in RNAPatternSet){
      gsea = suppressWarnings(fgsea::fgsea(TFGenes[[i]], genesRanks[,j], 50000))
      gseaResults[[k]] = gsea[order(gsea$pval),]
      print(head(gsea[order(gsea$pval), 1:3], n=10))
      k=k+1
    }
  }
  
  return(gseaResults)
}
