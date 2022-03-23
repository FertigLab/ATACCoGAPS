#' Find the accessibility of the peaks overlapping a set of genes and their
#' promoters
#'
#' The accessibility of a particular set of interest genes is checked by testing
#' overlap of peaks with the genes and gene promoters and then returning the
#' binarized accesibility data for those peaks
#'
#' @param geneList vector of HGNC gene symbols to find overlapping peaks for in
#'   the data
#' @param peakGranges a GRanges object corresponding to the peaks in the
#'   atacData matrix, in the same order as the rows of the atacData matrix
#' @param atacData a single-cell ATAC-seq count matrix peaks by cells
#' @param genome TxDb object to produce gene GRanges from
#' @return List of matrices corresponding to the accessible peaks overlapping
#'   with each gene across all cells in the data
#' @examples library(Homo.sapiens)
#' geneList <- c("TAL1", "IRF1")
#' data(schepGranges)
#' data("subsetSchepData")
#' accessiblePeaks <- geneAccessibility(geneList = geneList, peakGranges = schepGranges, atacData = subsetSchepData, genome = Homo.sapiens)
#' @export
geneAccessibility <- function(geneList, peakGranges, atacData, genome) {
  genes <- geneRanges(genome) #uses geneRanges function from genePeakMatching.R file
  symbols <- mcols(genes)[["SYMBOL"]]
  #get promoters
  prms = GenomicFeatures::promoters(GenomicFeatures::genes(genome, columns = "SYMBOL"),
                                    upstream = 1500, downstream = 500)
  
  geneInds <- which(symbols %in% geneList)
  
  if(length(geneList) != length(geneInds)) {
    warning(paste("Some genes in list did not have matching symbol in input database.", setdiff(geneList, symbols[geneInds]), "is not found"))
  }
  
  listedGenes <- genes[geneInds]
  
  #promoters
  prmsymbols = unlist(mcols(prms)[["SYMBOL"]])
  prmInds = which(prmsymbols %in% geneList)
  listedPromoters <- prms[prmInds]
  
  olaps <- vector("list", length(listedGenes))
  toRemove <- NULL
  for(i in seq_along(listedGenes)) {
    olap1 <- GenomicRanges::findOverlaps(listedGenes[i], peakGranges,
                                        ignore.strand = TRUE)
    olapInds1 <- subjectHits(olap1)
    olap2 <- GenomicRanges::findOverlaps(listedPromoters[i], peakGranges,
                                         ignore.strand = TRUE)
    olapInds2 <- subjectHits(olap2)
    
    if(length(olapInds1) == 0 & length(olapInds2) == 0) {
      print(paste("No overlapping peaks with", mcols(listedGenes)[["SYMBOL"]][i]))
      toRemove <- c(toRemove, i)
    }
    else{
      olaps[[i]] <- c(olapInds1,olapInds2)
    }
  }
  
  if(!is.null(toRemove)){
    olaps[[toRemove]] <- NULL
  }
  
  binaryATAC <- (atacData > 0) + 0
  
  geneAccessibile <- vector("list", length(listedGenes))
  geneAccessibile <- lapply(olaps, function(x){
    binaryATAC[x,]
  })
  
  if(!is.null(toRemove)){
    names(geneAccessibile) <- symbols[geneInds][-toRemove]
  }
  else{
    names(geneAccessibile) <- symbols[geneInds]
  }
  
  return(geneAccessibile)
  
}

