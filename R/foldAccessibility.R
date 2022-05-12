#' Estimate fold Accessibility of a Gene Relative to Average
#'
#' Compares the accessibility of peaks overlapping with a gene, as returned by
#' the geneAccessibility function to the average accessibility of peaks within a
#' given cell population. Meant to provide a rough estimate of how accessible a
#' gene is with values higher than 1 providing evidence of differential
#' accessibility (and thus implying possible transcription), with values lower
#' than 1 indicating the opposite.
#'
#' @param peaksAccessibility the binarized accessibility of a set of peaks; one
#'   value returned from the geneAccessibility function
#' @param cellTypeList list of celltypes grouping cells in the data
#' @param cellType the particular cell type of interest from within cellTypeList
#' @param binaryMatrix binarized scATAC data matrix
#' @return Fold accessibility value as compared to average peaks for a given
#'   cell type
#' @examples data("subsetSchepData")
#' data(schepCellTypes)
#' library(Homo.sapiens)
#' geneList <- c("TAL1", "IRF1")
#' data(schepGranges)
#' binarizedData <- (subsetSchepData > 0) + 0
#' accessiblePeaks <- geneAccessibility(geneList = geneList, peakGranges = schepGranges,
#'  atacData = subsetSchepData, genome = Homo.sapiens)
#' foldAccessibility(peaksAccessibility = accessiblePeaks$TAL1, cellTypeList = schepCellTypes,
#'  cellType = "K562 Erythroleukemia", binaryMatrix = binarizedData)
#' @export
foldAccessibility <- function(peaksAccessibility, cellTypeList, cellType, binaryMatrix) {
  
  if(is.list(peaksAccessibility)) {
    stop("Only one element of the list returned by the geneAccessibility function should be input as the peaksAccessibility parameter")
  }
  
  if(is.null(nrow(peaksAccessibility))) {
    percentAccessibility <- (sum(peaksAccessibility[which(cellTypeList == cellType)]))/sum(cellTypeList == cellType)
  }
  else {
    percentAccessibility <- (sum(apply(peaksAccessibility[,which(cellTypeList == cellType)], 1, sum))/nrow(peaksAccessibility))/sum(cellTypeList == cellType)
  }
  averageAccessibility <- mean(apply(binaryMatrix[,which(cellTypeList == cellType)],
                                     1, sum))/sum(cellTypeList == cellType)
  foldAcc <- percentAccessibility/averageAccessibility
  return(foldAcc)
}
