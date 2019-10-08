#' Transfer Learning between ATACseq data sets using projectR
#'
#' Wrapper function for projectR which finds overlaps between the peaks of the
#' atac data CoGAPS was run on and maps them to new data set the user wishes to
#' project learned patterns into.
#'
#' @param newData the ATAC data to project into
#' @param CoGAPSResult result from CoGAPS run on original ATAC data
#' @param originalPeaks peaks from the ATAC data Cogaps was run on
#' @param originalGranges granges of the peaks for the data set Cogaps was run
#'   on
#' @param newGranges granges of the peaks for the new data set
#' @return A matrix of the projected patterns in the input data as well as
#'   p-values for each element of that matrix.
#' @export


ATACTransferLearning = function(newData, CoGAPSResult, originalPeaks,
                                originalGranges, newGranges) {

  olaps = as.matrix(GenomicRanges::findOverlaps(originalGranges, newGranges))
  olap_peaks=originalPeaks[olaps[,1]]
  newDataSub = newData[olaps[,2],]

  pR = projectR::projectR(newDataSub, CoGAPSResult, olap_peaks, originalPeaks,
                          full = T)

  return(pR)
}
