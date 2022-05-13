#' Small subset of the scATAC-seq data from Schep et al, 2017, Nature Methods paper. 
#' 
#' Subset from the Schep et al data, used for examples in this package. 
#' 
#' @format A matrix with 5036 peaks and 600 cells in the order of the schepPeaks,
#'  schepCellTypes, and schepGranges data objects.
#'  
#' @source \url{10.1038/nmeth.4401}
"subsetSchepData"

#' GRanges corresponding to subsetSchepData
#' 
#' GRanges in the order of the peaks of the subsetSchepData object from the
#' Schep et al, 2017, Nature Methods paper.
#' 
#' @format GRanges of length 5036
#' 
#' @source \url{10.1038/nmeth.4401}
"schepGranges"

#' Peaks corresponding to subsetSchepData
#' 
#' Character vector of peaks in the order of the peaks of the subsetSchepData object from the
#' Schep et al, 2017, Nature Methods paper.
#' 
#' @format Character vector of length 5036
#' 
#' @source \url{10.1038/nmeth.4401}
"schepPeaks"

#' Cell types corresponding to subsetSchepData
#' 
#' Factor of cell types in the order of the subsetSchepData object from the
#' Schep et al, 2017, Nature Methods paper.
#' 
#' @format Factor of length 600 with 12 levels
#' 
#' @source \url{10.1038/nmeth.4401}
"schepCellTypes"

#' CogapsResult from the subsetSchepData object
#' 
#' Output from applying the CoGAPS algorithm to the subsetSchepData object.
#' 
#' @format Large CogapsResult
"schepCogapsResult"

#' List of human TFs and motifs from cisBP database
#' 
#' Information on human TFs and their correpsonding DNA motifs
#' 
#' @format Data frame with 95413 rows and 28 columns.
#' 
#' @source \url{http://cisbp.ccbr.utoronto.ca/}
"tfData"


#' Example list of motifs for examples
#' 
#' PWMMatrixList used for examples with functions based on DNA motifs.
#' Each entry contains the motif ID and the probability of each nucleotide
#' at each position, as a matrix.
#' 
#' @format PWMMatrixList of length 100
"exampleMotifList"


