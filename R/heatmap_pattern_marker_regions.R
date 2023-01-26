#' Create Heatmap of PatternMarker Peaks
#'
#' Function to make a heatmap of the accessibility of the most differentially
#' accessible regions as discovered by CoGAPS.
#'
#' @section Note: If you get the error:  "Error in plot.new() : figure margins
#'   too large" while using this function in RStudio just make the plotting pane
#'   in Rstudio larger and run the code again; this error only means the legend
#'   is being cut off in any case, the main plot will still appear correctly
#'
#' @param cgaps_result CogapsResult object from CoGAPS run
#' @param atac_data a numeric matrix of the ATAC data input to CoGAPS
#' @param celltypes a list or factor of celltypes corresponding to the positions
#'   of those cells in the atac_data matrix
#' @param numregions number of chromosomal regions/peaks to plot for each CoGAPS
#'   pattern. Default is 50. Plotting very large numbers of regions can cause
#'   significant slowdown in runtime
#' @param colColors column-wise colors for distinguishing celltypes. If NULL,
#'   will be generated randomly
#' @param rowColors row-wise colors for distinguishing patterns. If NULL will be
#'   generated randomly
#' @param patterns which patterns should be plotted, if NULL all will be plotted
#' @param order option whether to sort the data by celltype before plotting,
#'   TRUE by default
#' @param ... additional arguments to the heatmap.2 function
#' @return heatmap of the accessibility for numregions for each pattern
#' @examples data("schepCogapsResult")
#' data(schepCellTypes)
#' data("subsetSchepData")
#'
#' heatmapPatternMarkers(schepCogapsResult, atac_data = subsetSchepData,
#'                       celltypes = schepCellTypes, numregions = 50)
#' @export

heatmapPatternMarkers = function(cgaps_result, atac_data, celltypes, numregions = 50,
                                 colColors = NULL, rowColors = NULL, patterns = NULL,
                                 order = TRUE,...) {

  #convert atac read data to binary matrix
  binary_atac <- (atac_data > 0) + 0 #this line works because R represents TRUE/FALSE as 1/0

  #get regions corresponding to most elevated PatternMarker results
  patMarkers <- CoGAPS::patternMarkers(cgaps_result, threshold = "cut")
  patRanks <- as.data.frame(patMarkers[2])
  chr_regions <- rownames(patRanks)
  regionPatList <- vector(mode=  "list", length = ncol(patRanks))
  for(i in seq(ncol(patRanks))) {
    topPeaksPat <- order(patRanks[,i])[seq(numregions)]
    regionPatList[[i]] <- topPeaksPat
  }
  
  if(is.null(patterns)){
    patterns <- seq_along(regionPatList)
  }
  #subset the regions most differentially accessible for each pattern from the original data
  allPatSubset <- data.frame()
  for(i in patterns) {
    temppatsub <- as.matrix(binary_atac[unlist(regionPatList[i]),])
    allPatSubset <- rbind(allPatSubset, temppatsub)
  }

  allPatSubset <- as.matrix(allPatSubset)

  if(is.null(rowColors)){
    #produce vectors of colors for visualizing patterns in the heatmap
    rowColors <- rainbow(ncol(patRanks))
  }
  rowCols <- character(length(patterns)*numregions)
  j<-1
  for(i in patterns){
    color <- rep(rowColors[i], numregions)
    rowCols[j:(j+numregions-1)] <- color
    j<-j+numregions
  }

  celltypes <- as.factor(celltypes)

  if(order == TRUE) {
  #order the data by celltype
  allPatSubset <- rbind(allPatSubset, celltypes)
  ind <- (numregions*length(patterns))+1
  allPatSubset <- allPatSubset[,order(allPatSubset[ind,])]
  allPatSubset <- allPatSubset[-c(ind),]
  }

  if(is.null(colColors)){
  #produce vectors of colors for visualizing celltypes in the heatmap
  colsgen <- rainbow(length(levels(celltypes)))
  colsgen <- gtools::permute(colsgen)
  if(order == TRUE) {
  colColors <- colsgen[as.numeric(sort(celltypes))]
  }
  else{
    colColors <- colsgen[as.numeric(celltypes)]
  }
  }

  else{
    if(order == TRUE) {
      colColors <- colColors[as.numeric(sort(celltypes))]
    }
    else{
      colColors <- colColors[as.numeric(celltypes)]
    }
  }

  if(order == TRUE) {
    gplots::heatmap.2(allPatSubset, density.info="none", trace="none",
                      dendrogram='none', Rowv=FALSE, Colv=FALSE,
                      ColSideColors = colColors, RowSideColors = rowCols,
                      labRow = NA, labCol = as.character(sort(celltypes)), ...)
  }
  else {
    gplots::heatmap.2(allPatSubset, density.info="none", trace="none",
                      dendrogram='none', Rowv=FALSE, Colv=FALSE,
                      ColSideColors = colColors, RowSideColors = rowCols,
                      labRow = NA, labCol = as.character(celltypes), ...)
  }
}
