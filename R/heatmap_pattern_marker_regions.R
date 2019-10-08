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
#' @param numregions number of chromosomal regions to plot for each CoGAPS
#'   pattern
#' @param colColors column-wise colors for distinguishing celltypes. If NULL,
#'   will be generated randomly
#' @param seed random seed for generating colors for plot; for reproducibility
#' @param order option whether to sort the data by celltype before plotting,
#'   TRUE by default
#' @param ... additional arguments to the heatmap.2 function
#' @return heatmap of the accessibility for numregions for each pattern
#' @examples data(schepCogapsResult)
#' data(schepCelltypes)
#' data(schepFilteredData)
#'
#' heatmapPatternMarkers(schepCogapsResult, atac_data = schepFilteredData,
#'                       celltypes = schepCelltypes, numregions = 50)
#' @export

heatmapPatternMarkers = function(cgaps_result, atac_data, celltypes, numregions,
                                 colColors = NULL, seed = 42, order = TRUE,...) {

  #convert atac read data to binary matrix
  binary_atac = (atac_data > 0) + 0 #this line works because R represents TRUE/FALSE as 1/0

  #get regions corresponding to most elevated PatternMarker results
  patMarkers = CoGAPS::patternMarkers(cgaps_result,  threshold = "cut")
  patRanks = as.data.frame(patMarkers[2])
  chr_regions = rownames(patRanks)
  regionPatList = vector(mode=  "list", length = ncol(patRanks))
  for(i in seq(ncol(patRanks))) {
    topPeaksPat = order(patRanks[,i])[seq(numregions)]
    regionPatList[[i]] = topPeaksPat
  }

  #subset the regions most differentially accessible for each pattern from the original data
  allPatSubset = data.frame()
  for(i in seq_along(regionPatList)) {
    temppatsub = as.matrix(binary_atac[unlist(regionPatList[i]),])
    allPatSubset = rbind(allPatSubset, temppatsub)
  }

  allPatSubset = as.matrix(allPatSubset)


  #produce vectors of colors for visualizing patterns in the heatmap
  rainCols = rainbow(ncol(patRanks))
  rowColors = character(ncol(patRanks)*numregions)
  j=1
  for(i in seq(ncol(patRanks))){
    color = rep(rainCols[i], numregions)
    rowColors[j:(j+numregions-1)] = color
    j=j+numregions
  }

  celltypes = as.factor(celltypes)

  if(order == TRUE) {
  #order the data by celltype
  allPatSubset = rbind(allPatSubset, celltypes)
  ind = (numregions*ncol(patRanks))+1
  allPatSubset = allPatSubset[,order(allPatSubset[ind,])]
  allPatSubset = allPatSubset[-c(ind),]
  }

  if(is.null(colColors)){
  #produce vectors of colors for visualizing celltypes in the heatmap
  set.seed(seed)
  colsgen = rainbow(length(levels(celltypes)))
  colsgen = gtools::permute(colsgen)
  if(order == TRUE) {
  colColors = colsgen[as.numeric(sort(celltypes))]
  }
  else{
    colColors = colsgen[as.numeric(celltypes)]
  }
  }

  else{
    if(order == TRUE) {
      colColors = colColors[as.numeric(sort(celltypes))]
    }
    else{
      colColors = colColors[as.numeric(celltypes)]
    }
  }

  #plot heatmap using gplots package
  gplots::heatmap.2(allPatSubset, density.info="none", trace="none",
                    dendrogram='none', Rowv=FALSE, Colv=FALSE,
                    RowSideColors = rowColors, ColSideColors = colColors,
                    labRow = NA, labCol = NA, ...)
}
