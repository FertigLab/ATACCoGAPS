#' Heatmap Gene Accessibility
#'
#' Use the output from geneAccessibility function to plot a heatmap of the
#' accessible peaks for a particular gene.
#'
#' @param genePeaks The peaks corresponding to a singular gene; one element of
#'   the list ouptut by geneAccessibility()
#' @param celltypes List or factor of celltypes corresponding to the cells in
#'   the scATAC-seq data set the peaks were found in
#' @param colColors A vector of colors to color the celltypes by, if NULL a
#'   random vector of colors is generated
#' @param order should the data be ordered by the celltype classifier? TRUE by
#'   default
#' @param seed random seed to generate colors if colColors is NULL
#' @param ... additional arguments to the heatmap.2 function from the gplots
#'   package
#' @return A plot of the peaks overlapping with a particular gene of interest
#' @examples library(Homo.sapiens)
#' geneList <- c("TAL1", "EGR1")
#' data(schepGranges)
#' data(schepFilteredData)
#' data(schepCelltypes)
#' accessiblePeaks <- geneAccessibility(geneList = geneList, peakGranges = schepGranges, atacData = schepFilteredData, genome = Homo.sapiens)
#' heatmapGeneAccessibility(genePeaks = accessiblePeaks$TAL1, celltypes = schepCelltypes)
#' @export
heatmapGeneAccessibility <- function(genePeaks, celltypes, colColors = NULL, order = TRUE, seed = 42, ...) {
  
  if(is.list(genePeaks)) {
    stop("Only one element of the list returned by the geneAccessibility function should be input as the genePeaks parameter")
  }
  
  celltypes = as.factor(celltypes)
  
  if(order == TRUE) {
    #order the data by celltype
    genePeaks = rbind(genePeaks, celltypes)
    ind = nrow(genePeaks)
    genePeaks = genePeaks[,order(genePeaks[ind,])]
    genePeaks = genePeaks[-c(ind),]
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
  
  if(is.null(nrow(genePeaks))) {
    warning("There is only one peak matched to this gene, adding a row of zeroes so the peak can be plotted using a heatmap.")
    zeroRow <- rep(0, length(genePeaks))
    genePeaks = rbind(genePeaks, zeroRow)
  }
  
  if(order == TRUE) {
    gplots::heatmap.2(genePeaks, density.info="none", trace="none",
                      dendrogram='none', Rowv=FALSE, Colv=FALSE, ColSideColors = colColors,
                      labRow = NA, labCol = as.character(sort(celltypes)), ...)
  }
  else {
    gplots::heatmap.2(genePeaks, density.info="none", trace="none",
                      dendrogram='none', Rowv=FALSE, Colv=FALSE, ColSideColors = colColors,
                      labRow = NA, labCol = as.character(celltypes), ...)
  }
}
