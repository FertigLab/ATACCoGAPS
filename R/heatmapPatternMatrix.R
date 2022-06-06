#' Plot the patternMatrix as a heatmap
#'
#' Selects the pattermMatrix (patterns by cells) from the CoGAPSResult and plots
#' the data as a heatmap. Intended to visualize the celltypes distinguished by
#' the patterns found by CoGAPS.
#'
#' @param cgaps_result CoGAPSResult object from a CoGAPS run or the pattern
#'   matrix (matrix must be set equal to TRUE in the latter case)
#' @param patterns numerical vector of patterns to be plotted; if null all
#'   patterns are plotted
#' @param sample.classifier factor of sample classifications for all cells for
#'   the data to be plotted by (e.g. celltypes)
#' @param cellCols vector of colors to be used for the cell classes; should have
#'   the same number of colors as levels of the sample.classifier factor. If
#'   left null a list of colors is produced
#' @param sort TRUE if samples will be sorted according to sample.classifier
#'   prior to plotting
#' @param matrix if false cgaps_result is interpreted as a CoGAPSResult object,
#'   if true it is interpreted as the pattern matrix being input directly
#' @param rowColors vector of colors to plot along patterns, if NULL generated automatically
#' @param ... additional arguments to the heatmap.2 function
#' @return Heatmap of patternMatrix with color labels for samples
#' @examples data("schepCogapsResult")
#' data(schepCellTypes)
#'
#' heatmapPatternMatrix(schepCogapsResult, sample.classifier = schepCellTypes)
#' @export
heatmapPatternMatrix <- function(cgaps_result, sample.classifier, cellCols = NULL,
                                sort = TRUE, patterns = NULL, matrix = FALSE,
                                rowColors = NULL, ...) {

  #if inputting the pattern matrix rather than full cogaps object
  if(matrix == TRUE) {
    samplePatterns <- cgaps_result
  }
  else{
    #get the pattern matrix from the cogaps result
    samplePatterns <- CoGAPS::getSampleFactors(cgaps_result)
  }


  #if no color vector is input generate a random list of colors to plot by
  if(is.null(cellCols)){
    #get the number of cell classifications to be plotted
    cellnum <- length(levels(sample.classifier))
    cellCols <- rainbow(cellnum)
    cellCols <- gtools::permute(cellCols)
  }

  #record if a specific number of patterns to plot is input
  if(!is.null(patterns)){patternsGiven <- TRUE}
  else{patternsGiven <- FALSE}

  #if the desired patterns to plot is not input, find the total number of patterns
  #to plot all patterns
  if(is.null(patterns)){
    nPatterns <- ncol(samplePatterns)
    patterns <- seq(nPatterns)
  }
  
  if(is.null(rowColors)) {
    #produce vectors of colors for visualizing patterns in the heatmap
    rowColors <- rainbow(length(patterns))
  }

  #sort all samples so samples of the same factor are adjacent in the plot
  #if sort == TRUE
  if(sort == TRUE) {
    samplePatterns <- apply(samplePatterns, 2, function(x){x[order(sample.classifier)]})

    #plot heatmap
    if(patternsGiven == TRUE) {samplePatterns<- samplePatterns[,patterns]}
    gplots::heatmap.2(t(samplePatterns), density.info="none", trace="none",
                      dendrogram='none', Rowv=NULL, Colv=NULL,
                      RowSideColors = rowColors,
                      ColSideColors = cellCols[as.numeric(sample.classifier[order(sample.classifier)])],
                      labRow = NA, labCol = as.character(sort(sample.classifier)), ...)
  }

  else{
    #plot heatmap
    if(patternsGiven == TRUE) {samplePatterns<- samplePatterns[,patterns]}
    gplots::heatmap.2(t(samplePatterns), density.info="none", trace="none",
                      dendrogram='none', Rowv=NULL, Colv=NULL,
                      RowSideColors = rowColors,
                      ColSideColors = cellCols[as.numeric(sample.classifier)],
                      labRow = NA, labCol = sort(sample.classifier), ...)

  }

}

