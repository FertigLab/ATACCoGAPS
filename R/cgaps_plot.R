#' Plot Individual CoGAPS Patterns
#'
#' Function to plot each pattern of the pattern matrix from a cogapsResult and
#' color by cell classifier information to identify which patterns identify
#' which cell classes.
#'
#' @param cgaps_result CoGAPSResult object from a CoGAPS run or the pattern
#'   matrix (matrix must be set equal to TRUE in the latter case)
#' @param sample.classifier factor of sample classifications for all cells for
#'   the data to be plotted by (e.g. celltypes)
#' @param cols  vector of colors to be used for the cell classes; should have
#'   the same number of colors as levels of the sample.classifier factor. If
#'   left null a list of colors is produced
#' @param sort TRUE if samples will be sorted according to sample.classifier
#'   prior to plotting
#' @param patterns numerical vector of patterns to be plotted; if null all
#'   patterns are plotted
#' @param matrix if false cgaps_result is interpreted as a CoGAPSResult object,
#'   if true it is interpreted as the pattern matrix
#' @param ... addition arguments to plot function
#' @return Series of plots of pattern matrix patterns colored by cell
#'   classifications
#' @examples data("schepCogapsResult")
#' data(schepCellTypes)
#'
#' cgapsPlot(schepCogapsResult, schepCellTypes)
#' @export

cgapsPlot <- function(cgaps_result, sample.classifier, cols = NULL, sort = TRUE,
                     patterns = NULL, matrix = FALSE, ...){

  #if inputting the pattern matrix rather than full cogaps object
  if(matrix == TRUE) {
    samplePatterns <- cgaps_result
  }
  else{
  #get the pattern matrix from the cogaps result
  samplePatterns <- CoGAPS::getSampleFactors(cgaps_result)
  }
  #get the number of cell classifications to be plotted
  cellnum <- length(levels(sample.classifier))



  #if no color vector is input generate a random list of colors to plot by
  if(is.null(cols)){
    cols <- rainbow(cellnum)
    cols <- gtools::permute(cols)
  }
  #if the desired patterns to plot is not input, find the total number of patterns
  #to plot all patterns
  if(is.null(patterns)){
    nPatterns <- ncol(samplePatterns)
    patterns <- seq_len(nPatterns)
  }

  #sort all samples so samples of the same factor level are adjacent in the plot
  #if sort == TRUE
  if(sort == TRUE) {
    #set the graphical parameters
    old <- options(stringsAsFactors = FALSE)
    on.exit(options(old), add = TRUE)
    par(oma=c(0, 0, 0, 5))

    for(i in patterns) {
      pattern <- samplePatterns[,i]
      sorted_pat <- pattern[order(sample.classifier)]

    #loop through each column (pattern) in the pattern matrix
      plotname <- paste("Pattern", i)

      #plot the pattern, with legend outside of the graphing region
      plot(sorted_pat, col =cols[as.numeric(sample.classifier[order(sample.classifier)])],
           main = plotname, ...)
      legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
             legend = levels(sample.classifier), col= cols, cex = 0.6, pch = 16)
    }
  }

  else {
    #set the graphical parameters
    old <- options(stringsAsFactors = FALSE)
    on.exit(options(old), add = TRUE)
    par(oma=c(0, 0, 0, 5))

    #loop through each column (pattern) in the pattern matrix
    for(i in patterns) {
      pattern <- samplePatterns[,i]
      plotname <- paste("Pattern", i)
      #plot the pattern, with legend outside of the graphing region
      plot(pattern, col =cols[as.numeric(sample.classifier)], main = plotname, ...)
      legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
             legend = levels(sample.classifier), col= cols, cex = 0.6, pch = 16)
    }
  }
}
