#' Map Peaks to DNA motifs in scATAC-seq Data
#'
#' Provides functionality to summarize scATAC-seq data by motifs from peak
#' summary. Uses motifmatchr to prepare data for CoGAPS run using motif
#' summarization
#' 
#' @param motifList PWMatrixList object of motifs (from the TFBS tools package)
#' @param scATACData matrix of scATACseq data, peaks (rows) by cells (columns)
#' @param granges GenomicRanges object corresponding to all peaks used to summarize scATACData
#' @param genome The UCSC Genome to use for input to motifmatchr (e.g "hg19")
#' @param cellNames List of cellnames corresponding to the cells in scATACData
#' @param pCutoff p-value cutoff for motifmatchr, 5e-09 by default to identify only matches with high confidence
#' @return matrix for input to CoGAPS with summary to motifs; motifs by cells
#' @examples \donttest{
#' motifSummTest = motifSummarization(motifList = motifs, scATACData = scatac, granges = peakGranges, genome = "hg19", cellNames = cells, pCutoff = 5e-09)
#' }
#' @export 

motifSummarization = function(motifList, scATACData, granges, genome,
                              cellNames, pCutoff = 5e-09) {
  #call motifmatchr to find DNA motifs in peak regions
  matchedmotifs = motifmatchr::matchMotifs(motifList, granges, 
                                           genome = genome, out = "scores",
                                           p.cutoff = pCutoff)
  region_motif_matches = motifmatchr::motifMatches(matchedmotifs)
  
  #finding motifs that are never matched
  unmatched=apply(region_motif_matches, 2, function(x) {sum(x)==0})
  unmatchedMotifs = which(unmatched == TRUE)
  
  #finding indices of motifmatches in the data
  mminds = apply(region_motif_matches, 2, function(x)
  {
    motifmatches = which(x ==TRUE)
    return(list(motifmatches))
  })
  
  #preallocate the matrix
  motif_scATACData_matrix = matrix(nrow = length(mminds),
                                   ncol = ncol(scATACData))
  #aggregate motif matches in individual cells to produce motif x cells matrix
  for(i in seq_along(mminds)) {
    for(j in seq(ncol(scATACData))){
      chrRegions = scATACData[mminds[[i]][[1]], j]
      readsum = sum(chrRegions)
      motif_scATACData_matrix[i, j] = readsum
    }
  }
  
  #name rows and columns
  colnames(motif_scATACData_matrix) = cellNames
  rownames(motif_scATACData_matrix) = TFBSTools::ID(motifList)
  
  #remove motifs that are never matched
  if(length(unmatchedMotifs) < 0) {
    motif_scATACData_matrix = motif_scATACData_matrix[-unmatchedMotifs,]
  }
  
  return(motif_scATACData_matrix)
}
