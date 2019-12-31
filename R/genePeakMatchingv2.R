#function to produce the gene objects for comparison to the genomic ranges
#from the data
#txdb = organism TxDb with gene symbols
geneRanges = function(txdb) {
  #get granges with gene symmbol metadata
  g = GenomicFeatures::genes(txdb, columns="SYMBOL")
  col = mcols(g)[["SYMBOL"]]
  #create GRanges without metadata to make gene symbols character vector
  genes = GenomicRanges::granges(g)[rep(seq_along(g), elementNROWS(col))]
  mcols(genes)[["SYMBOL"]] = as.character(unlist(col))
  genes
}

#function to find genes within the top genomic regions
findOverlap = function(query, subject) {
  #find the overlaps between the GRanges containing peak regions and the GRanges
  #of known genes
  olaps = GenomicRanges::findOverlaps(query, subject, ignore.strand = TRUE)
  #get symbols
  symbols = mcols(query)[["SYMBOL"]][queryHits(olaps)]

  return(symbols)
}

#function that matches genes for one pattern
geneMatch = function(regionIndex, generanges, genome) {
  
  #get the specific top ranges needed and concatenate into one GRanges object
  patRanges = generanges[regionIndex]

  #get gene ranges for all genes in the genome
  gns = geneRanges(genome)
  
  #get overlaps with genes
  patGenes = findOverlap(gns, patRanges)
  
  #get promoters
  prms = GenomicFeatures::promoters(GenomicFeatures::genes(genome, columns = "SYMBOL"),
            upstream = 1500, downstream = 500)
  
  #overlaps with promoters
  patPromoters = findOverlap(prms, patRanges)
  patPromoters = unlist(patPromoters)
  
  #combine lists of genes with exonic or intronic overlap with a PatternMArker peak
  #to those with promoter overlap with a peak
  patOvGenes = c(patGenes, patPromoters)
  #filter duplicate genes
  patOvGenes = patOvGenes[-which(duplicated(patOvGenes))]
  return(patOvGenes)
}


#' Match genes to pattern differentiating peaks
#'
#' Function to take as input CoGAPS results for ATAC-seq data and find genes
#' within the most "pattern-defining" regions (as identified by cut thresholded
#' pattern Marker statistic from the CoGAPS package), as well as the nearest
#' gene and the nearest gene following the region. Note: a TxDb object for the
#' genome of interest must be loaded prior to running this function.
#'
#' @param cogapsResult the CogapsResult object produced by a CoGAPS run
#' @param generanges GRanges object corresponding to the genomic regions
#'   identified as peaks for the ATAC-seq data that CoGAPS was run on
#' @param genome A TxDb object for the genome of interest, it must be loaded
#'   prior to calling this function
#' @param scoreThreshold threshold for the most pattern defining peaks as per
#'   the PatternMarker statistic from the CoGAPS package. Default is NULL,
#'   returning all PatternMarker peaks. Useful to reduce computational time, as
#'   top results are reasonably robust to using more stringent thresholds
#' @return double nested list containing lists of the genes in, nearest, and
#'   following the peaks matched each pattern
#' @examples data(schepCogapsResult)
#' data(schepGranges)
#'
#' genes = genePatternMatch(cogapsResult = schepCogapsResult, generanges = schepGranges, genome = Homo.sapiens)
#' @export
genePatternMatch = function(cogapsResult, generanges, genome, scoreThreshold = NULL) {
  
  #get PatternMarker peak indices
  patMarkers = CoGAPS::patternMarkers(cogapsResult)
  if(is.null(scoreThreshold)){
    peaks = patMarkers$PatternMarkers
    nPeaks = lapply(peaks, length)
    PMRanks = patMarkers$PatternMarkerRanks
    regionPatList = vector(mode = "list", length = ncol(PMRanks))
    for(i in seq(ncol(PMRanks))) {
      topPeaksPat = which(PMRanks[,i] %in% seq(nPeaks[[i]]))
      regionPatList[[i]] = topPeaksPat
    }
  }
  else{
    patScores = as.data.frame(patMarkers$PatternMarkerScores)
    chr_regions = rownames(patScores)
    regionPatList = vector(mode=  "list", length = ncol(patScores))
    for(i in seq(ncol(patScores))) {
      topPeaksPat = which(patScores[,i] < scoreThreshold)
      regionPatList[[i]] = topPeaksPat
    }
  }
  
  #print number of peaks used based on patternMarker score threshold
  numPeaks = unlist(lapply(regionPatList, length))
  names(numPeaks) = lapply(seq(length(regionPatList)),
                           function(x) {paste("Pattern", x)})
  print("Number of peaks used for each pattern:", quote = FALSE)
  print(numPeaks)
  
  #run the geneMatch function for every pattern
  filenames = vector(mode = "list", length = length(regionPatList))
  for(i in seq_along(regionPatList)) {
    patgenetmp = suppressMessages(suppressWarnings(geneMatch(regionPatList[[i]],
                                                             generanges,
                                                             genome = genome)))
    nam <- paste("pattern", i, "genes", sep = "")
    filenames[i] = nam
    assign(nam, patgenetmp)
  }
  
  #put all patterns into a double nested list to be returned as output
  ind =paste(filenames, collapse = ",")
  geneslist = eval(parse(text = paste("list(", ind, ")")))
  
  return(geneslist)
}
