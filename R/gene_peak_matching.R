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
  #make list of genes found, with format allowing for non-matches as recorded entries
  f1 = factor(subjectHits(olaps), levels=seq_along(subjectLength(olaps)))
  IRanges::splitAsList(mcols(query)[["SYMBOL"]][queryHits(olaps)], f1)
}

#function to find nearest gene not within region
nearestGenes = function(genes, patternranges){
  #find nearest genes to Granges objects from most pattern differentiating peaks
  nearid = GenomicRanges::nearest(patternranges, genes, ignore.strand = T)
  neargrs = genes[nearid]
  neargenes = mcols(neargrs)[["SYMBOL"]]
  return(neargenes)
}

#function to find nearest gene following the region
followingGenes = function(genes, patternranges){
  #find following genes to Granges objects from most pattern differentiating peaks
  followid = GenomicRanges::follow(patternranges, genes, ignore.strand = T)
  #check for NA values and only record cases with following genes available
  if(sum(is.na(followid))==0) {
    followgrs = genes[followid]
    followgenes = mcols(followgrs)[["SYMBOL"]]
  }
  else {
    naInds = which(is.na(followid))
    followid = followid[-naInds]
    followgrs = genes[followid]
    followgenes = mcols(followgrs)[["SYMBOL"]]
  }
  return(followgenes)
}


#function that matches genes for one pattern
geneMatch = function(regionIndex, generanges, genome) {

  #get the specific top ranges needed and concatenate into one GRanges object
  patRanges = generanges[regionIndex]

  gns = geneRanges(genome)
  genesInPatTmp = findOverlap(gns, patRanges)

  #convert output of findOverlap into a normal list
  #add names to each element for the patternMarker peak it's within
  genesInPat = c()
  for(i in seq_along(genesInPatTmp)){
    tmp = unlist(genesInPatTmp[i], use.names = F)
    if(length(tmp)==0) {
      genesInPat = c(genesInPat,NA)
      namelength = length(genesInPat)
      names(genesInPat)[namelength] = i
    }
    else if (length(tmp) > 1) {
      genesInPat = c(genesInPat,tmp)
      namelength = length(genesInPat)
      names(genesInPat)[(namelength-length(tmp)+1):namelength] = i
    }
    else{
      genesInPat = c(genesInPat,tmp)
      namelength = length(genesInPat)
      names(genesInPat)[namelength] = i
    }
  }

  #get nearest genes and genes following the peaks
  genesNearPat = nearestGenes(gns, patRanges)
  genesFollowPat = followingGenes(gns, patRanges)

  #return a list of lists with the three sets of gene information
  return(list(genesWithinRegion = genesInPat, genesNearRegion = genesNearPat,
              genesFollowingRegion = genesFollowPat))

}



#' Match genes to pattern differentiating peaks
#'
#' Function to take as input CoGAPS results for ATAC-seq data and find genes
#' within the most "pattern-defining" regions (as identified by cut thresholded
#' pattern Marker statistic from the CoGAPS package), as well as the nearest
#' gene and the nearest gene following the region
#'
#' @param cogapsResult the CogapsResult object produced by a CoGAPS run
#' @param generanges GRanges object corresponding to the genomic regions
#'   identified as peaks for the ATAC-seq data that CoGAPS was run on
#' @param genome Homo.sapiens or Mus.musculus currently supported
#' @param scoreThreshold threshold for the most pattern defining peaks as per
#'   the PatternMarker statistic from the CoGAPS package
#' @return double nested list containing lists of the genes in, nearest, and
#'   following the peaks matched each pattern
#' @examples data(schepCogapsResult)
#' data(schepGranges)
#'
#' genes = genePatternMatch(cogapsResult = schepCogapsResult, generanges = schepGranges, genome = Homo.sapiens, scoreThreshold = 0.03)
#' @export
genePatternMatch = function(cogapsResult, generanges, genome, scoreThreshold = 0.03) {

  #get PatternMarker peak indices
  patMarkers = CoGAPS::patternMarkers(cogapsResult)
  patScores = as.data.frame(patMarkers$PatternMarkerScores)
  chr_regions = rownames(patScores)
  regionPatList = vector(mode=  "list", length = ncol(patScores))
  for(i in seq(ncol(patScores))) {
    topPeaksPat = which(patScores[,i] < scoreThreshold)
    regionPatList[[i]] = topPeaksPat
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
