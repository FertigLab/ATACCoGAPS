#function to find motifs in each accessible region for a CoGAPS pattern
patternAccessibleMotifs = function(pattern, generanges, motifs, genome,
                                   motifsPerRegion = 1) {

  #get the GRanges for the differentially accessible regions
  patRanges = generanges[pattern]

  #use motifmatchr to match motifs to regions
  matchedmotifs = motifmatchr::matchMotifs(motifs, patRanges, genome = genome,
                                           out = "scores")
  motif.scores = motifmatchr::motifScores(matchedmotifs)
  motifIDs = TFBSTools::ID(motifs)

  #find the mostly closely matching motif for each region
  topMotifs = vector(mode = "list", nrow(motif.scores))
  for(i in seq(nrow(motif.scores))) {
    topScores = sort(motif.scores[i,], decreasing = TRUE)[seq(motifsPerRegion)]
    if(topScores==0){
      topMotifs[i] = "NoMatchedMotif_NoMatchedTF"
      names(topMotifs)[i] = "NoMatchedMotif_NoMatchedTF"
    }
    else{
    inds = which(motif.scores[i,] %in% topScores)
    tmpmotifs = motifIDs[inds]
    topMotifs[i] = list(tmpmotifs)
    }
  }

  return(topMotifs)
}


#' Find Motifs and TFs from PatternMarker Peaks
#'
#' Function that takes CoGAPS result and list of DNA motifs as input and returns
#' motifs which match to the most pattern-defining peaks for each pattern.
#'
#' @param cogapsResult the result object from a CoGAPS run
#' @param generanges GRanges objects corresponding to the genomic regions which
#'   form the rows of the ATAC-seq data that CoGAPS was run on
#' @param motiflist a PWMlist of motifs to search the regions for
#' @param genome the ucsc genome version to use e.g. "hg19", "mm10"
#' @param scoreThreshold threshold for the most pattern defining peaks as per
#'   the PatternMarker statistic from the CoGAPS package. By default is NULL, in
#'   which case all Pattern defining peaks will be used for motif matching. Used
#'   to reduce compute time, as results are quite robust across thresholds
#' @param motifsPerRegion number of top motifs to return from each peak
#' @return motifPatternMatch: nested list of the top motif for each region for x
#'   number of regions for each pattern
#' @examples data(exampleMotifList)
#' data(schepGranges)
#' data(schepCogapsResult)
#' 
#' motifsByPattern = motifPatternMatch(schepCogapsResult, schepGranges, exampleMotifList, "hg19")
#' @export
motifPatternMatch = function(cogapsResult, generanges, motiflist,
                             genome, scoreThreshold = NULL, motifsPerRegion = 1) {

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
    regionPatList = vector(mode = "list", length = ncol(patScores))
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

  #run patternAccessibleMotifs for all patterns
  filenames = vector("list", length(regionPatList))
  for(i in seq_along(regionPatList)) {
    motifstmp = suppressMessages(suppressWarnings(patternAccessibleMotifs(regionPatList[[i]],
                                                                          generanges, motiflist,
                                                                          genome, motifsPerRegion)))
    nam <- paste("pattern", i, "motifs", sep = "")
    filenames[i] = nam
    assign(nam, motifstmp)
  }

  #return nested list of motifs by pattern
  ind =paste(filenames, collapse = ",")
  motifslist = eval(parse(text = paste("list(", ind, ")")))
  return(motifslist)

}

#' @describeIn motifPatternMatch Match motifs to TFs based on the list of motifs
#'   returned by motifPatternMatch
#'
#' @param motifList list produced by the motifPatternMatch function
#' @param tfData dataframe of motifs and TFs from cisBP database
#' @return getTFs: list containing list of dataframes of tfData subset to matched TFs
#'   and list of how many times each TF was matched to a motif/peak
#' @examples data(exampleMotifList)
#' data(schepGranges)
#' data(schepCogapsResult)
#' data(tfData)
#' 
#' motifsByPattern = motifPatternMatch(schepCogapsResult, schepGranges, exampleMotifList, "hg19")
#' motifTFs = getTFs(motifsByPattern, tfData)
#' @export
getTFs = function(motifList, tfData) {
  #lists to iteratively record TF info for each pattern
  matchInfoList = vector("list", length(motifList))
  summaryList = vector("list", length(motifList))
  #iterate over motifs for each pattern and find matches in tfData
  for(i in seq_along(motifList)) {
    match_info = tfData[which(tfData$Motif_ID %in% unlist(motifList[[i]])),]
    TfNames = match_info$TF_Name
    TF.fct = as.factor(TfNames)
    #summary of the family that the TFs belong to
    TF_summary = summary(TF.fct)
    matchInfoList[[i]] = match_info
    summaryList[[i]] = TF_summary
  }

  summaryList = lapply(summaryList, sort, decreasing = TRUE)

  return(list(matchInfoList, summaryList))
}

#' @describeIn motifPatternMatch function to match TFs identified by getTFs
#'   function to a list of regulatory networks of genes known for those TFs
#'
#' @param TFs result from getTFs function
#' @param networks a list of regulatory networks of genes corresponding to TFs;
#'   we include humanRegNets and mouseRegNets, downloaded from the TTrust
#'   database (Han et al Nucleic Acid Res. 2018)
#' @return findRegulatoryNetworks: list of TFs for which we have annotations and
#'   the corresponding gene networks for each pattern
#' @examples data(exampleMotifList)
#' data(schepGranges)
#' data(schepCogapsResult)
#' data(tfData)
#' 
#' motifsByPattern = motifPatternMatch(schepCogapsResult, schepGranges, exampleMotifList, "hg19")
#' motifTFs = getTFs(motifsByPattern, tfData)
#' 
#' regNets = findRegulatoryNetworks(motifTFs, ATACCoGAPS:::humanRegNets)
#' @export
findRegulatoryNetworks = function(TFs, networks) {
  tfOut = TFs[[1]]

  regNetList = vector("list", length(tfOut))
  for(i in seq_along(tfOut)) {
    patList=networks[which(names(networks) %in% tfOut[[i]]$TF_Name)]
    regNetList[[i]] = patList
  }
  return(regNetList)
}

#' @describeIn motifPatternMatch function to match functional annotation to a
#'   list of TFs from the getTFs function
#'
#' @param TFs object of TF info returned from the getTFs function
#' @return getTFDescriptions: list of functional annotations for all TFs in each
#'   pattern
#' @examples data(exampleMotifList)
#' data(schepGranges)
#' data(schepCogapsResult)
#' data(tfData)
#' 
#' motifsByPattern = motifPatternMatch(schepCogapsResult, schepGranges, exampleMotifList, "hg19")
#' motifTFs = getTFs(motifsByPattern, tfData)
#' 
#' tfDesc = getTFDescriptions(motifTFs)
#' @export
getTFDescriptions = function(TFs) {
  tfOut = TFs[[1]]

  geneDescriptionList = vector("list", length(tfOut))
  for(i in seq_along(tfOut)) {
    patAnn = humanGeneAnnotations[which(humanGeneAnnotations$Symbol %in% tfOut[[i]]$TF_Name),
                                  c("summary", "Symbol")]
    patAnn =cbind(patAnn$Symbol, patAnn$summary)
    geneDescriptionList[[i]] = patAnn
  }
  return(geneDescriptionList)
}

