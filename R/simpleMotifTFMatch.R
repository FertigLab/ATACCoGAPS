#' Motif/TF Matching in a Single Function
#'
#' If the user does not have a specific set of motifs, transcription factors, or
#' regulatory networks that they want to match against, simply uses the core
#' motifs from the JASPAR database to find motifs and TFs in the most Pattern
#' differentiating peaks, as well as regulatory networks from TTrust database
#' corresponding to the identified TFs. This is used to provide transcription
#' factors with functional annotation which may suggest plausible unknown
#' regulatory mechanisms operating in the cell types of interest within the
#' data.
#'
#' @param cogapsResult result object from CoGAPS run
#' @param generanges GRanges object corresponding to peaks in ATACseq data
#'   CoGAPS was run on
#' @param organism organism name (e.g. "Homo sapiens")
#' @param genome genome version to use (e.g. hg19, mm10)
#' @param scoreThreshold threshold for the most pattern defining peaks as per
#'   the PatternMarker statistic from the CoGAPS package
#' @param motifsPerRegion number of motifs to attempt to find within each peak
#' @return list containing list of matched motifs, list of transciption factors,
#'   regulatory gene networks known for those TFs, functional annotations,
#'   summary showing how many times each TF was matched to a peak, and the
#'   downloaded set of motifs for the user to save for reproducibility
#' @examples data(schepCogapsResult)
#' data(schepGranges)
#'
#' motifResults = simpleMotifTFMatch(schepCogapsResult,
#'  generanges = schepGranges, organism = "Homo sapiens", genome = "hg19", scoreThreshold = 0.03)
#' @export
simpleMotifTFMatch = function(cogapsResult, generanges, organism,
                            genome, scoreThreshold = 0.03, motifsPerRegion = 1) {

  if(organism == "Homo sapiens") {
    networks = humanRegNets
  }
  else if(organism == "Mus musculus") {
    networks = mouseRegNets
  }
  else{print("Only Homo Sapiens and Mus Musculus are supported for fastMotifTFMatch. Use motifPatternMatch and downstream functions instead.", quote = F)}

  jMotifs = chromVAR::getJasparMotifs(organism)
  #call motifPatternMatch
  patternJMotifs = motifPatternMatch(cogapsResult, generanges, jMotifs,
                                     genome, scoreThreshold, motifsPerRegion)
  #make lists of all relevant info to return
  tfNameList = vector("list", length(patternJMotifs))
  motifNameList = vector("list", length(patternJMotifs))
  for(i in seq_along(patternJMotifs)) {
    motifNames = names(unlist(patternJMotifs[[i]]))
    splitNames = unlist(lapply(motifNames, stringr::str_split, "_"))
    tfNames = splitNames[which(c(1:length(splitNames))%%2==0)]
    motifNames = splitNames[which(c(1:length(splitNames))%%2==1)]
    tfNameList[[i]] = tfNames
    motifNameList[[i]] = motifNames
  }

  summs = vector("list", length(tfNameList))
  for(i in seq_along(tfNameList)) {
    patTFs = as.factor(tfNameList[[i]])
    summs[[i]]=summary(patTFs)
  }

  summs= lapply(summs, sort, decreasing = T)

  regNetList = vector("list", length(tfNameList))
  for(i in seq_along(tfNameList)) {
    patList=networks[which(names(networks) %in% tfNameList[[i]])]
    regNetList[[i]] = patList
  }

  geneDescriptionList = vector("list", length(patternJMotifs))
  for(i in seq_along(tfNameList)) {
    patAnn = humanGeneAnnotations[which(humanGeneAnnotations$Symbol %in% tfNameList[[i]]), 3:4]#make this create cleaner object
    geneDescriptionList[[i]] = patAnn
  }

  return(list(tfNames = tfNameList, motifNames = motifNameList,
              regulatoryNetworks = regNetList,
              tfDescriptions = geneDescriptionList, tfMatchSummary = summs,
              downloadedMotifs = jMotifs))

}
