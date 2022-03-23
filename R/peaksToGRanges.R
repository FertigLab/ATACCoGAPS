#' List of peaks to GRanges
#'
#' Wrapper function for makeGrangesFromDataFrame() from the GenomicRanges
#' package to build GRanges objects from character list of chromosomal regions
#' because this is a common format to receive peak information.
#'
#' @param region_list character list or vector of chromosomal regions/peaks in
#'   form chromosomenumber(sep)start(sep)end eg. Chr1-345678-398744
#' @param sep separator between information pieces of string (conventionally "-"
#'   or ".")
#' @return GRanges corrsponding to input list of region information
#' @section Note: If region_list is a dataframe you should use the GenomicRanges
#'   function makeGRangesFromDataFrame which this function applies
#' @examples data(schepPeaks)
#'
#' schepGranges = peaksToGRanges(schepPeaks, sep = "-")
#' @export
peaksToGRanges = function(region_list, sep) {

  #split strings into chr, start location, and end location
  splitNames = unlist(lapply(region_list, stringr::str_split, stringr::coll(sep)))
  chrs = splitNames[seq(1,length(splitNames), by = 3)]
  starts = splitNames[seq(2,length(splitNames), by = 3)]
  ends = splitNames[seq(3,length(splitNames), by = 3)]

  #convert starting and ending location info to numeric values
  starts = as.numeric(starts)
  ends = as.numeric(ends)

  #build GRanges by creating dataframe to supply to the makeGRangesFromDataFrame function
  df = data.frame(chrs, starts, ends)
  granges = GenomicRanges::makeGRangesFromDataFrame(df, ignore.strand = TRUE, seqnames.field = "chrs",
                                                    start.field = "starts",
                                                    end.field = "ends")

  return(granges)
}
