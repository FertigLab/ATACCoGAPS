#' Filter scATACseq by sparsity
#'
#' Function to filter a set of scATACseq data by sparsity and return a subset of
#' filtered data, as well as list of the remaining cells and peaks.
#'
#' @param data matrix of read counts peaks x cells
#' @param cell_list list of cell names/identifiers for the data
#' @param peak_list list of peaks from the data
#' @param cell_cut threshold of sparsity to filter at (eg. 0.99 filters all
#'   cells with more than 99 percent zero values)
#' @param peak_cut threshold of sparsity to filter at for peaks
#' @return nested list containing the subset data, a list of peaks, and list of
#'   cells
#' @examples data("subsetSchepData")
#' data("schepPeaks")
#' data("schepCellTypes")
#' 
#' outData = dataSubsetBySparsity(subsetSchepData, schepCellTypes, schepPeaks)
#'   
#' @export
dataSubsetBySparsity = function(data, cell_list, peak_list, cell_cut = 0.99, peak_cut = 0.99) {

  #find the sparsity for each peak in the data
  peaks_sparsity = apply(data, 1, function(x){sum(x==0)/length(x)})

  #print the number of peaks being filtered from the data
  print(paste("Peaks to be filtered:", sum(peaks_sparsity > peak_cut)), quote = FALSE)


  #find the sparsity for each cell in the data
  cells_sparsity = apply(data, 2, function(x){sum(x==0)/length(x)})

  #print the number of cells being filtered from the data
  print(paste("Cells to be filtered:", sum(cells_sparsity > cell_cut)), quote = FALSE)

  #create new list of peaks that will be kept in output data
  peaks_to_remove = which(peaks_sparsity > peak_cut)
  if(length(peaks_to_remove) == 0) {
    peaks_sub = peak_list
  }
  else{
  peaks_sub = peak_list[-peaks_to_remove]
  }

  #create list of celltypes that will be kept in subset data
  cells_to_remove = which(cells_sparsity > cell_cut)
  if(length(cells_to_remove) == 0){
    celltypes_sub = cell_list
  }
  else{
    celltypes_sub = cell_list[-cells_to_remove]
  }
  celltypes_list_sub = as.character(celltypes_sub)

  #subset the input data
  if(length(cells_to_remove) == 0 & length(peaks_to_remove) == 0) {
    data_sub = data
  }
  else if(length(cells_to_remove) == 0) {
    data_sub = data[-peaks_to_remove,]
  }
  else if(length(peaks_to_remove) == 0) {
    data_sub = data[,-cells_to_remove]
  }
  else{
    data_sub = data[-peaks_to_remove, -cells_to_remove]
  }

  return(list(subset_data = data_sub, subset_cells = celltypes_list_sub,
              subset_peaks = peaks_sub))
}
