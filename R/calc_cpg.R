#' Calculate the cells per gene
#'
#' This function will calculate how many cells have non-zero expression values of each gene.
#'
#' @param ex_sc an expression set
#' @export
#' @details
#'
#' @examples
#' calc_cpg(input_norm)
function(ex_sc){
  cells_per_gene <- c()
  for (i in 1:nrow(exprs(ex_sc))){
    cells_per_gene[i] <- length(which(exprs(ex_sc)[i,] != 0))
  }
  names(cells_per_gene) <- rownames(exprs(ex_sc))
  return(cells_per_gene)
}
