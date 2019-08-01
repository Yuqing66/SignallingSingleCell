#' create table of DE results
#'
#' This function will create a table with each gene as a row and each DE output file as a column.
#'
#' @param folders storing the output files of findDEgenes(), with cell type as facet_by.
#' @param celltype the cell type to select
#' @param input the ex_sc used for findDEgenes()
#' @param keeps the column in the result table that want to keep
#' @param filter_by the column to filter by
#' @param number rank the 'filter_by', keep the top 'number'
#' @param fc_cutoff used for logFC filter
#' @export
#' @details
#' Use large 'number' and small 'fc_cutoff' to get all the genes kept.
#' @examples
#'
function(folders, celltype, input, keeps = c("logFC", "FDR"), filter_by = "FDR", number = 30000, fc_cutoff = 0.000001){
  types <- sort(unique(pData(input)$CellType))
  lists <- vector(mode = "list", length = length(types))
  for (l in 1:length(folders)) {
    DEfolder <- folders[l]

    detabs <- vector(mode = "list", length = length(types))
    filelist <- list.files(DEfolder)
    filelist <- paste0(DEfolder, "/", filelist)

    detab <- matrix(ncol = length(filelist) / length(types) * (length(keeps) + 1),nrow = nrow(input))

    colname_detab <- c()
    comparisons <- filelist[grep(types[1], filelist)]
    for (i in 1:length(comparisons)) {
      comparison <- comparisons[i]
      comparison <- strsplit(comparison, paste0(types[1], "_"))[[1]][2]
      comparison <- gsub(".txt", "", comparison)
      colname_detab <- c(colname_detab, comparison, keeps)
    }
    colnames(detab) <- colname_detab
    rownames(detab) <- rownames(input)

    for (k in 1:length(detabs)) {
      detabs[[k]] <- detab
    }

    names(detabs) <- types

    for (i in 1:length(filelist)) {
      int <- filelist[i]
      interested <- c()
      # Get the CellType of the file
      for (j in 1:length(types)) {
        match <- grep(types[j], int)
        if (length(match > 0)) {
          if (match == 1) {
            interested <- types[j]
          }
        }
      }

      comparison <- strsplit(int, paste0(interested, "_"))[[1]][2]
      comparison <- gsub(".txt", "", comparison)
      #comparison <- gsub(paste0("_",folder), "", comparison)
      list_index <- match(interested, names(detabs))
      col_index <- match(comparison, colnames(detab))
      tmp <- read.table(int,
                        sep = "\t",
                        row.names = 1,
                        header = T)

      # Filter the significant genes
      tmp <- tmp[abs(tmp[, "logFC"]) > fc_cutoff, ]
      if (nrow(tmp) < number) {
        int_gene <- rownames(tmp)[order(tmp[, filter_by])][1:nrow(tmp)]
        print(paste0(
          "high fold change genes < top genes required:",
          nrow(tmp),
          " in ",
          int
        ))
      } else{
        int_gene <- rownames(tmp)[order(tmp[, filter_by])][1:number]
      }

      final_ind <- match(int_gene, rownames(detabs[[list_index]]))
      detabs[[list_index]][, col_index] <- 0
      detabs[[list_index]][final_ind, col_index] <- 1
      for (k in 1:length(keeps)) {
        detabs[[list_index]][final_ind, col_index + k] <-
          tmp[int_gene, keeps[k]]
      }
    }

    # delete non-significant genes
    for (i in 1:length(detabs)) {
      inter <- detabs[[i]]
      if (ncol(inter) == 1 + length(keeps)) {
        sums <- inter[, 1]
      } else{
        sums <- apply(inter[, grep(keeps[1], colnames(inter)) - 1], 1, sum)
      }
      remove <- which(sums == 0)
      inter <- inter[-remove, ]
      detabs[[i]] <- inter
    }

    # 'bin' column is deleted from the original filterDEgenes.R

    lists[[l]] <- detabs
  }
  names(lists) <- folders

  ###### rearrange the files ######
  filenames <- folders
  listnames <- names(lists[[1]])
  lists2 <- c()
  for (i in 1:length(listnames)) {
    for (j in 1:length(filenames)) {
      lists2[[listnames[i]]][j] <- lists[[j]][listnames[i]]
    }
    names(lists2[[listnames[i]]]) <- filenames
    #saveRDS(lists2[[listnames[i]]], file = paste0("DEgenes_", listnames[i], ".rds"))
  }

  ###### merge into a table ######
  #celltype <- "TC"
  data <- lists2[[celltype]]
  genes <- c()
  for (i in 1:length(data)) {
    genes <- union(genes, rownames(data[[i]]))
  }
  table <-
    matrix(nrow = length(genes), ncol = length(data) * length(keeps))
  rownames(table) <- genes
  colnames(table) <-
    paste0(rep(names(data), each = length(keeps)), "_", rep(keeps, length(data)))
  for (i in 1:length(data)) {
    table[rownames(data[[i]]), (length(keeps) * (i - 1) + 1):(length(keeps) *i)] <- data[[i]][, 2:(1 + length(keeps))]
  }
  return(table)
}
