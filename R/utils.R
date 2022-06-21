#' Parse GMT format to get list of pathways
#'
#' @param gmt GMT file
#' @param rm characters to remove from entry names
#' @param nmin minimum number of genes to keep pathway
#' @return named list of vectors, each representing a pathway
#' @export
calc_module_scores <- function(mat,
                        paths,
                        nbin = 24,
                        nsample = 100,
                        nthread = 4) {
  gs <- rownames(mat)
  res <- furrr::future_map(paths, function(x) {
    calc_modulescore(mat,
                     x,
                     gs,
                     nbin,
                     nsample,
                     nthread)
  })
  res <- do.call(rbind, res)
  colnames(res) <- colnames(mat)
  rownames(res) <- names(paths)
  res
}

fast_bind_rows_dt <- function(l) {
  l2 <- map(l, function(x) {as.list(x)})
  data.table::rbindlist(l2)
}
