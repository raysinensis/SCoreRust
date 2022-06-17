#' Caluculate Module/Pathway Scores , Parallel on R side
#'
#' @param mat expression matrix
#' @param paths named lists, each with genes of one pathway
#' @param nbin number of bins to arrange all genes
#' @param nsample number of genes to sample for each target gene
#' @param nthread number of threads for parallel processing of many pathways
#' @param precal option on how background order is calculated, set to TRUE to use same calculations for every pathway, pass vector to use specific order
#' @return matrix of pathways x cells
#' @export
calc_module_scores <- function(mat,
                        paths,
                        nbin = 24,
                        nsample = 100,
                        nthread = 4,
                        precal = T) {
  gs <- rownames(mat)
  paths <- map(paths, function(x) {
    intersect(gs, x)
  })
  paths[map(paths, function(x) {length(x) > 5}) %>% unlist()]
  if (precal != FALSE) {
    if (length(precal) >1 ) {
      message("passing precal argument as order...")
      res <- furrr::future_map(paths, function(x) {
        calc_modulescore_orderin(mat,
                                 x,
                                 gs,
                                 nbin,
                                 nsample,
                                 nthread,
                                 precal)
      })
    } else {
      gorder <- order_expr(mat, gs, 4)
      res <- furrr::future_map(paths, function(x) {
        calc_modulescore_orderin(mat,
                                 x,
                                 gs,
                                 nbin,
                                 nsample,
                                 nthread,
                                 gorder)
      })
    }
  } else {
    res <- furrr::future_map(paths, function(x) {
      calc_modulescore(mat,
                             x,
                             gs,
                             nbin,
                             nsample,
                             nthread)
    })
  }

  res <- do.call(rbind, res)
  colnames(res) <- colnames(mat)
  rownames(res) <- names(paths)
  res
}

fast_bind_rows_dt <- function(l) {
  l2 <- map(l, function(x) {as.list(x)})
  data.table::rbindlist(l2)
}
