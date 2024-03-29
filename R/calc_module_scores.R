#' Calculate Module/Pathway Scores , Parallel on R side
#'
#' @param mat expression matrix
#' @param paths named lists, each with genes of one pathway
#' @param nbin number of bins to arrange all genes
#' @param nsample number of genes to sample for each target gene
#' @param nthread number of threads for parallel processing of many pathways
#' @param precal option on how background order is calculated, set to TRUE to use same calculations for every pathway, pass vector to use specific order
#' @param nrep number of repeats to calculate sampling of bins
#' @param rust use rust (faster but no sparse matrix support) implementations
#' @return matrix of pathways x cells
#' @export
calc_module_scores <- function(mat,
                        paths,
                        nbin = 24,
                        nsample = 100,
                        nthread = 4,
                        precal = T,
                        nrep = 20,
                        rust = T) {
  options(future.globals.maxSize= 89128960000)
  gs <- rownames(mat)
  paths <- map(paths, function(x) {
    intersect(gs, x)
  })
  paths <- paths[map(paths, function(x) {length(x) > 5}) %>% unlist()]
  if (precal != FALSE) {
    if (length(precal) > 1) {
      message("passing precal argument as order...")
      gorder <- precal
    } else {
      if (rust) {
        gorder <- order_expr(mat, gs, 4)
      } else {
        gorder <- order_expr_r(mat, gs)
      }
    }
    if (length(paths) >= 50) {
      message(">= 50 pathways detected, calculating background first to avoid repeat calculations...")
      avgs <- 1:length(gorder)
      names(avgs) <- gorder
      bins <- cut_number(x = avgs + rnorm(n = length(avgs))/1e30, n = nbin, labels = FALSE, right = FALSE)
      names(bins) <- names(avgs)
      bac <- calc_background_sampling(mat,
                                      bins,
                                      nsample = nsample,
                                      nrep = nrep)
      pb <- progressr::progressor(steps = length(paths))
      res <- furrr::future_map(paths,.options = furrr::furrr_options(seed = TRUE),  function(x) {
        pb()
        pos <- Matrix::colMeans(mat[x,])
        ctrl <- get_background(bac, x, bins, nrep)
        pos - ctrl
      })
      
    } %>% progressr::with_progress() else {
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

#' Calculate Module/Pathway background, to be used
#'
#' @param mat expression matrix
#' @param bins gene bin vector to use specific order
#' @param nsample number of genes to sample for each target gene
#' @param nrep number of repeats to calculate sampling of bins
#' @return matrix of background bin reps x cells
#' @export
calc_background_sampling <- function(mat,
                                     bins,
                                     nsample = 100,
                                     nrep = 20) {
  data.cut <- bins
  bac <- matrix(
    data = numeric(length = 1L),
    nrow = nrep * max(data.cut),
    ncol = ncol(mat)
  )
  for (x in 1:max(data.cut)) {
    for (i in 1:nrep) {
      temp_bin <- data.cut[which(data.cut == x)]
      temp <- names(x = sample(
        x = temp_bin,
        size = nsample,
        replace = FALSE
      ))
      bac[i + nrep * (x - 1), ] <- Matrix::colMeans(x = mat[temp, ])
    }
  }
  bac
}

#' Construct a Module/Pathway background matrix from gene list
#'
#' @param background calculated background expression matrix
#' @param features vector of genes in pathway
#' @param bins gene bin vector to use specific order
#' @param nrep number of repeats to calculate sampling of bins
#' @return matrix of sampled background # of feature genes x cells
#' @export
get_background <- function(background,
                           features,
                           bins,
                           nrep) {
  bac1 <- matrix(
    data = numeric(length = 1L),
    nrow = length(features),
    ncol = ncol(background)
  )
  for (j in 1:length(x = features)) {
    ind <- sample(nrep, 1) + nrep * (bins[features[j]] - 1)
    bac1[j, ] <- background[ind,]
  }
  Matrix::colMeans(bac1)
}

#' Order genes by overall expression
#'
#' @param mat expression matrix
#' @return vector of gene names, high to low
#' @export
order_expr_r <- function(mat, 
                         gs) {
  data.avg <- Matrix::rowMeans(x = mat[gs,])
  data.avg[order(data.avg)] %>% names()
}