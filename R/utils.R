#' Parse GMT file to gene lists
#'
#' @param path gmt file (from http://www.gsea-msigdb.org/gsea/downloads.jsp)
#' @param sep remove string to keep path and genes
#' @param rm remove string from path name
#' @return named list of genes
#' @export
gmt_to_list <- function(path,
                        sep = "\thttp.*?\t",
                        rm = "^.+?_") {
  df <- readr::read_lines(path) %>% as.data.frame() %>% setNames("X1")
  df <- tidyr::separate(df,
                        X1,
                        sep = sep,
                        into = c("path", "genes")
  ) %>% dplyr::mutate(path = stringr::str_remove(path, rm))
  
  l <- str_split(df$genes, "\t")
  names(l) <- df$path
  return(l)
}