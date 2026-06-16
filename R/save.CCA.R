
#' Save CCA3 results to (human readable) file
#'
#' Store results of current canonical variable in a file. 
#' First row: information about lambdas and correlation coefficients. 
#' Next rows for all three data sets for each weight > 0 weight and feature name. 
#' Weights are sorted in decreasing order, so most important features are on top.
#' 
#' @param CCA3 list with getCCA3() information to save
#' @param filename CCA-project name as base output filename (_numcv.txt be appended)
#' @param omic.names names of biological datasets (X) to separate selected features list in text file
#' 
#'  @importFrom utils tail
#'
#' @return
#' @export
#'
#' @author Andrea Thum
#' @examples
#' TRUE
save.CCA <- function(CCA3, filename, omic.names) {
  for (numcv in c(1:dim(CCA3$cc3.weight.x[[1]])[2])) {
    Sxj <- lapply(CCA3$cc3.weight.x, function(x) sort(abs(x[, numcv]), decreasing = T)) # # decreasing by weight vectors
    oS <- lapply(CCA3$cc3.weight.x, function(x) order(abs(x[, numcv]), decreasing = T))
    Sxj <- lapply(seq_along(CCA3$cc3.weight.x), function(i) {Sxj[[i]] * sign(CCA3$cc3.weight.x[[i]][, numcv])[oS[[i]]]})
    
    Szj <- sort(abs(CCA3$cc3.weight.z[, numcv]), decreasing = T)
    oS <- order(abs(CCA3$cc3.weight.z[, numcv]), decreasing = T)
    Szj <- Szj * sign(CCA3$cc3.weight.z[, numcv])[oS]
    
    nX <- length(CCA3$cc3.weight.x)
    lambdas <- paste(paste0(
      c(paste0("lambda.x", seq_len(nX)), "lambda.z"),
      ": ",
      c(CCA3$lambda[seq_len(nX), numcv],
        CCA3$lambda[nX + 1, numcv])
    ),
    collapse = " | "
    )
    
    corr <- CCA3$corr[numcv]
    fname = paste(filename, "_", numcv, ".txt", sep = "")
    
    write(file = fname,
          paste("Results spCCA, lambdas:", lambdas, "\n",
                paste(numcv, ". Variable", sep = ""), "correlation:", corr, "\n"),
          append = F
    )
    for (S in seq_along(Sxj)){
      out <- paste('### ',omic.names[S],' ###')
      write(out, file = fname, append = TRUE)
      for (j in 1:length(Sxj[[S]])) {
        if (Sxj[[S]][j] == 0.0)
          break
        out <- paste(names(Sxj[[S]])[j], Sxj[[S]][j], sep = "\t")
        write(out, file = fname, append = TRUE)
      }}
    
    write("\n", file = fname, append = TRUE)
    
    out <- paste('### Design Variables ###')
    write(out, file = fname, append = TRUE)
    for (j in 1:length(Szj)) {
      if (Szj[j] == 0.0)
        break
      out <- paste(names(Szj)[j], Szj[j], sep = "\t")
      write(out, file = fname, append = TRUE)
    }
  }
}
