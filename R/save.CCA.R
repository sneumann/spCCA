
#' Save CCA3 results to (human reabable) file
#'
#' @param CCA3 list with CCA3 information to save
#' @param filename base output filename (_numcv.txt be appended)
#'
#' @return
#' @export
#'
#' @examples
#' TRUE
save.CCA <- function(CCA3, filename) {
  for (numcv in c(1:dim(CCA3$cc3.weight.x)[2])) {
    Sxj <- sort(abs(CCA3$cc3.weight.x[, numcv]), decreasing = T) # # decreasing by weight vectors
    oS <- order(abs(CCA3$cc3.weight.x[, numcv]), decreasing = T)
    Sxj <- Sxj * sign(CCA3$cc3.weight.x[, numcv])[oS]
    
    Syj <- sort(abs(CCA3$cc3.weight.y[, numcv]), decreasing = T)
    oS <- order(abs(CCA3$cc3.weight.y[, numcv]), decreasing = T)
    Syj <- Syj * sign(CCA3$cc3.weight.y[, numcv])[oS]
    
    Szj <- sort(abs(CCA3$cc3.weight.z[, numcv]), decreasing = T)
    oS <- order(abs(CCA3$cc3.weight.z[, numcv]), decreasing = T)
    Szj <- Szj * sign(CCA3$cc3.weight.z[, numcv])[oS]
    
    lambdas <- paste("lambda.x:", CCA3$lambda[1, numcv],
                     "lambda.y:", CCA3$lambda[2, numcv],
                     "lambda.z:", CCA3$lambda[3, numcv])
    corr <- CCA3$corr[numcv]
    fname = paste(filename, "_", numcv, ".txt", sep = "")
    
    write(file = fname,
      paste("Results spCCA, lambdas:", lambdas, "\n",
        paste(numcv, ". Variable", sep = ""), "correlation:", corr, "\n"),
      append = F
    )
    for (j in 1:length(Sxj)) {
      if (Sxj[j] == 0.0)
        break
      out <- paste(names(Sxj)[j], Sxj[j], sep = "\t")
      write(out, file = fname, append = TRUE)
    }
    
    write("\n", file = fname, append = TRUE)
    for (j in 1:length(Syj)) {
      if (Syj[j] == 0.0)
        break
      out <- paste(names(Syj)[j], Syj[j], sep = "\t")
      write(out, file = fname, append = TRUE)
    }
    write("\n", file = fname, append = TRUE)
    for (j in 1:length(Szj)) {
      if (Szj[j] == 0.0)
        break
      out <- paste(names(Szj)[j], Szj[j], sep = "\t")
      write(out, file = fname, append = TRUE)
    }
  }
}
