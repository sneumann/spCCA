#' compare patterns of two eigenvectors
#' 
#' return manhattan (default) or euklidean distance
#'
#' @param x,y vectors to compare  
#' @param dm distance measure, default is manhattan
#'
#' @return manhattan distance by default; other distances possible
#' @export
#'
#' @examples
comparePattern <- function(x, y, dm = c("manhattan", "euklidean"))
{
  v <- x - y
  if (dm == "manhattan")
    len <- sum(abs(v))
  if (dm == "euklidean")
    len <- sum(v * v)
  return(len)
}

#' cluster weight vectors by assigning them to one pattern
#' The addPattern function clusters weight vectors by assigning them to one pattern 
#' in the given list. It checks if the mean cluster vector is similar to a new z.vector, 
#' and makes a new cluster if not similar to any existing cluster, 
#' counts the number of vectors associated with one cluster, 
#' and saves the correlation for the cluster. 
#' The function returns the updated list of patterns.
#'
#' The addPattern function is called within the get.best.lambdas
#'
#' @param List.Pattern The list of existing patterns.
#' @param corr (STN: LLM says:) The correlation value.
#' @param xyz.vector (STN: LLM says:) The combined vector of x, y, and z-weight.
#' @param z.vector The vector for comparison to find or create clusters
#'
#' @return 
#'   List of cluster items:
#'   item [[1]]: correlation
#'   item [[2]]: combined vector of x, y, and z-weight
#'   item [[3]]: counter for number of vectors in cluster
#' @export
#'
#' @examples

addPattern <- function(List.Pattern, corr, xyz.vector, z.vector) {
  if (length(List.Pattern) == 0) {
    List.Pattern <- list(list(corr, 
                              matrix(xyz.vector, 
                                     ncol = 1, 
                                     nrow = length(xyz.vector)), 
                              1))
  } else {
    found <- FALSE
    for (i in 1:length(List.Pattern)) {
      items <- List.Pattern[[i]]
      s <- dim(items[[2]])[1] - length(z.vector) + 1
      e <- dim(items[[2]])[1]
      if (comparePattern(z.vector, rowMeans(as.matrix(items[[2]][s:e, ]))) <= 0.5) {
        # same pattern
        found <- TRUE
        List.Pattern[[i]][[1]] = List.Pattern[[i]][[1]] + corr
        List.Pattern[[i]][[2]] <- cbind(List.Pattern[[i]][[2]], xyz.vector)
        List.Pattern[[i]][[3]] <- List.Pattern[[i]][[3]] + 1
        break
      }
    }
    if (!found) {
      List.Pattern <- append(List.Pattern, list(list(corr, 
                                                     matrix(xyz.vector,
                                                            ncol = 1,
                                                            nrow = length(xyz.vector)), 
                                                     1)))
    }
  }
  return(List.Pattern)
}
