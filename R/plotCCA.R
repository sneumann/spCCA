#' Plot CCA results to PDF file(s)
#'
#' Plot Canonical Variables and top two features for each canonical variable
#' 
#' 
#' @param CCA3 list with getCCA3() information to plot
#' @param X original biological data matrices without subtraction of canonical variables, for plots of top features
#' @param Experiments Sample names for X-axis. If NULL, rownames from CCA3 are considered
#' @param filename If NULL, figures are plotted on screen, otherwise CCA-project name as base output filename (_CV_numcv.pdf and _Top_numcv.pdf will be appended)
#' @param legend names of data sets X, Z for legend, e.g. c("Genes", "Metabolites","Design")
#'
#' @importFrom graphics mtext par
#'
#' @return
#' @export
#' @author Andrea Thum
#' @examples
#' TRUE
plotCCA <- function(CCA3, X, Experiments = NULL, filename=NULL, legend) {
  if (length(legend) != length(X)+1)
    stop("Length of legend must match the total number of datasets")
  CV.X = CCA3$cc3.CV.x
  CV.Z = CCA3$cc3.CV.z

  numberCanVar <- dim(CV.X[[1]])[2] # number of can. variables
  if(is.null(Experiments))
    x_label <- rownames(CV.X[[1]])
  else
    x_label <- Experiments
  cols <- c("blue", "darkgreen", "purple", "orange","black","yellow","pink")
  max.lab <- max(nchar(x_label))
  xlab.line <- max.lab * 0.06 + 4
  par(mar = c(8, 4, 4, 2))
  
  for (numcv in 1:numberCanVar) {
    # plot latent variable
    x.cv <- lapply(CV.X, function(y) {
      y <- y[, numcv]
      y.length <- as.numeric(sqrt(t(y) %*% y))
      y <- y / y.length
      y })
    
    z.cv <- CV.Z[, numcv]
    z.cv.length <- as.numeric(sqrt(t(z.cv) %*% z.cv))
    z.cv <- z.cv / z.cv.length
    
    if (cor(z.cv, x.cv[[1]]) < 0)
      z.cv <- -z.cv
    
    maxy <- max(abs(c(unlist(x.cv), z.cv)))
    miny <- min(abs(c(unlist(x.cv), z.cv)))
    
    maxyy <- ceiling(maxy / 10 ^ (floor(log(maxy, 10)))) * (10 ^ floor(log(maxy, 10)))
    if (min(c(unlist(x.cv), z.cv)) > 0)
      minyy <- 0
    else
      minyy <- -maxyy
    
    pvx <- c()
    for (i in 1:dim(CV.X[[1]])[1])
      pvx <- c(pvx, i - 1, i)
    pvy <- c()
    for (i in 1:dim(CV.X[[1]])[1])
      pvy <- c(pvy, x.cv[[1]][i], x.cv[[1]][i])
    
    if (!is.null(filename)) {
      pdf(paste(filename, "_CV_", numcv, ".pdf", sep = ""))
    }
    
    plot(pvx + 0.4,
         pvy,
         ylim = c(minyy, maxyy),
         main = paste(numcv, ". canonical variable", sep = ""),
         type = "l",
         col = "green",
         xaxt ="n",
         xlab = "",
         ylab = "Normalized Intensity",
         axes = F,
         frame.plot = T)
    
    axis(1,
         at = pvx[seq(1, length(pvx), by = 2)] + 0.4,
         labels = x_label,
         las = 2,
         cex.axis = 0.5)
    mtext("Experiments",side=1, line = xlab.line)
    
    for(k in 2:length(x.cv)) {
      
      pvy <- rep(x.cv[[k]], each = 2)
      
      lines(pvx + 0.3 + 0.1*k,
            pvy,
            col = cols[k-1])
    }
    
    if (length(z.cv) > 0) {
      pvy <- rep(z.cv, each = 2)
      lines(pvx + 0.4 + 0.1*length(x.cv), pvy, col = "red")
    }
    
    
    axis(2, at = (c(0:6) * (maxyy - minyy)) / 5 + minyy)
    legend("topright" ,
           legend,
           cex = 0.5,
           y.intersp = 0.7,
           col = c("green", cols[1:length(x.cv)-1], "red"),
           pch = rep(1, 3),
           lty = 1:2)

    if (!is.null(filename)) {
      dev.off()
    }
    
    
    ######################################################
    #plot top features
    
    Sxj <- lapply(CCA3$cc3.weight.x,function(x) sort(abs(x[, numcv]), decreasing = T)) # decreasing by weight vectors
    
    plot1 <- lapply(Sxj, function(x) names(x)[1:2])
    dims <- sapply(plot1, function(x) length(x))
    
    Data1 <- lapply(X, function(x) t(x))
    D3_list <- lapply(seq_along(Data1), function(i)
      as.matrix(Data1[[i]][rownames(Data1[[i]]) %in% plot1[[i]], ]))
    D3 <- do.call(rbind, D3_list)
    
    starts <- c(1, cumsum(dims)[-length(dims)] + 1)
    ends <- cumsum(dims)
    
    for(k in seq_along(dims)) {
      plot1[[k]] <- rownames(D3)[starts[k]:ends[k]]
    }
    
    maxy <- max(D3)
    miny <- abs(min(0, min(D3)))
    
    maxyy <- ceiling(maxy / 10 ^ (floor(log(maxy, 10)))) * (10 ^ floor(log(maxy, 10)))
    minyy <- -ceiling(miny / 10 ^ (floor(log(miny, 10)))) * (10 ^ floor(log(miny, 10)))
    if (is.nan(minyy))
      minyy <- 0
    
    colors1 <- rainbow(length(unlist(plot1)), start = 0, end = 0.65)
    
    pvy <- rep(D3[1, ], each = 2)
    pvx <- rep(seq_len(ncol(D3)), each = 2)
    pvx[seq(1, length(pvx), by = 2)] <- pvx[seq(1, length(pvx), by = 2)] - 1
    
    if (!is.null(filename)) {
      pdf(paste(filename, "_Top_", numcv, ".pdf", sep = ""))
    }
    
    plot(pvx + 0.3,
         pvy,
         ylim = c(minyy, maxyy),
         xaxt = "n",
         xlab = "",
         ylab = "Intensity",
         col = colors1[1],
         axes = F,
         type = "l",
         frame.plot = T)
    axis(
      side = 1,
      at = pvx[seq(1, length(pvx), by = 2)] + 0.3,
      labels = x_label,
      las = 2,
      cex.axis = 0.5      
    )
    mtext("Experiments",side=1, line = xlab.line)
    axis(2, at = (c(0:6) * (maxyy - minyy)) / 5 + minyy)
    
    if (length(unlist(plot1)) > 1) {
      for (j in 2:length(unlist(plot1))) {
        pvy <- rep(D3[j, ], each = 2)
        lines(pvx + 0.3, pvy, col = colors1[j], type = "l")
      }
    }
    
    legend("topright" ,
           rownames(D3),
           cex = 0.5,
           y.intersp = 0.7,
           col = colors1,
           pch = rep(1, length(unlist(plot1))),
           lty = 1:2)
    
    if (!is.null(filename)) {
      dev.off()
    }
  }
}
