

#legend: e.g. c("Genes","Metabolites","Patterns")
#




#' Title
#'
#' @param CCA3 
#' @param X 
#' @param Y 
#' @param filename 
#' @param legend 
#'
#' @return
#' @export
#'
#' @examples
#' TRUE
plotCCA <- function(CCA3, X, Y, filename, legend) {
  CV.X = CCA3$cc3.CV.x
  CV.Y = CCA3$cc3.CV.y
  CV.Z = CCA3$cc3.CV.z
  numberCanVar <- dim(CV.X)[2]#number of can. variables
  for (numcv in 1:numberCanVar) {
    #plot latent variable
    x.cv <- CV.X[, numcv]
    x.cv.length <- as.numeric(sqrt(t(x.cv) %*% x.cv))
    x.cv <- x.cv / x.cv.length
    
    y.cv <- CV.Y[, numcv]
    y.cv.length <- as.numeric(sqrt(t(y.cv) %*% y.cv))
    y.cv <- y.cv / y.cv.length
    
    z.cv <- CV.Z[, numcv]
    z.cv.length <- as.numeric(sqrt(t(z.cv) %*% z.cv))
    z.cv <- z.cv / z.cv.length
    
    if (cor(z.cv, x.cv) < 0)
      z.cv <- -z.cv
    
    maxy <- max(abs(c(x.cv, y.cv, z.cv)))
    miny <- min(abs(c(x.cv, y.cv, z.cv)))
    
    
    maxyy <- ceiling(maxy / 10 ^ (floor(log(maxy, 10)))) * (10 ^ floor(log(maxy, 10)))
    if (min(c(x.cv, y.cv, z.cv)) > 0)
      minyy <- 0
    else
      minyy <- -maxyy
    
    pvx <- c()
    for (i in 1:dim(CV.X)[1])
      pvx <- c(pvx, i - 1, i)
    pvy <- c()
    for (i in 1:dim(CV.X)[1])
      pvy <- c(pvy, x.cv[i], x.cv[i])
    
    pdf(paste(filename, "_CV_", numcv, ".pdf", sep = ""))
    plot(
      pvx + 0.4,
      pvy,
      ylim = c(minyy, maxyy),
      main = paste(numcv, ". canonical variable", sep = ""),
      type = "l",
      col = "green",
      xlab = "Experiments",
      ylab = "Normalized Intensity",
      axes = F,
      frame.plot = T
    )
    
    pvy <- c()
    for (i in 1:dim(CV.X)[1])
      pvy <- c(pvy, y.cv[i], y.cv[i])
    
    lines(pvx + 0.5, pvy, col = "blue")
    if (length(z.cv) > 0) {
      pvy <- c()
      for (i in 1:dim(CV.X)[1])
        pvy <- c(pvy, z.cv[i], z.cv[i])
      lines(pvx + 0.6, pvy, col = "red")
    }
    
    axis(2, at = (c(0:6) * (maxyy - minyy)) / 5 + minyy)
    legend(
      "topright" ,
      legend,
      cex = 0.8,
      col = c("green", "blue", "red"),
      pch = rep(1, 3),
      lty = 1:2
    )
    
    dev.off()
    ######################################################
    
    
    #plot top features
    
    
    Sxj <- sort(abs(CCA3$cc3.weight.x[, numcv]), decreasing = T)#absteigend nach Gewicht die Vektoren
    Syj <- sort(abs(CCA3$cc3.weight.y[, numcv]), decreasing = T)
    
    plot2 <- names(Syj)[1:2]
    plot1 <- names(Sxj)[1:2]
    
    Data1 <- t(X)
    Data2 <- t(Y)
    D3 <- as.matrix(Data1[rownames(Data1) %in% plot1, ])
    D3 <- rbind(D3, as.matrix(Data2[rownames(Data2) %in% plot2, ]))
    plot1 <- rownames(D3)[1:length(plot1)]
    plot2 <- rownames(D3)[(length(plot1) + 1):(length(plot2) + length(plot1))]
    
    maxy <- max(D3)
    miny <- abs(min(0, min(D3)))
    
    maxyy <- ceiling(maxy / 10 ^ (floor(log(maxy, 10)))) * (10 ^ floor(log(maxy, 10)))
    minyy <- -ceiling(miny / 10 ^ (floor(log(miny, 10)))) * (10 ^ floor(log(miny, 10)))
    if (is.nan(minyy))
      minyy <- 0
    
    colors1 <- rainbow(length(plot1), start = 0.4, end = 0.65)
    colors2 <- rainbow(length(plot2), start = 0, end = 0.15)
    
    pvy <- c()
    for (i in 1:dim(D3)[2])
      pvy <- c(pvy, D3[1, ][i], D3[1, ][i])
    pvx <- c()
    for (i in 1:dim(D3)[2])
      pvx <- c(pvx, i - 1, i)
    
    pdf(paste(filename, "_Top_", numcv, ".pdf", sep = ""))
    plot(
      pvx + 0.3,
      pvy,
      ylim = c(minyy, maxyy),
      xlab = "Experiments",
      ylab = "Intensity",
      col = colors1[1],
      axes = F,
      type = "l",
      frame.plot = T
    )
    
    axis(2, at = (c(0:6) * (maxyy - minyy)) / 5 + minyy)
    if (length(plot1) > 1) {
      for (j in 2:length(plot1)) {
        pvy <- c()
        for (i in 1:dim(D3)[2])
          pvy <- c(pvy, D3[j, ][i], D3[j, ][i])
        lines(pvx + 0.3, pvy, col = colors1[j], type = "l")
      }
    }
    if (length(plot2) > 1) {
      for (j in 1:length(plot2)) {
        pvy <- c()
        for (i in 1:dim(D3)[2])
          pvy <- c(pvy, D3[j + length(plot1), ][i], D3[j + length(plot1), ][i])
        lines(pvx + 0.3, pvy, col = colors2[j], type = "l")
      }
    }
    legend(
      "topright" ,
      rownames(D3),
      cex = 0.8,
      col = c(colors1[1:length(plot1)], colors2[1:length(plot2)]),
      pch = rep(1, length(c(plot1, plot2))),
      lty = 1:2
    )
    
    dev.off()
  }
  
}
