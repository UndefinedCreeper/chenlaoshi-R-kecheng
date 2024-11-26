# Build own functions here ----

# calculate chi-squared statistic using R base
# Oi is observed vector, and Ei is expected proportional value

#首先

cal_x2 <- function(Oi, Ei_prop){
  # Oi is the observed value
  # Ei is expected proportional value
  t <- sum(Oi) # total number of observation
  r <- Ei_prop/sum(Ei_prop) # ratio/probability
  exp <- r*t # expected value
  chi_squared <- sum((Oi - exp)^2 / exp)# Calculate chi-squared statistic
  
  # Print the result
  return(chi_squared)
}


# panelutils.R 
#
# License: GPL-2 
# Author:  Francois Gillet
#          23 August 2012

## Put Pearson, Spearman or Kendall correlations on the upper panel
panel.cor <- function(x, y, method="pearson", digits=3, cex.cor=1.2, no.col=FALSE)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method=method)
  ra <- cor.test(x, y, method=method)$p.value
  txt <- round(r, digits)
  prefix <- ""
  if(ra <= 0.1) prefix <- "."
  if(ra <= 0.05) prefix <- "*"
  if(ra <= 0.01) prefix <- "**"
  if(ra <= 0.001) prefix <- "***"
  if(no.col)
  {
    color <- 1
    if(r < 0) { if(ra <= 0.001) sig <- 4 else sig <- 3 }
    else { if(ra <= 0.001) sig <- 2 else sig <- 1 }
  }
  else
  {
    sig <- 1
    if(ra <= 0.001) sig <- 2
    color <- 2
    if(r < 0) color <- 4
  }
  txt <- paste(txt, prefix, sep="\n")
  text(0.5, 0.5, txt, cex = cex.cor, font=sig, col=color)
}


## Put histograms on the diagonal
panel.hist <- function(x, no.col=FALSE, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  his <- hist(x, plot=FALSE)
  breaks <- his$breaks; nB <- length(breaks)
  y <- his$counts
  y <- y/max(y)
  if(no.col) rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
  else rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


## Add black lowess curves to scatter plots
panel.smoothb <- function (x, y, col=par("col"), bg=NA, pch=par("pch"), 
                           cex=1, col.smooth="black", span=2/3, iter=3, ...) 
{
  points(x, y, pch=pch, col=col, bg=bg, cex=cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f=span, iter=iter), col=col.smooth, ...)
}


#Usage:
#pairs(num.mat, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)
#pairs(num.mat, lower.panel=panel.smooth, upper.panel=panel.cor, method="kendall")

# coldiss()
# Color plots of a dissimilarity matrix, without and with ordering
#
# License: GPL-2
# Author:  Francois Gillet
#          23 August 2012 - rev. 07 June 2016

"coldiss" <- function(D,
                      nc = 4,
                      byrank = TRUE,
                      diag = FALSE) {
  require(gclus)
  
  D <- as.dist(as.matrix(D))
  
  if (max(D) > 1)
    D <- D / max(D)
  
  if (byrank) {
    spe.color <- dmat.color(1 - D, cm.colors(nc))
  }
  else {
    spe.color <- dmat.color(1 - D, byrank = FALSE, cm.colors(nc))
  }
  
  spe.o <- order.single(1 - D)
  speo.color <- spe.color[spe.o, spe.o]
  
  op <- par(mfrow = c(1, 2), pty = "s")
  
  if (diag) {
    plotcolors(
      spe.color,
      rlabels = attributes(D)$Labels,
      main = "Dissimilarity Matrix",
      dlabels = attributes(D)$Labels
    )
    plotcolors(
      speo.color,
      rlabels = attributes(D)$Labels[spe.o],
      main = "Ordered Dissimilarity Matrix",
      dlabels = attributes(D)$Labels[spe.o]
    )
  }
  else {
    plotcolors(spe.color, rlabels = attributes(D)$Labels,
               main = "Dissimilarity Matrix")
    plotcolors(speo.color,
               rlabels = attributes(D)$Labels[spe.o],
               main = "Ordered Dissimilarity Matrix")
  }
  
  par(op)
}

# Usage:
# coldiss(D = dissimilarity.matrix, nc = 4, byrank = TRUE, diag = FALSE)
# If D is not a dissimilarity matrix (max(D) > 1), then D is divided by max(D)
# nc 							number of colours (classes)
# byrank = TRUE		equal-sized classes
# byrank = FALSE	equal-length intervals
# diag = TRUE			print object labels also on the diagonal

# Example:
# coldiss(spe.dj, nc = 9, byrank = FALSE, diag = TRUE)


# Map of the Doubs river (see Chapter 2)
# License: GPL-2
# Author:  Francois Gillet

drawmap <-
  function(xy = spa, clusters, main = "Clusters along the Doubs river", tcol = "white") {
    
    # Draw the Doubs river
    plot(
      xy,
      asp = 1,
      type = "n",
      main = main,
      xlab = "x coordinate (km)",
      ylab = "y coordinate (km)"
    )
    lines(xy, col = "light blue")
    text(65, 20, "Upstream", cex = 1.2)
    text(15, 32, "Downstream", cex = 1.2)
    
    # Add the clusters
    k <- length(levels(factor(clusters)))
    for (i in 1:k)
    {
      points(
        xy[clusters == i, 1],
        xy[clusters == i, 2],
        #        pch = i + 20,
        pch = i + 15,
        cex = 3,
        col = i + 1,
        bg = i + 1
      )
    }
    text(xy,
         row.names(xy),
         cex = 0.8,
         col = tcol,
         font = 2)
    legend(
      "bottomright",
      paste("Cluster", 1:k),
      #      pch = (1:k) + 20,
      pch = (1:k) + 15,
      col = 2:(k + 1),
      pt.bg = 2:(k + 1),
      pt.cex = 2,
      bty = "n"
    )
    
  }


# Map of the Doubs river (see Chapter 4)
# License: GPL-2
# Author:  Francois Gillet

drawmap3 <-
  function(xy = spa, clusters, main = "Clusters along the Doubs river", 
           colors = palette()[-1], pch = 21, tcol = "black") {
    
    # Draw the Doubs river
    plot(
      xy,
      asp = 1,
      type = "n",
      main = main,
      xlab = "x coordinate (km)",
      ylab = "y coordinate (km)"
    )
    lines(xy, col = "light blue")
    text(65, 20, "Upstream", cex = 1.2)
    text(15, 32, "Downstream", cex = 1.2)
    
    # Add the clusters
    k <- length(levels(factor(clusters)))
    for (i in 1:k)
    {
      points(
        xy[clusters == i, 1],
        xy[clusters == i, 2],
        pch = pch,
        cex = 3,
        col = "white",
        bg = colors[i]
      )
    }
    text(xy,
         row.names(xy),
         cex = 0.8,
         col = tcol,
         font = 2)
    legend(
      "bottomright",
      paste("Cluster", 1:k),
      pch = 22,
      col = "white",
      pt.bg = colors,
      pt.cex = 2,
      bty = "n"
    )
    
  }

# Function hcoplot()
# Reorder and plot dendrogram with colors for groups and legend
#
# Usage:
# hcoplot(tree = hclust.object, diss = dissimilarity.matrix, k = nb.clusters, 
#	title = paste("Reordered dendrogram from",deparse(tree$call),sep="\n"))
#
# License: GPL-2 
# Author:  Francois Gillet, 23 August 2012
# Revised: Daniel Borcard, 31 August 2017

"hcoplot" <- function(tree, 
                      diss, 
                      lab = NULL,
                      k, 
                      title = paste("Reordered dendrogram from", 
                                    deparse(tree$call), 
                                    sep="\n"))
{
  require(gclus)
  gr <- cutree(tree, k=k)
  tor <- reorder.hclust(tree, diss)
  plot(tor, 
       labels = lab,
       hang=-1, 
       xlab=paste(length(gr),"sites"), 
       sub=paste(k,"clusters"), 
       main=title)
  so <- gr[tor$order]
  gro <- numeric(k)
  for (i in 1 : k)
  {
    gro[i] <- so[1]
    if (i<k) so <- so[so!=gro[i]]
  }
  rect.hclust(tor, 
              k = k, 
              border = gro + 1, 
              cluster = gr)
  legend("topright", 
         paste("Cluster", 1 :k ), 
         pch = 22, 
         col = 2 : (k + 1), 
         bty = "n")
}


# Function to compute and test the statistic 'a' by permutation.
# 'a' is the co-occurrence of two species across the sites.
#
# Parameters:
# mat = site-by-species matrix. Sites are rows, species are columns.
# nperm = number of permutations. Choose this number so that the smallest
#         p-values will remain significant after correction for multiple
#         testing, e.g. Holm correction
#
# License: GPL-2
# Author:  Pierre Legendre
#          2010

test.a <-
  function (mat, nperm = 999) {
    A <- system.time({
      mat <- as.matrix(mat)
      site.names <- rownames(mat)
      sp.names <- colnames(mat)
      
      # Transform the 'pa' or abundance data to presence-absence
      mat <- decostand(mat, "pa")
      
      n <- nrow(mat)
      p <- ncol(mat)
      a <- t(mat) %*% mat
      
      # Permutation tests for a
      p.a = matrix(1, p, p, dimnames = list(sp.names, sp.names))
      for (iperm in 1:nperm) {
        perm.mat = mat[sample(n), 1]
        for (j in 2:p) {
          vec <- mat[sample(n), j]
          perm.mat <- cbind(perm.mat, vec)
        }
        #
        a.perm <- t(perm.mat) %*% perm.mat
        for (j in 2:p) {
          for (jj in 1:(p - 1)) {
            if (a.perm[j, jj] >= a[j, jj])
              p.a[j, jj] <- p.a[j, jj] + 1
          }
        }
      }
      p.a <- p.a / (nperm + 1)
      
      for (j in 1:(p - 1)) {
        for (jj in (j + 1):p) {
          p.a[j, jj] <- NA
        }
      }
      diag(p.a) <- NA
      
    })
    A[3] <- sprintf("%2f", A[3])
    cat("Computation time =", A[3], " sec", '\n')
    
    out <- list(a = a,
                p.a = p.a,
                p.a.dist = as.dist(p.a))
    out
  }


bartlett.perm <- function(y, fact, centr="MEDIAN", nperm=999, alpha=0.05)
  
  # Computation of parametric, permutational and bootstrap versions of the
  # Bartlett test of homogeneity of variances.
  #
  # The data are centred to their within-group medians (default) or means 
  # before the tests.
  #
  # Prior to the computation of the test of homogeneity of variances,
  # a Shapiro-Wilk test of normality of residuals is computed. If the residuals
  # are not normally distributed, a warning is issued because this nonnormality
  # influences the type I error of the parametric test, which will very likely 
  # have an incorrect type I error.
  #
  # USAGE
  # bartlett.perm(y, fact, centr, nperm, alpha)
  #
  # ARGUMENTS
  # y       a numeric vector of data values
  # fact    a vector or factor object giving the group for the corresponding
  #         elements of y
  # centr   should the data, within groups, be centred on their medians ("MEDIAN")
  #         or on their means ("MEAN")?
  # nperm   number of permutations
  # alpha   level of rejection of the H0 hypothesis of normality of residuals
  #         in the Shapiro-Wilk test
  #
  # RESULT
  # Bartlett          Bartlett's K-squared test statistic
  # Param.prob        Parametric probability (P-value) of the test
  # Permut.prob       Permutational probability (P-value) of the test
  # Bootstrap.prob    Bootstrap probability (P-value) of the test
  #
  # DETAILS
  #
  # Centring the groups on their median or mean is very important for permutation
  # and bootstrap tests to be correct when the groups do not share the same 
  # position. Permuting groups with unequal mean or median artificially increases 
  # the within-group variances of the permuted data.
  #
  # License: GPL-2
  # Author:  Daniel Borcard
  #          1 February 2016
  
  # ------------------------------------------------------
# EXAMPLE
#
# Species abundance type data:
# y1 <- log1p(round(rlnorm(5,0.2,2)))
# y2 <- log1p(round(rlnorm(5,1,2)))
# y3 <- log1p(round(rlnorm(5,2,5)))
# yy <- c(y1,y2,y3)
#
# Factor
# fac <- gl(3,5, labels=c("groupe1","groupe2","groupe3"))
#
# Bartlett test with centring on the group medians
# bartlett.perm(yy, fac, centr="MEDIAN", nperm=999, alpha=0.05)
# ------------------------------------------------------


{
  
  fact <- as.factor(fact)
  
  normal <- shapiro.test(resid(lm(y~fact)))
  if(normal[2]<=alpha){
    cat("\n-------------------------------------------------------------------")
    cat("\nWARNING") 
    cat("\nThe residuals of the ANOVA model are not normally distributed.")
    cat("\nThis is likely to change the rate of type I error of the test")
    cat("\nof homogeneity of variances.")
    cat("\n-------------------------------------------------------------------")
    cat("\n")
  }
  
  # Trap for groups with 0 dispersion
  y.sd <- tapply(y, fact, sd)
  if(any(y.sd==0)) {
    cat("\n-------------------------------------------------------------------")
    cat("\nPROBLEM ENCOUNTERED") 
    cat("\nOne or more groups have zero variance. Please chek and correct.")
    cat("\nThe computations are meaningless if a group is made of observations")
    cat("\nthat all have the same value.")
    cat("\n-------------------------------------------------------------------")
    cat("\n")
    stop
  }
  
  
  CENTRE <- c("MEDIAN", "MEAN")
  centr <- match.arg(centr, CENTRE)
  
  
  # Within-group centring of data
  
  if(centr == "MEDIAN"){
    meds <- tapply(y, fact, median, na.rm=TRUE)
    y.c <- y - meds[fact]
  }
  else{
    means <- tapply(y, fact, mean, na.rm=TRUE)
    y.c <- y - means[fact]
  }
  
  
  bart <- bartlett.test(y.c,fact)
  
  # Permutation tests
  
  cat("\nPermutations running...")
  cat("\n")
  
  compt.perm <- 1
  
  for(i in 1:nperm) {
    
    yprime <- sample(y.c)
    
    bart.perm <- bartlett.test(yprime,fact)
    if(bart.perm[[1]] >= bart[[1]]){
      compt.perm=compt.perm+1}
    
  }
  
  # Bootstrap tests
  # Difference with permutation test: resampling is done with replacement
  
  cat("\nBootstrap running...")
  cat("\n")
  cat("\n")
  
  compt.boot <- 1
  
  for(i in 1:nperm) {
    
    yboot <- sample(y.c, replace=TRUE)
    
    bart.boot <- bartlett.test(yboot,fact)
    if(bart.boot[[1]] >= bart[[1]]){
      compt.boot=compt.boot+1}
    
  }
  
  
  Result <- matrix(0,1,4)
  colnames(Result) <- c("Statistic", "Param.prob", "Permut.prob", "Bootstrap.prob")
  rownames(Result) <- "Bartlett" 
  
  Result[1,1] <- round(bart[[1]],4)
  Result[1,2] <- round(bart[[3]],4)
  Result[1,3] <- compt.perm/(nperm+1)
  Result[1,4] <- compt.boot/(nperm+1)
  
  Result
  
}


### Function to perform (unpaired) Kruskal-Wallis or paired Wilcoxon test with 
### post-hoc multiple comparisons and boxplots with letters 
### (using agricolae::kruskal)

# Arguments:
# Y = a numeric response variable
# X = a factor (groups, treatments...) or a qualitative variable that can be 
# converted to a factor
# p.adj = correction of p-values for multiple comparisons 
# (default=none, bonferroni, holm...)
# paired = if TRUE and if 2 groups only, then Wilcoxon paired test is performed
# (objects are supposed to be ordered twice the same way in the data frame)

# Value:
# A summary table ($comparison) is printed and boxplots are drawn with result
# of the test ($p.value) and, if it is significant, of post-hoc tests as letters
# (in decreasing order)

# # Example:
# library(agricolae)
# source("boxplerk.R")
# library(stats)
# data(InsectSprays)
# boxplerk(
#   Y = InsectSprays$count,
#   X = InsectSprays$spray,
#   ylab = "count",
#   xlab = "spray",
#   bcol = "bisque",
#   p.adj = "holm"
# )

# License: GPL-2
# Author:  Francois Gillet
#          2021-08-22 (adapted to agricolae 1.3-5)


boxplerk <-
  function(Y,
           X,
           main = NULL,
           xlab = "X factor",
           ylab = "Y value",
           bcol = "bisque",
           p.adj = "none",
           cexy = 1,
           varwidth = TRUE,
           las = 1,
           paired = FALSE) {
    
    if (!is.factor(X)) {
      X <- as.factor(X)
    }
    aa <- levels(X)
    
    tt1 <- matrix(nrow = length(aa), ncol = 7)
    for (i in 1:length(aa)) {
      temp <- Y[X == aa[i]]
      tt1[i, 1] <- mean(temp, na.rm = TRUE)
      tt1[i, 2] <- sd(temp, na.rm = TRUE) / sqrt(length(temp))
      tt1[i, 3] <- sd(temp, na.rm = TRUE)
      tt1[i, 4] <- min(temp, na.rm = TRUE)
      tt1[i, 5] <- max(temp, na.rm = TRUE)
      tt1[i, 6] <- median(temp, na.rm = TRUE)
      tt1[i, 7] <- length(temp)
    }
    tt1 <- as.data.frame(tt1)
    row.names(tt1) <- aa
    colnames(tt1) <- c("mean", "se", "sd", "min", "max", "median", "n")
    
    boxplot(
      Y ~ X,
      main = main,
      xlab = xlab,
      ylab = ylab,
      las = las,
      col = bcol,
      cex.axis = cexy,
      cex.lab = cexy,
      varwidth = varwidth)
    
    require(agricolae)
    comp <- kruskal(Y, X, p.adj = p.adj)
    gror <- comp$groups[aa, ]
    tt1$rank <- gror$Y
    tt1$group <- gror$groups
    sig <- "ns"
    
    if (paired == TRUE & length(aa) == 2) {
      coms <- wilcox.test(Y ~ X, paired = TRUE)
      pp <- coms$p.value
      tt1$group <- rep("a", 2)
      if (pp <= 0.05) 
        tt1$group[which(tt1$rank == min(tt1$rank))] <- "b"
    }
    else {
      pp <- comp$statistics$p.chisq
    }
    
    if (pp <= 0.1)
      sig <- "."
    if (pp <= 0.05)
      sig <- "*"
    if (pp <= 0.01)
      sig <- "**"
    if (pp <= 0.001)
      sig <- "***"
    if (pp <= 0.0001)
      sig <- "****"
    mtext(
      sig,
      side = 3,
      line = 0.5,
      adj = 0,
      cex = 2,
      font = 1
    )
    
    if (pp <= 0.1)
      mtext(
        tt1$group,
        side = 3,
        at = c(1:length(aa)),
        line = 0.5,
        cex = 1,
        font = 4
      )
    
    list(comparison = tt1, p.value = pp, p.adjust = p.adj)
  }


### Function to perform (unpaired) ANOVA or paired t-test with post-hoc
### multiple comparisons and boxplots with letters (using agricolae::LSD.test)

# Arguments:
# Y = a numeric response variable
# X = a factor (groups, treatments...) or a qualitative variable that can be 
# converted to a factor
# p.adj = correction of p-values for multiple comparisons
# (default=none, bonferroni, holm...)
# paired = if TRUE and if 2 groups only, then paired t-test is performed
# (objects are supposed to be ordered twice the same way in the data frame)

# Value:
# A summary table ($comparison) is printed and boxplots are drawn with result
# of the test ($p.value) and, if it is significant, of post-hoc tests as letters
# (in decreasing order)

# # Example:
# library(agricolae)
# data(sweetpotato)
# # Check ANOVA assumptions
# shapiro.test(resid(aov(sweetpotato$yield ~ sweetpotato$virus)))
# bartlett.test(sweetpotato$yield, sweetpotato$virus)
# source("boxplert.R")
# boxplert(
#   Y = sweetpotato$yield,
#   X = sweetpotato$virus,
#   ylab = "yield",
#   xlab = "virus",
#   bcol = "orange",
#   p.adj = "holm"
# )

# License: GPL-2
# Author:  Francois Gillet
#          2021-08-22 (adapted to agricolae 1.3-5)


boxplert <-
  function(Y,
           X,
           main = NULL,
           xlab = "X factor",
           ylab = "Y value",
           bcol = "bisque",
           p.adj = "none",
           cexy = 1,
           varwidth = TRUE,
           las = 1,
           paired = FALSE) {
    
    if (!is.factor(X)) {
      X <- as.factor(X)
    }
    aa <- levels(X)
    
    tt1 <- matrix(nrow = length(aa), ncol = 6)
    for (i in 1:length(aa)) {
      temp <- Y[X == aa[i]]
      tt1[i, 1] <- mean(temp, na.rm = TRUE)
      tt1[i, 2] <- sd(temp, na.rm = TRUE) / sqrt(length(temp))
      tt1[i, 3] <- sd(temp, na.rm = TRUE)
      tt1[i, 4] <- min(temp, na.rm = TRUE)
      tt1[i, 5] <- max(temp, na.rm = TRUE)
      tt1[i, 6] <- length(temp)
    }
    tt1 <- as.data.frame(tt1)
    row.names(tt1) <- aa
    colnames(tt1) <- c("mean", "se", "sd", "min", "max", "n")
    
    boxplot(
      Y ~ X,
      main = main,
      xlab = xlab,
      ylab = ylab,
      las = las,
      col = bcol,
      cex.axis = cexy,
      cex.lab = cexy,
      varwidth = varwidth
    )
    
    require(agricolae)
    sig <- "ns"
    model <- aov(Y ~ X)
    
    if (paired == TRUE & length(aa) == 2) {
      coms <- t.test(Y ~ X, paired = TRUE)
      pp <- coms$p.value
      tt1$group <- rep("a", 2)
      if (pp <= 0.05) 
        tt1$group[which(tt1$mean == min(tt1$mean))] <- "b"
    }
    else {
      pp <- anova(model)$Pr[1]
      comp <- LSD.test(model,
                       "X",
                       alpha = 0.05,
                       p.adj = p.adj,
                       group = TRUE)
      gror <- comp$groups[aa, ]
      tt1$group <- gror$groups
    }
    
    if (pp <= 0.1)
      sig <- "."
    if (pp <= 0.05)
      sig <- "*"
    if (pp <= 0.01)
      sig <- "**"
    if (pp <= 0.001)
      sig <- "***"
    if (pp <= 0.0001)
      sig <- "****"
    mtext(
      sig,
      side = 3,
      line = 0.5,
      adj = 0,
      cex = 2,
      font = 1
    )
    
    if (pp <= 0.1) {
      mtext(
        tt1$group,
        side = 3,
        at = c(1:length(aa)),
        line = 0.5,
        cex = 1,
        font = 4
      )
    }
    
    list(comparison = tt1, p.value = pp, p.adjust = p.adj)
  }


# Function to compute a binary dissimilarity matrix from clusters
grpdist <- function(X){
  require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  distgr
}
