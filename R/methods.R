print.bcatest <- function(x, ...){
  cat("\tbca permutation test\n\n")
  cat("statistic:\n")
  print(x$statistic)
  cat("\np-values:\n")
  print(x$pvalues)  
}

print.mbca <- function(x, ...){
  pcaeig <- x$pcaeig
  pcarank <- sum((pcaeig/pcaeig[1]) > 1e-07)
  inertot <- sum(pcaeig[1:pcarank])
  eig <- x$eig
  rank <- sum((eig/eig[1]) > 1e-07)
  gdist <- sum(eig[1:rank])  
  cat("\tBetween Component Analysis\n")
  cat("\nEigenvalues:\n")
  print(eig[1:rank])  
  cat("\nInertia percentage:\n")
  print(gdist/inertot)  
}

plot.mbca <- function(x, ...){
  def.par <- par(no.readonly = TRUE)
  on.exit(par(def.par))
  layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3), respect = TRUE)
  s.arrow(x$c1*2)
  scatterutil.eigen(x$eig, nf=2)
  s.class(x$ls, x$group, col=1:length(levels(x$group)))
}
