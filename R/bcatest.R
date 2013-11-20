bcatest <- function(x, Klist, nf=NULL, nperm=999){
  X <- x$mat
  plw <- x$plw
  lw <- x$lw
  group <- x$group
  
  kinert <- function(tab, group, plw, lw, Klist, nf){
    mm <- apply(tab, 2, function(x, plw, lw, group){ 
      tapply(x * plw, group, sum)/lw
    }, plw=plw, lw=lw, group=group) 
    eig <- sapply(Klist, function(K, mm, plw, lw, nf){
      wm <- K %*% mm
      wm <- scale(wm, scale=FALSE)
      mtab <- wm * sqrt(as.vector(K %*% as.vector(lw)))
      if (nrow(mtab) < ncol(mtab)){
        evlist <- .Call("tcrossprodeigen", mtab, package="MultiCA")
      } else {
        evlist <- .Call("crossprodeigen", mtab, package="MultiCA")
      }
      eig <- sort(evlist$eigenvalues, decreasing=TRUE)
      rank <- sum((eig/eig[1]) > 1e-07)
      if (!is.null(nf)){
        if (nf < rank) rank <- nf
      }
      eig <- eig[1:rank]
      sum(eig)
    }, mm=mm, plw=plw, lw=lw, nf=nf)
    return(eig)
  }  
  
  ppi <- replicate(nperm, {
    pid <- sample(1:nrow(X))
    kinert(X[pid,], group, plw, x$lw, Klist, nf)
  })  
  
  pcaeig <- x$pcaeig
  pcarank <- sum((pcaeig/pcaeig[1]) > 1e-07)
  if (!is.null(nf)){
    if (nf < pcarank) pcarank <- nf
  }
  inertot <- sum(pcaeig[1:pcarank])
  
  oratio <- kinert(X, group, plw, x$lw, Klist, nf)/inertot
  permat <- ppi/inertot
  pvalues <- apply(apply(permat, 2, function(x, oratio) max(x) > oratio, oratio=oratio), 1, mean)
  if (!is.null(names(Klist))) rownames(permat) <- names(oratio) <- names(pvalues) <- names(Klist)
  
  out <- list()
  out$statistic <- oratio
  out$pvalues <- pvalues
  out$permutations <- t(permat)
  out$nf <- nf
  class(out) <- "bcatest"
  return(out)
}