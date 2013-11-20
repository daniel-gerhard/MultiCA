mbca <- function(mat, group){
  mat <- as.matrix(mat)
  nr <- nrow(mat)
  nc <- ncol(mat)
  if (is.null(rownames(mat))) rownames(mat) <- 1:nr
  if (is.null(colnames(mat))) colnames(mat) <- 1:nc
  row.w <-  rep(1, nr)/nr
  col.w <- rep(1, nc)
  if (any(is.na(mat))) stop("na entries in mat")
  # centering and scaling analog to ade4
  f1 <- function(v) sum(v * row.w)/sum(row.w)
  f2 <- function(v) sqrt(sum(v * v * row.w)/sum(row.w))
  center <- apply(mat, 2, f1)
  mat <- sweep(mat, 2, center)
  norm <- apply(mat, 2, f2)
  norm[norm < 1e-08] <- 1
  mat <- sweep(mat, 2, norm, "/")
  
  smat <- mat * sqrt(row.w)
  smat <- sweep(smat, 2, sqrt(col.w), "*")
  storage.mode(smat) <- "double"
  if (nr < nc){
    pevlist <- .Call("tcrossprodeigen", smat, package="MultiCA")
  } else {
    pevlist <- .Call("crossprodeigen", smat, package="MultiCA")
  }
  oe <- order(pevlist$eigenvalues, decreasing=TRUE)
  pevlist$eigenvalues <- pevlist$eigenvalues[oe]
  peig <- pevlist$eigenvalues
  prank <- sum((peig/peig[1]) > 1e-07)
  
  # group PCA
  cla.w <- tapply(row.w, group, sum)
  mean.w <- function(x, w, group, cla.w) {
    z <- x * w
    z <- tapply(z, group, sum)/cla.w
    return(z)
  }
  tab <- apply(mat, 2, mean.w, w = row.w, group = group, cla.w = cla.w)
  rownames(tab) <- levels(group)  
  
  tnr <- nrow(tab)  
  stab <- tab * sqrt(as.vector(cla.w))
  stab <- sweep(stab, 2, sqrt(col.w), "*")
  storage.mode(stab) <- "double"
  if (tnr < nc){
    evlist <- .Call("tcrossprodeigen", stab, package="MultiCA")
  } else {
    evlist <- .Call("crossprodeigen", stab, package="MultiCA")
  }
  oe <- order(evlist$eigenvalues, decreasing=TRUE)
  evlist$eigenvalues <- evlist$eigenvalues[oe]
  evlist$eigenvectors <- evlist$eigenvectors[,oe]
  
  eig <- evlist$eigenvalues
  rank <- sum((eig/eig[1]) > 1e-07)
  
  ratio <- sum(eig[1:rank])/sum(peig[1:prank])  
  dval <- sqrt(eig[1:rank])    
  
  if (tnr >= nc) {
    c1 <- evlist$eigenvectors[,1:rank] * 1/sqrt(col.w)
    auxi2 <- sweep(tab, 2, col.w, "*")
    li <- auxi2 %*% c1
    colnames(c1) <- paste("CS", (1:ncol(c1)), sep = "")    
    rownames(c1) <- make.names(colnames(tab), unique = TRUE)
    colnames(li) <- paste("Axis", 1:ncol(li), sep = "")
    rownames(li) <- rownames(tab)
    co <- sweep(c1, 2, dval, "*")
    colnames(co) <- paste("Comp", (1:ncol(co)), sep = "")
    l1 <- sweep(li, 2, dval, "/")
    colnames(l1) <- paste("RS", (1:ncol(l1)), sep = "")
  } else {
    l1 <- evlist$eigenvectors[,1:rank] * 1/sqrt(as.vector(cla.w))
    auxi2 <- t(sweep(tab, 1, cla.w, "*"))
    co <- data.frame(auxi2 %*% l1)
    colnames(l1) <- paste("RS", (1:ncol(l1)), sep = "")
    rownames(l1) <- rownames(tab)
    colnames(co) <- paste("Comp", (1:ncol(co)), sep = "")
    rownames(co) <- make.names(colnames(tab), unique = TRUE)
    li <- sweep(l1, 2, dval, "*")
    colnames(li) <- paste("Axis", (1:ncol(li)), sep = "")
    c1 <- sweep(co, 2, dval, "/")
    colnames(c1) <- paste("CS", (1:ncol(c1)), sep = "")
  }
  
  ls <- mat %*% as.matrix((c1 * unlist(col.w)))
  rownames(ls) <- rownames(mat)
  colnames(ls) <- colnames(c1)
  
  out <- list()
  out$eig <- eig
  out$pcaeig <- peig
  out$lw <- cla.w
  out$plw <- row.w
  out$cw <- col.w
  out$tab <- tab
  out$mat <- mat
  out$group <- group
  out$li <- li
  out$l1 <- l1
  out$co <- co
  out$c1 <- c1
  out$ls <- ls
  class(out) <- "mbca"
  return(out)  
}