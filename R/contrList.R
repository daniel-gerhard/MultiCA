contrList <- function(n, type=c("Dunnett", "Tukey", "AVE"), base=1){
  type <- match.arg(type)
  k <- length(n)
  knames <- names(n)  
  if (is.null(knames)) knames <- as.character(1:k)
  switch(type, Dunnett = {
    v0 <- numeric(length=k)
    cvec <- v0
    cvec[base] <- 1
    CL <- lapply(c(1:k)[-base], function(i){
      nvec <- v0
      nvec[i] <- 1
      m <- rbind(cvec, nvec)
      rownames(m) <- knames[c(base, i)]
      colnames(m) <- knames
      return(m)
    }) 
    names(CL) <- paste(knames[base], "vs.", knames[c(1:k)[-base]])
  }, Tukey = {
    combi <- combn(1:k, 2)
    v0 <- numeric(length=k)
    CL <- lapply(1:ncol(combi), function(i){
      vec1 <- vec2 <- v0
      vec1[combi[1,i]] <- 1
      vec2[combi[2,i]] <- 1
      m <- rbind(vec1, vec2)
      rownames(m) <- knames[combi[,i]]
      colnames(m) <- knames
      return(m)
    })
    names(CL) <- apply(combi, 2, function(x) paste(knames[x], collapse=" vs. "))
  }, AVE = {
    CL <- lapply(1:k, function(i){
      vec1 <- numeric(length=k)
      vec1[i] <- 1
      vec2 <- n/sum(n[-i])
      vec2[i] <- 0
      m <- rbind(vec1, vec2)
      rownames(m) <- c(knames[i], paste("AVE", paste(knames[-i], collapse=",")))
      colnames(m) <- knames
      return(m)
    })
    names(CL) <- paste("C", 1:k, sep="")
  })  
  class(CL) <- "contrList"
  return(CL)
}