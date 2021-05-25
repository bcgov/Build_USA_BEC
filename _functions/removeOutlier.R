removeOutlier <- function(dat, alpha,numIDvars){
  out <- foreach(curr = unique(as.character(dat$BGC)), .combine = rbind) %do% {
    temp <- dat[dat$BGC == curr,]
    md <- tryCatch(mahalanobis(temp[,-c(1:numIDvars)],
                               center = colMeans(temp[,-c(1:numIDvars)]),
                               cov = cov(temp[,-c(1:numIDvars)])), error = function(e) e)
    if(!inherits(md,"error")){
      ctf <- qchisq(1-alpha, df = ncol(temp)-1)
      outl <- which(md > ctf)
      cat("Removing", length(outl), "outliers from",curr, "; ")
      if(length(outl) > 0){
        temp <- temp[-outl,]
      }
    }
    temp
    
  }
  return(out)
}