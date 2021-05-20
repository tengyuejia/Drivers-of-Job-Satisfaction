cormatrix <- function(x){
  cormatrix <- rcorr(as.matrix(x))
  cormatrix.r <- cormatrix[[1]]
  cormatrix.p <- cormatrix[[3]]
  cormatrix.N <- cormatrix[[2]]
  cormatrix.rANDp<-matrix(NA,nrow=nrow(cormatrix.r),ncol=ncol(cormatrix.r))
  colnames(cormatrix.rANDp)<-colnames(cormatrix.r)
  for (i in seq_len(nrow(cormatrix.r))) {
    for (j in seq_len(ncol(cormatrix.r))){
      if (
        is.na(cormatrix.p[i,j])){
        a <- c(" ")
      }else if (round(cormatrix.p[i,j],digits = 2) <= 0.001){
        a <- c("***")
      }else if (round(cormatrix.p[i,j],digits = 2) <= 0.01 && round(cormatrix.p[i,j],digits = 2) > 0.001){
        a<-c("**")
      }else if (round(cormatrix.p[i,j],digits = 2) <= 0.05 && round(cormatrix.p[i,j],digits = 2) > 0.01){
        a<-c("*")
      }else if (round(cormatrix.p[i,j],digits = 2) > 0.05){
        a<-c(" ")
      }
      cormatrix.rANDp[i,j]<-paste0(round(cormatrix.r[i,j], digits = 2), a , collapse = "")
    }
  }
  cormatrix.rANDp <- as.data.frame(cormatrix.rANDp)
  cormatrix.rANDp
}