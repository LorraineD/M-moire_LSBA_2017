RMSE <- function(Ychap,Y,k=1){
  PRESS <- sum((Y - Ychap)^2)
  if (is.matrix(Y)) {n <- nrow(Y) ; m <- ncol(Y)}   else {n <- 1 ; m <- length(Y)}
  MSE1 <- PRESS / (n*m)
  MSE0 <- sum(scale(Y,scale=F)^2) / (n*m)
  MSE <- MSE1 + 0.03 * k * MSE0
  RMSE <- sqrt(MSE)
  return(RMSE)
  
}

VE <- function(Ychap,Y,k=1){
  PRESS <- sum((Y - Ychap)^2)
  if (is.matrix(Y)) {n <- nrow(Y) ; m <- ncol(Y)}   else {n <- 1 ; m <- length(Y)}
  MSE1 <- PRESS / (n*m)
  MSE0 <- sum(scale(Y,scale=F)^2) / (n*m)
  MSE <- MSE1 + 0.03 * k * MSE0
  VE <- 100*(MSE0 - MSE)/MSE 
  return(VE)
  
}



Q2 <- function(Ychap,Y){
  if(is.vector(Y)){
    PRESS <- sum((Y - Ychap)^2)
    TSS <- sum(scale(Y,scale=F)^2)
    Q2 <- 1 - PRESS/TSS
    return(Q2)
  }
  
  if(is.matrix(Y)){
    ncol <- ncol(Y)
    PRESS <- sum((Y - Ychap)^2)
    TSS <-  sum(scale(Y,scale=F)^2)
    Q2 <- 1 - PRESS/TSS
    
    Q2j <- numeric(ncol)
    for(i in 1:ncol){
      PRESS <- sum((Y[,i] - Ychap[,i])^2)
      TSS <- sum(scale(Y[,i],scale=F)^2)
      Q2j[i] <- 1 - PRESS/TSS    }
   return(list(Q2,Q2j))
  }
}

normD <- function(x, D=D) {
  if(is.vector(x)) sqrt(t(x) %*% D %*% x) 
  if(is.matrix(x)) tr(t(x)%*%D%*%x)^(1/2)
  else "ERROR"
}


PRESS_GOMCIA <- function(Y, Ychap){
  EY <- Y - Ychap
  sum(apply(EY,1,function(x) sqrt( sum(x^2)))) * 1 / nrow(Y)
}



R2 <- function(res,A = NULL,gen=1, bloc="Y"){
  if (is.null(A)) A <- res$A
  nb <- length(res$Xini)
  
  if (bloc %in% c("Y","ALL")) {
    Y <- res$Yini[[1]]
    R.y <- 0 
    for(k in 1:A){
      if(gen==1) T_k <- res$Tk[[k]]
      if(gen==3) {
        T_k <- 0
        for(i in 1:nb) {T_k  <- T_k + res$MU[i,k] * res$Tk[[k]][,i]}
        T_k <- matrix(T_k,ncol=1)
        }
     
      val <- normD( T_k %*% ginv(t(T_k) %*% T_k) %*% t(T_k) %*% Y ,D = res$D)^2 / normD(Y,res$D)^2
      R.y <- R.y + val}}      
  
  if (bloc %in% c("X","ALL")) {
    R.x <- numeric(nb)
    for(k in 1:A){
      if(gen==1) T_k <- res$Tk[[k]]
      if(gen==3) {
        T_k <- 0
        for(j in 1:nb) {T_k  <- T_k + res$MU[j,k] * res$Tk[[k]][,j]}
        T_k <- matrix(T_k,ncol=1)
      }
    
    for(i in 1:nb){
        X <- res$Xini[[i]]
        R.x[i] <- 0 
      val <- normD( T_k %*% ginv(t(T_k) %*% T_k) %*% t(T_k) %*% X ,D = res$D)^2 / normD(X,res$D)^2
      R.x[i] <- R.x[i] + val}}}
  
  if (bloc == "Y") return(R.y)
  if (bloc == "X") return(R.x)
  if (bloc == "ALL") return(list(R.y,R.x))
  
  }  
      
      

tr <- function (m){
  total_sum <- 0
  if(is.matrix(m))  {
    row_count <- nrow(m)
    col_count <- ncol(m)
    if(row_count == col_count)    {
      total_sum <-sum(diag(m))
      total_sum    }
    else    {      message ('Matrix is not square')    }  }
  else  {    message( 'Object is not a matrix')  }}



VE_PCA <- function(x,tg,K){
  nb <- length(x)
  VE <- matrix(0,nb,K)
  for(b in 1:nb){
  XD <- x[[b]]
  TOT <- sum(diag(cov(NCI[[b]])))
  for(i in 1:K){
    a <- t(XD) %*% tg[,i] %*% ginv(t(tg[,i]) %*% tg[,i])
    tt <- tg[,i] %*% t(a)
    B <- ginv(t(tg[,i]) %*% tg[,i]) %*% t(tg[,i]) %*%XD
    XD <- XD - tg[,i]%*%B
    VE[b,i] <- VE[b,i] + norm(tt,"2")/TOT
  }
  }
  VEB <- VE
  for(i in 2:K){
    VEB[,i] <- rowSums(VE[,1:i])
  }
  return(VEB)
  }



  
