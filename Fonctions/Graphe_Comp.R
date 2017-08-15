Plot.TK <- function(res, comp=NULL, bloc.x=NULL, bloc.y=NULL, Y.Class=NULL, main=NULL, ylab=NULL, xlab=NULL){
  
  if(is.null(comp)) comp <- c(1,2)

  if ("gen"%in%names(res)) {
    nb.x <- ncol(res$Tk[[1]])
    nb.y <- ncol(res$Ul[[1]])
    if(is.null(bloc.x)) bloc.t <- seq(1,nb.x)
    if(is.null(bloc.y)) bloc.u <- seq(1,nb.y)
    TK.1 <- data.frame(res$Tk[[comp[1]]][,bloc.t])
    TK.2 <- data.frame(res$Tk[[comp[2]]][,bloc.t])
    UK.1 <- data.frame(res$Ul[[comp[1]]][,bloc.u])
    UK.2 <- data.frame(res$Ul[[comp[2]]][,bloc.u])
    ID <- rownames(TK.1)
    if(is.null(names(res$Xinitial))) {NameX <- paste0("BlocX",seq(1,length(res$Xinitial)))}
      else NameX <- names(res$Xinitial)
    if(is.null(names(res$Yinitial))) {NameY <- paste0("BlocY",seq(1,length(res$Yinitial)))}
      else NameY <- names(res$Yinitial)
    
    MAT <- NULL
    for(i in 1:nb.x) MAT <- rbind(MAT,data.frame(D1 = TK.1[,i], D2 = TK.2[,i], BLOC=NameX[i]))
    for(i in 1:nb.y) MAT <- rbind(MAT,data.frame(D1 = UK.1[,i], D2 = UK.2[,i], BLOC=NameY[i]))
    
    nb <- NB <- nb.x + nb.y 
    YEND <- XEND <- NULL 
    for(i in 2:nb.x) {
      XEND <- c(XEND,TK.1[,i])
      YEND <- c(YEND,TK.2[,i])}
    for(i in 1:nb.y) {
      XEND <- c(XEND,UK.1[,i])
      YEND <- c(YEND,UK.2[,i])}
    
    SEG <- data.frame(X=rep(TK.1[,1],(NB-1)),Y=rep(TK.2[,1],(NB-1)),XEND=XEND,YEND=YEND)
  }
  
  
  if ("AVE"%in%names(res)) {
    nb <- length(res$X)
    if(is.null(bloc.x)) bloc.t <- seq(1,nb) else bloc.t <- bloc
    nb <- length(bloc.t)
    if(is.null(names(res$X))) NameX <- paste0("BlocX",seq(1,length(res$Y))) else NameX <- names(res$X)
    XEND <- YEND <- MAT <- NULL
    for(k in 1:nb) {
      i <- bloc.t[k]
      MAT <- rbind(MAT,data.frame(D1 = res$variates[[i]][,comp[1]], D2 =  res$variates[[i]][,comp[2]], BLOC=NameX[i]))}
    for(k in 2:nb) {
      i <- bloc.t[k]
      XEND <- c(XEND,res$variates[[i]][,comp[1]])
      YEND <- c(YEND,res$variates[[i]][,comp[2]])
    }
    SEG <- data.frame(X=rep(res$variates[[bloc.t[1]]][,comp[1]],(nb-1)),Y=rep(res$variates[[bloc.t[1]]][,comp[2]],(nb-1)),XEND=XEND,YEND=YEND)
      }
  
  if("tb"%in%names(res)){
    nb <- length(res$tb)
    if(is.null(bloc.x)) bloc.t <- seq(1,nb) else bloc.t <- bloc
    if(is.null(names(res$tb))) NameX <- paste0("BlocX",seq(1,length(res$Y))) else NameX <- names(res$tb)
    XEND <- YEND <- MAT <- NULL
    for(k in 1:nb) {
      i <- bloc.t[k]
      MAT <- rbind(MAT,data.frame(D1 = res$tb[[i]][,comp[1]], D2 =  res$tb[[i]][,comp[2]], BLOC=NameX[i]))
      }
    for(k in 2:nb) {
      i <- bloc.t[k]
      XEND <- c(XEND,res$tb[[i]][,comp[1]])
      YEND <- c(YEND,res$tb[[i]][,comp[2]])
    }
    SEG <- data.frame(X=rep(res$tb[[bloc.t[1]]][,comp[1]],(nb-1)),Y=rep(res$tb[[bloc.t[1]]][,comp[2]],(nb-1)),XEND=XEND,YEND=YEND)
  }
    

  if(is.null(main)) main <- paste("Graphique des scores pour les composantes",comp[1],"et",comp[2])
  if(is.null(xlab)) xlab <- paste("Composante",comp[1])
  if(is.null(ylab)) ylab <- paste("Composante",comp[2])
  
  ggplot()  + geom_segment(aes(X,Y,xend=XEND,yend=YEND,col=rep(ClassY,(nb-1))),data=SEG) + geom_point(data=MAT,aes(D1,D2,pch=BLOC,col=rep(ClassY,nb))) + guides(color=guide_legend("Type de cancers:")) + labs(title=main) + ylab(ylab) + xlab(xlab) + theme_bw() + theme(text = element_text(size=18))
    
  }
  
Plot.W <- function(res, comp=NULL, bloc.x=NULL, bloc.y=NULL, main=NULL, ylab=NULL, xlab=NULL, type=NULL){
  
  if (length(type%in%c("cor","cov"))==0) stop("Erreur : choisir cor ou cov pour le type")
  if (is.null(comp)) comp <- c(1,2)
  
  if ("gen"%in%names(res)) {
    nb.x <- ncol(res$Tk[[1]])
    nb.y <- ncol(res$Ul[[1]])
    if(is.null(bloc.x)) bloc.t <- seq(1,nb.x)
    if(is.null(bloc.y)) bloc.u <- seq(1,nb.y)}
    NB <- nb.x + nb.y 
  
    VCV <- list()
    for(b in 1:length(bloc.t)){
      nv <- ncol(res$Xinitial[[bloc.t[b]]])
      VCV[[b]] <- matrix(NA,ncol=2,nrow=nv)
      for(a in 1:2){
        for(j in 1:nv){
        VCV[[b]][,a] <- cor(res$Tk[[comp[a]]][,bloc.t[b]],res$Xinitial[[bloc.t[[b]]]])[1,]
        }}}
        
    for(b in 1:length(bloc.u)){
      nv <- ncol(res$Yinitial[[bloc.u[b]]])
      VCV[[(length(bloc.t) + b)]] <- matrix(NA,ncol=2,nrow=nv)
      for(a in 1:2){
        for(j in 1:nv){
          if(type=="cor") {VCV[[(length(bloc.t) + b)]][,a] <- cor(res$Ul[[comp[a]]][,bloc.u[b]],res$Yinitial[[bloc.u[[b]]]])[1,]}
          if(type=="cov") {VCV[[(length(bloc.t) + b)]][,a] <- cov(res$Ul[[comp[a]]][,bloc.u[b]],res$Yinitial[[bloc.u[[b]]]])[1,]}
        }}}
    MAT <- NULL
    Name <- c(names(res$Xinitial[bloc.t]),names(res$Yinitial[bloc.u]))
    for(i in 1:length(VCV)) MAT <- rbind(MAT,data.frame(C1=VCV[[i]][,1],C2=VCV[[i]][,2],Bloc = Name[i]))
    
    
    
    if(is.null(main)) main <- paste("Cercle de corrélation pour les composantes",comp[1],"et",comp[2])
    if(is.null(xlab)) xlab <- paste("Composante",comp[1])
    if(is.null(ylab)) ylab <- paste("Composante",comp[2])
        
   ggplot(MAT,aes(C1,C2,col=Bloc))  + geom_point() + guides(color=guide_legend("Bloc")) + labs(title=main) + ylab(ylab) + xlab(xlab)
      return(invisible(list(Graphe,MAT)))
      }
      
      
Var.Covar <- function(res, comp=NULL, bloc.x=NULL, bloc.y=NULL, type=NULL){
  
  if (length(type%in%c("cor","cov"))==0) stop("Erreur : choisir cor ou cov pour le type")

  if ("gen"%in%names(res)) {
    if (is.null(comp)) comp <- seq(1,res$A)
    
    nb.x <- ncol(res$Tk[[1]])
    nb.y <- ncol(res$Ul[[1]])
    if(is.null(bloc.x)) bloc.t <- seq(1,nb.x)
    if(is.null(bloc.y)) bloc.u <- seq(1,nb.y)
      NB <- nb.x + nb.y 
  Tk.Super <- data.frame(lapply(res$Tk,rowSums))
  colnames(Tk.Super) <- paste("Composante",comp)
  VCV <- list()
  for(b in 1:length(bloc.t)){
    nv <- ncol(res$Xinitial[[bloc.t[b]]])
    VCV[[b]] <- matrix(NA,ncol=length(comp),nrow=nv)
    for(a in 1:length(comp)){
      for(j in 1:nv){
        VCV[[b]][,a] <- cor(Tk.Super[,comp[a]],res$Xinitial[[bloc.t[[b]]]])[1,]
      }}}
  
  for(b in 1:length(bloc.u)){
    nv <- ncol(res$Yinitial[[bloc.u[b]]])
    VCV[[(length(bloc.t) + b)]] <- matrix(NA,ncol=length(comp),nrow=nv)
    for(a in 1:length(comp)){
        if(type=="cor") {VCV[[(length(bloc.t) + b)]][,a] <- cor(res$Ul[[comp[a]]][,bloc.u[b]],res$Yinitial[[bloc.u[[b]]]])[1,]}
        if(type=="cov") {VCV[[(length(bloc.t) + b)]][,a] <- cov(res$Ul[[comp[a]]][,bloc.u[b]],res$Yinitial[[bloc.u[[b]]]])[1,]}
      }}

  Name <- c(names(res$Xinitial[bloc.t]),names(res$Yinitial[bloc.u]))}
  
  if ("AVE"%in%names(res)) {
    if (is.null(comp)) comp <-max(res$ncomp)
    
    nb <- length(res$X)
    if(is.null(bloc.x)) bloc.t <- seq(1,nb) else bloc.t <- bloc
    nb <- length(bloc.t)
    if(is.null(names(res$X))) NameX <- paste0("BlocX",seq(1,length(res$Y))) else NameX <- names(res$X)
    
    Tk.Super <- matrix(0, ncol=comp,nrow=nrow(res$X[[1]]))
    for(i in 1:A) {
      for(j in 1:nb) Tk.Super[,i] <- Tk.Super[,i] + res$variates[[j]][,i] 
          }
    colnames(Tk.Super) <- paste("Composante",seq(1,A))
    VCV <- list()
    for(b in 1:length(bloc.t)){
    nv <- ncol(res$X[[bloc.t[b]]])
    VCV[[b]] <- matrix(NA,ncol=comp,nrow=nv)
    for(a in 1:A){
      if(type=="cor")  VCV[[b]][,a] <- cor(Tk.Super[,a],res$X[[bloc.t[b]]])[1,]
      if(type=="cov")  VCV[[b]][,a] <- cov(Tk.Super[,a],res$X[[bloc.t[b]]])[1,]
    }
    Name <- names(res$X)
    }
    }
    
  if("tb"%in%names(res)){
    if (is.null(comp))   comp <- ncol(res$t)
    A <- comp
    nb <- length(res$X)
    if(is.null(bloc.x)) bloc.t <- seq(1,nb) else bloc.t <- bloc
    nb <- length(bloc.t)
    if(is.null(names(res$X))) NameX <- paste0("BlocX",seq(1,nb)) else NameX <- names(res$X)
    Tk.Super <- mod$t[,1:comp]
    colnames(Tk.Super) <- paste("Composante",seq(1,comp))
    
    VCV <- list()
    for(b in 1:length(bloc.t)){
      nv <- ncol(res$X[[bloc.t[b]]])
      VCV[[b]] <- matrix(NA,ncol=comp,nrow=nv)
      for(a in 1:A){
        if(type=="cor")  VCV[[b]][,a] <- cor(Tk.Super[,a],res$X[[bloc.t[b]]])[1,]
        if(type=="cov")  VCV[[b]][,a] <- cov(Tk.Super[,a],res$X[[bloc.t[b]]])[1,]
      }
      Name <- names(res$X)
      
    }}
  
  
      MAT <- NULL
  for(i in 1:length(VCV)) MAT <- rbind(MAT,data.frame(VCV[[i]],Bloc = Name[i]))
  return(MAT)}
  
      
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


Nmin <- function(x,N){
  vec <- as.vector(x)
  svec <- sort(vec)
  target <- svec[N]
  which(x==target,arr.ind = TRUE)
}
