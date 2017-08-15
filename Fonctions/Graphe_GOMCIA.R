# Fonction représentation GOMCIA
options("encoding" = "UTF-8")
Plot.MU <- function(res,comp1=1,comp2=2,label=rownames(res$MU),main=paste("Representation de mu pour les composantes",comp1,"et",comp2),xlab=paste("Composante",comp1),ylab=paste("Composante",comp2),lim=NULL){
  lim <- ifelse(is.null(lim),max(apply(res$MU,2,max)) + 0.5,lim)
  xx <- data.frame(res$MU)
  ggplot(xx,aes(x=0,y=0,xend=xx[,comp1],yend=xx[,comp2],label=label)) + geom_segment() + xlim(0,lim) + ylim(0,lim)  + geom_text(aes(x=xx[,comp1],y=xx[,comp2]),check_overlap = FALSE, nudge_x = 0.05,hjust=0) + theme_bw() + labs(title=main) + ylab(ylab) + xlab(xlab) + theme(text = element_text(size=18),legend.position="bottom")
}



Plot.TK <- function(res,comp=c(1,2),bloc=1,main=NULL,xlab=paste("Composante",comp1),ylab=paste("Composante",comp2),ClassY=NULL,ellipse=FALSE){
  comp1=comp[1]
  comp2=comp[2]
  
  if(bloc=="G"){
      Tk1 <- res$Tk[[comp1]] %*% res$MU[,comp1]
      Tk2 <- res$Tk[[comp2]] %*% res$MU[,comp2]
      xx  <- data.frame(Tk1,Tk2)
      if(is.null(main)) main <- paste("Projection des observations sur les composantes",comp1,"et",comp2,"-  Global",sep=" ")
    }
  else {  
    xx <- data.frame(res$Tk[[comp1]][,bloc],res$Tk[[comp2]][,bloc]) 
    if(is.null(main))  main=paste("Representation des observations pour les composantes",comp1,"et",comp2)
    }
    
  ClassY <- as.factor(ClassY)
  if(nlevels(ClassY) > 6){
    if (ellipse) {
      ggplot(xx,aes(xx[,1],xx[,2],col=ClassY,shape=ClassY))  + scale_shape_manual(values=1:nlevels(ClassY)) + geom_point(size=2)  + labs(title=main) + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :"),shape=guide_legend("Classe des observations :")) + theme(text = element_text(size=18),legend.position="bottom") + stat_ellipse(type="norm",level=0.9) 
      }
    
    else {  ggplot(xx,aes(xx[,1],xx[,2],col=ClassY)) + geom_point(size=2) + labs(title=main) + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :")) +   theme(legend.position="bottom",plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=16) , legend.text=element_text(size=16))}}
  
  else {  
    if (ellipse) {
    ggplot(xx,aes(xx[,1],xx[,2],col=ClassY,pch=ClassY)) + geom_point(size=2) + labs(title=main) + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :"),pch=guide_legend(("Classe des observations :"))) +   theme(legend.position="bottom",plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=16) , legend.text=element_text(size=16)) + stat_ellipse() }
    else {    ggplot(xx,aes(xx[,1],xx[,2],col=ClassY,pch=ClassY)) + geom_point(size=2) + labs(title=main) + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :"),pch=guide_legend(("Classe des observations :"))) +   theme(legend.position="bottom",plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=16) , legend.text=element_text(size=16))}}
  
  
}


Plot.UK <-   function(res,comp1=1,comp2=2,bloc=1,main=paste("Projection des observations sur les composantes",comp1,"et",comp2,"- bloc Y",sep=" "),xlab=NULL,ylab=NULL,ClassY=NULL){
    x1 <- data.frame(res$Ul[[comp1]])
    x2 <- data.frame(res$Ul[[comp2]])
    xx <- data.frame(x1[,bloc],x2[,bloc])
    ggplot(xx,aes(xx[,1],xx[,2],col=ClassY)) + geom_point() + labs(title=main) + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') +theme_bw()  
  }

Plot.Var <- function(res,comp1=1,comp2=2,bloc=1,label=rownames(res$Ak[[bloc]]),ylab=paste("Composante",comp2),xlab=paste("Composante",comp1),main=ifelse(type=="cov",paste("Covariance entre les variables du bloc",bloc,"et les composantes",comp1,"et",comp2),paste("Corrélation entre les variables du bloc",bloc,"et les composantes",comp1,"et",comp2)),type="cov",seuil=0,TOP=NULL){
  blocX <- res$Xini
  nb <- length(blocX)
  CovVar1 <- numeric(ncol(blocX[[bloc]]))
  CovVar2 <- numeric(ncol(blocX[[bloc]]))
  
  if (is.null(TOP)){
  if(type=="cov"){
      for(i in 1:ncol(blocX[[bloc]])){CovVar1[i] <- cov(blocX[[bloc]][,i],res$Tk[[comp1]][,bloc])}
      for(i in 1:ncol(blocX[[bloc]])){CovVar2[i] <- cov(blocX[[bloc]][,i],res$Tk[[comp2]][,bloc])}

  xx <- data.frame(CovVar1,CovVar2)
  KEEP <- apply(abs(xx),1,function(x) any(x >= seuil))
  rownames(xx) <- colnames(blocX[[bloc]])
  xx <- xx[KEEP,]
  pl<- ggplot(xx,aes(xx[,1],xx[,2],label=as.character(label[KEEP]))) + geom_point()  + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') +theme_bw() + geom_text(aes(x=xx[,1],y=xx[,2]),check_overlap = FALSE, nudge_x = 0.01,hjust=0) + labs(title=main)
  return(pl)
  }
  
  if(type=="cor"){
    for(i in 1:ncol(blocX[[bloc]])){CovVar1[i] <- cor(blocX[[bloc]][,i],res$Tk[[comp1]][,bloc])}
    for(i in 1:ncol(blocX[[bloc]])){CovVar2[i] <- cor(blocX[[bloc]][,i],res$Tk[[comp2]][,bloc])}

    xx <- data.frame(CovVar1,CovVar2)
    KEEP <- apply(abs(xx),1,function(x) any(x >= seuil))
    rownames(xx) <- colnames(blocX[[bloc]])
    xx <- xx[KEEP,]
    pl <- ggplot(xx,aes(xx[,1],xx[,2],label=label[KEEP])) + geom_point()  + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') +theme_bw() + geom_text(aes(x=xx[,1],y=xx[,2]),check_overlap = FALSE, nudge_x = 0.01,hjust=0) + labs(title=main)
    return(pl)
    }}

  if (!is.null(TOP)) {
    

    
    if(type=="cov"){
      for(i in 1:ncol(blocX[[bloc]])){CovVar1[i] <- cov(blocX[[bloc]][,i],res$Tk[[comp1]][,bloc])}
      
      for(i in 1:ncol(blocX[[bloc]])){CovVar2[i] <- cov(blocX[[bloc]][,i],res$Tk[[comp2]][,bloc])}
      
      xx <- data.frame(CovVar1,CovVar2)
      KEEP <- apply(abs(xx),1,function(x) any(x >= seuil))
      xx <- xx[KEEP,]
      rownames(xx) <- colnames(blocX[[bloc]])
      ID_1 <- rownames(xx[order(abs(xx[,1]),decreasing = TRUE),])[1:TOP]
      ID_2 <- rownames(xx[order(abs(xx[,2]),decreasing = TRUE),])[1:TOP]
      
      pl <- ggplot(xx,aes(xx[,1],xx[,2]),label=label) + geom_point(col=ifelse(label%in%c(ID_1,ID_2),'blue','black'))  + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') +theme_bw() + labs(title=main) + geom_text_repel(aes(label=ifelse(label%in%c(ID_1,ID_2),as.character(label),'')))
      return(pl)
    }
    
    if(type=="cor"){
      for(i in 1:ncol(blocX[[bloc]])){CovVar1[i] <- cov(blocX[[bloc]][,i],res$Tk[[comp1]][,bloc])}
      
      for(i in 1:ncol(blocX[[bloc]])){CovVar2[i] <- cov(blocX[[bloc]][,i],res$Tk[[comp2]][,bloc])}
      
      xx <- data.frame(CovVar1,CovVar2)
      KEEP <- apply(abs(xx),1,function(x) any(x >= seuil))
      xx <- xx[KEEP,]
      rownames(xx) <- colnames(blocX[[bloc]])
      ID_1 <- rownames(xx[order(abs(xx[,1]),decreasing = TRUE),])[1:TOP]
      ID_2 <- rownames(xx[order(abs(xx[,2]),decreasing = TRUE),])[1:TOP]
     pl <-  ggplot(xx,aes(xx[,1],xx[,2]),label=label) + geom_point(col=ifelse(label%in%c(ID_1,ID_2),'blue','black'))  + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') +theme_bw() + labs(title=main) + geom_text_repel(aes(label=ifelse(label%in%c(ID_1,ID_2),as.character(label),'')))
      return(pl)}
  }
  }
    
    
    
    
    
    
    



Plot.Ak <- function(res, comp=1,bloc=1,main=paste("Poids des variables du bloc",bloc)){
  
  if(is.vector(comp)){
    axesx <- seq(1,nrow(res$Ak[[bloc]])) 
    xx <- NULL
    for (i in 1:length(comp)){ 
      xx <- rbind(xx,cbind(axesx,res$Ak[[bloc]][,comp[i]],comp[i]))
    }
    xx <- data.frame(xx)
    xx$V3 <- paste("Dim",xx$V3)
    
    ggplot(xx) + geom_line(aes(axesx,V2,colour=V3,linetype=V3)) +  scale_x_continuous(breaks=axesx,label=rownames(res$Ak[[bloc]])) + xlab(" ") + ylab("Poids des variables") + theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.background = element_rect(fill = "white"),legend.position="top",axis.line=element_line (colour="black")) +ylim(-max(abs(xx$V2))-0.2,max(abs(xx$V2))+0.2) + labs(title=main) + guides(color=guide_legend("Composantes:"),linetype=guide_legend("Composantes:"))
    
    }

  else{
    xx <- data.frame(seq(1,length(res$Ak[[bloc]][,comp])),res$Ak[[bloc]][,comp])
  ggplot(xx,aes(xx[,1],xx[,2])) + geom_line()  +  scale_x_continuous(breaks=xx[,1],label= rownames(xx)) + xlab(" ") + ylab("Poids des variables")+theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.background = element_rect(fill = "white"),axis.line=element_line (colour="black")) +ylim(-max(abs(xx[,2]))-0.2,max(abs(xx[,2]))+0.2) + labs(title=main)
}}



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


normD <- function(x, D) {
  if(is.vector(x)) {sqrt(t(x) %*% D %*% x)}
  else tr(t(x)%*%D%*%x)^(1/2)
}


Nmin <- function(x,N){
  vec <- as.vector(x)
  svec <- sort(vec)
  target <- svec[N]
  which(x==target,arr.ind = TRUE)
}

