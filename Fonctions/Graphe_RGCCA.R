Plot.TK <- function(res,comp=c(1,2),bloc=1,main=NULL,xlab=paste("Composante",comp1),ylab=paste("Composante",comp2),ClassY=NULL,ellipse=FALSE){
  comp1=comp[1]
  comp2=comp[2]
  nb <- length(res$X)
  if (bloc=="G") {
    
  x1 <- x2 <- numeric(nrow(res$X[[1]]))
    for(i in 1:nb){
      x1 <- x1 + res$variates[[i]][,comp1]
      x2 <- x2 + res$variates[[i]][,comp2]
    }
    xx <- data.frame(x1,x2) 
    if(is.null(main)) main <- paste("Projection des observations sur les composantes",comp1,"et",comp2,"- global",sep=" ")
    }
  
  else {
    if(is.null(main)) main <- paste("Projection des observations sur les composantes",comp1,"et",comp2,"- bloc",bloc,sep=" ")
    if(!is.numeric(bloc)) {stop("le paramètre BLOC est incorrect \n Valeurs possibles : G ou valeur numérique")}
    xx <- data.frame(res$variates[[bloc]][,comp]) }
  ClassY<- factor(ClassY)
  
  if(length(table(ClassY)) > 6){
    if (ellipse) {
      ggplot(xx,aes(xx[,1],xx[,2],col=ClassY,shape=ClassY))  + scale_shape_manual(values=1:nlevels(ClassY)) + geom_point(size=2) + labs(title=main) + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :"),shape=guide_legend("Classe des observations :")) + theme(text = element_text(size=18),legend.position="bottom") + stat_ellipse(type="norm",level=0.9)
      }
    
    else {  ggplot(xx,aes(xx[,1],xx[,2],col=ClassY,shape=ClassY))  + scale_shape_manual(values=1:nlevels(ClassY)) + geom_point(size=2) + labs(title=main) + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :"),shape=guide_legend("Classe des observations :")) +   theme(legend.position="bottom",plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=16) , legend.text=element_text(size=16))}}
  
  else {  
    if (ellipse) {
      ggplot(xx,aes(xx[,1],xx[,2],col=ClassY,shape=ClassY))  + scale_shape_manual(values=1:nlevels(ClassY)) + geom_point(size=2) + labs(title=main) + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :"),shape=guide_legend(("Classe des observations :"))) +   theme(legend.position="bottom",plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=16) , legend.text=element_text(size=16)) + stat_ellipse(type="norm",level=0.9) }
    else {   ggplot(xx,aes(xx[,1],xx[,2],col=ClassY,shape=ClassY))  + scale_shape_manual(values=1:nlevels(ClassY)) + geom_point(size=2) + labs(title=main) + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :"),shape=guide_legend(("Classe des observations :"))) +   theme(legend.position="bottom",plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=16) , legend.text=element_text(size=16))
    }}
}

Plot.Var <- function(res,comp1=1,comp2=2,bloc=1,label=colnames(res$X[[bloc]]),ylab=paste("Composante",comp2),xlab=paste("Composante",comp1),main=ifelse(type=="cov",paste("Covariance entre les variables du bloc",bloc,"et les composantes",comp1,"et",comp2),paste("Corrélation entre les variables du bloc",bloc,"et les composantes",comp1,"et",comp2)),type="cov",seuil=0,TOP=NULL){
  blocX <- mod$X
  nb <- length(blocX)
  CovVar1 <- numeric(ncol(blocX[[bloc]]))
  CovVar2 <- numeric(ncol(blocX[[bloc]]))
  
  if (is.null(TOP)){
    if(type=="cov"){
      for(i in 1:ncol(blocX[[bloc]])){CovVar1[i] <- cov(blocX[[bloc]][,i],res$variates[[bloc]][,comp1])}
      for(i in 1:ncol(blocX[[bloc]])){CovVar2[i] <- cov(blocX[[bloc]][,i],res$variates[[bloc]][,comp2])}
      
      xx <- data.frame(CovVar1,CovVar2)
      KEEP <- apply(abs(xx),1,function(x) any(x >= seuil))
      rownames(xx) <- colnames(blocX[[bloc]])
      xx <- xx[KEEP,]
      pl<- ggplot(xx,aes(xx[,1],xx[,2],label=as.character(label[KEEP]))) + geom_point()  + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') +theme_bw() + geom_text(aes(x=xx[,1],y=xx[,2]),check_overlap = FALSE, nudge_x = 0.01,hjust=0) + labs(title=main)
      return(pl)
    }
    
    if(type=="cor"){
      for(i in 1:ncol(blocX[[bloc]])){CovVar1[i] <- cor(blocX[[bloc]][,i],res$variates[[bloc]][,comp1])}
      for(i in 1:ncol(blocX[[bloc]])){CovVar2[i] <- cor(blocX[[bloc]][,i],res$variates[[bloc]][,comp2])}
      
      xx <- data.frame(CovVar1,CovVar2)
      KEEP <- apply(abs(xx),1,function(x) any(x >= seuil))
      rownames(xx) <- colnames(blocX[[bloc]])
      xx <- xx[KEEP,]
      pl <- ggplot(xx,aes(xx[,1],xx[,2],label=label[KEEP])) + geom_point()  + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') +theme_bw() + geom_text(aes(x=xx[,1],y=xx[,2]),check_overlap = FALSE, nudge_x = 0.01,hjust=0) + labs(title=main)
      return(pl)
    }}
  
  if (!is.null(TOP)) {
    
    
    
    if(type=="cov"){
      for(i in 1:ncol(blocX[[bloc]])){CovVar1[i] <- cov(blocX[[bloc]][,i],res$variates[[bloc]][,comp1])}
      for(i in 1:ncol(blocX[[bloc]])){CovVar2[i] <- cov(blocX[[bloc]][,i],res$variates[[bloc]][,comp2])}
      
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
      for(i in 1:ncol(blocX[[bloc]])){CovVar1[i] <- cor(blocX[[bloc]][,i],res$variates[[bloc]][,comp1])}
      for(i in 1:ncol(blocX[[bloc]])){CovVar2[i] <- cor(blocX[[bloc]][,i],res$variates[[bloc]][,comp2])}
      
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
