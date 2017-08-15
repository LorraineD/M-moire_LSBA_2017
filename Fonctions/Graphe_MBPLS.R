Plot.TK <- function(res,comp=c(1,2),bloc=1,main=paste("Projection des observations sur les composantes",comp[1],"et",comp[2],"- bloc",bloc,sep=" "),xlab=paste("Composante",comp[1]),ylab=paste("Composante",comp[2]),ClassY=NULL,ellipse=FALSE){
  
  if (bloc=="G") {xx <- data.frame(res$scoreX[,comp])}    
  else {if(!is.numeric(bloc)) {stop("le paramètre BLOC est incorrect \n Valeurs possibles : G ou valeur numérique")}
    xx <- data.frame(res$blocscore[[bloc]][,comp]) }
  
  ClassY <- as.factor(ClassY)
  if(nlevels(ClassY) > 6){
    if (ellipse) {
      ggplot(xx,aes(xx[,1],xx[,2],col=ClassY,shape=ClassY))  + scale_shape_manual(values=1:nlevels(ClassY)) + geom_point(size=2)  + labs(title=main) + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :"),shape=guide_legend("Classe des observations :")) + theme(text = element_text(size=18),legend.position="bottom") + stat_ellipse(type="norm",level=0.9)
      }
    
    else { ggplot(xx,aes(xx[,1],xx[,2],col=ClassY,shape=ClassY))  + scale_shape_manual(values=1:nlevels(ClassY)) + geom_point(size=2) + labs(title=main) + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :"),shape=guide_legend("Classe des observations :")) +   theme(legend.position="bottom",plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=18) , legend.text=element_text(size=16))}}
  
  else {  
    if (ellipse) {
      ggplot(xx,aes(xx[,1],xx[,2],col=ClassY,pch=ClassY)) + geom_point(size=2) + labs(title=main) + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :"),pch=guide_legend(("Classe des observations :"))) +   theme(legend.position="bottom",plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=16) , legend.text=element_text(size=16)) + stat_ellipse(type="norm",level=0.9) }
    else {    ggplot(xx,aes(xx[,1],xx[,2],col=ClassY,pch=ClassY)) + geom_point(size=2) + labs(title=main) + ylab(ylab) + xlab(xlab) + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :"),pch=guide_legend(("Classe des observations :"))) +   theme(legend.position="bottom",plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=16) , legend.text=element_text(size=16))
    }}
}