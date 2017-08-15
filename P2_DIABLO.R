# Analyse RGCCA - prédiction (DIABLO) sur les données NCI-60
# LDA - Mémoire septembre 2017

# Chargement données & packages -----------------------------------

setwd("C:/Users/Lorraine/Desktop/Code_propre")

require(mixOmics)
source("Settings_Data.R")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/RGCCA/Graphe_RGCCA.R")

options("encoding" = "UTF-8")


C <- matrix( c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,1,1,1,0),ncol=5)
Tau <- rep(0,5)

##  1. Déterminer le nombre de composantes ----------------------------------------------

# 1.1 RMSE et taux de bonnes classification ------------



Pred_CV <- EY <- list()
for(i in 1:2) {
  EY[[i]] <- list()
  Pred_CV[[i]] <- matrix(NA,ncol=10,nrow=n)
  for(a in 1:10) EY[[i]][[a]] <- matrix(NA,ncol=7,nrow=n)}

for(i in 1:n){
  cat("i = ",i, "\n")
  XX <- list(agilent = NCI[[1]][-i,],hgu133=NCI[[2]][-i,],hgu133p2=NCI[[3]][-i,],hgu95=NCI[[4]][-i,])
  YY <- ClassY[-i]
  test <- list(agilent = as.matrix(t(NCI[[1]][i,])),hgu133= as.matrix(t(NCI[[2]][i,])),hgu133p2=as.matrix(t(NCI[[3]][i,])),hgu95=as.matrix(t(NCI[[4]][i,])))
  
  
  mod_NR <- block.plsda(XX,YY,ncomp=10, scale=FALSE, scheme="horst")
  mod_R <- block.plsda(XX,YY,ncomp=10, scale=TRUE, scheme="horst")
  for(a in 1:10){
    
    
    cat("   a = ",a, "\n")
    PP <- predict(mod_NR,newdata=test)
    EY[[1]][[a]][i,] <- PP$AveragedPredict[1,,a]
    Pred_CV[[1]][i,a] <- PP$AveragedPredict.class[[1]][,a]
    
    PP <- predict(mod_R,newdata=test)
    EY[[2]][[a]][i,] <- PP$AveragedPredict[1,,a]
    Pred_CV[[2]][i,a] <- PP$AveragedPredict.class[[1]][,a]
  }}


Press_k <- RMSE_k <- BC <- Q2_k <- matrix(NA,nrow=10, ncol=2)
for(i in 1:2){
  for(a in 1:10){
    Q2_k[a,i] <- Q2(EY[[i]][[a]],Class.C)[[1]]
    RMSE_k[a,i] <- RMSE(EY[[i]][[a]],Class.C,k=a)
    Press_k[a,i] <- PRESS_GOMCIA(Class.C,EY[[i]][[a]])
    BC[a,i] <- mean(Pred_CV[[i]][,a]==ClassY)*100
  }
}

RMSE_k <- data.frame(RMSE_k)
BC <- data.frame(BC)
Q2_k <- data.frame(Q2_k)

colnames(RMSE_k) <- colnames(BC) <- c("Données Non-réduites","Données Réduites")


png(file=paste0(pathI,"/NCI_RGCCA_RMSE.png"), width = 777, height = 463)
ggplot(RMSE_k) + geom_point(aes(seq(1,10),RMSE_k[,1],col="Non-réduites")) + geom_line(aes(seq(1,10),RMSE_k[,1],col="Non-réduites",lty="Non-réduites"))  + geom_point(aes(seq(1,10),RMSE_k[,2],color="Réduites")) + geom_line(aes(seq(1,10),RMSE_k[,2],color="Réduites",lty="Réduites")) +  labs(title="RMSE selon le nombre de composantes",x="Nombre de composantes",y="RMSE")  + theme_bw() + guides(color=guide_legend("Modèle :"),lty=guide_legend("Modèle :")) + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10))  + theme(text = element_text(size=18),legend.position="bottom")
dev.off()

png(file=paste0(pathI,"/NCI_RGCCA_BC.png"), width = 777, height = 463)
ggplot(BC) + geom_point(aes(seq(1,10),BC[,1],color="Non-réduites")) + geom_line(aes(seq(1,10),BC[,1],color="Non-réduites",lty="Non-réduites")) + ylim(15,100)  + geom_point(aes(seq(1,10),BC[,2],color="Réduites")) + geom_line(aes(seq(1,10),BC[,2],color="Réduites",lty="Réduites"))  +  labs(title="Pourcentage de bonnes classifications selon le nombre de composantes",x="Nombre de composantes",y="% de bonnes classifications")  + theme_bw() + guides(color=guide_legend("Modèle :"),lty=guide_legend("Modèle :")) + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10))  + theme(text = element_text(size=18),legend.position="bottom")
dev.off()

print(xtable(t(RMSE_k),caption="Valeur du RMSE obtenu selon le nombre de composantes et selon le modèle",digits=2), floating = getOption("xtable.floating", TRUE), table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_RGCCA_RMSE.tex"))

print(xtable(t(BC),caption="Pourcentage de bonnes classifications selon le nombre de composantes et selon le modèle",digits=2), floating = getOption("xtable.floating", TRUE), table.placement = getOption("xtable.table.placement", "H"),scalebox='0.95',file=paste0(pathR,"/NCI_RGCCA_BC.tex"))


Ncomp = c(6,6)


## 1.2 Capacité de prédiction des modèles --------------------------------------

Comp.mod <- matrix(NA,ncol=2,nrow=6) 
colnames(Comp.mod) <- colnames(BC) 
rownames(Comp.mod) <- c("PRESS-CV","PRESS","RMSE-CV","RMSE","Bonne classification - CV","Bonne classification")

Comp.mod[1,1] <- Press_k[Ncomp[1],1]
Comp.mod[3,1] <- RMSE_k[Ncomp[1],1]
Comp.mod[5,1] <- BC[Ncomp[1],1]

Comp.mod[1,2] <- Press_k[Ncomp[2],2]
Comp.mod[3,2] <- RMSE_k[Ncomp[2],2]
Comp.mod[5,2] <- BC[Ncomp[2],2]

mod_NR <- block.plsda(NCI,ClassY,ncomp=6, scale=FALSE, scheme="horst")
mod_R <- block.plsda(NCI,ClassY,ncomp=6, scale=TRUE, scheme="horst")

Pred_NR <- predict(mod_NR,newdata=NCI)
Pred_R <- predict(mod_R, newdata=NCI)

Comp.mod[2,] <- c(PRESS_GOMCIA(Y = Class.C,Pred_NR$AveragedPredict[,,6]),PRESS_GOMCIA(Y = Class.C,Pred_R$AveragedPredict[,,6]))
Comp.mod[4,] <- c(RMSE(Y = Class.C,Pred_NR$AveragedPredict[,,6]),RMSE(Y = Class.C,Pred_R$WeightedPredict[,,6]))
Comp.mod[6,] <- c(mean(Pred_NR$AveragedPredict.class[[1]][,6]==ClassY),mean(Pred_R$AveragedPredict.class[[1]][,6]==ClassY))*100

print(xtable(Comp.mod,caption="Capacité de prédiction pour les deux modèles avec 6 composantes",digits=3,scientific=T,label="NCI:RGCCA_COMP"), floating = getOption("xtable.floating", TRUE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_RGCCA_ModComp.tex"))


# 2. Modèle DIABLO - 6 composantes et données réduites --------------------------------------------------

A = 6
mod <- block.plsda(NCI,ClassY,ncomp=A, scale=TRUE, scheme="horst")

# 2.1 Classe prédite vs classe observée -----------------------------------------------

Pred <- predict(mod, newdata=NCI)
table(Pred$WeightedPredict.class[[1]][,A],ClassY)

table(Pred_CV[[2]][,A],ClassY)

# 2.2 Projections des observations ---------------------------------------

Plot.TK(mod, comp=c(1,2),bloc="G", ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(3,4),bloc="G", ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(5,6),bloc="G", ClassY=ClassY, ellipse=TRUE)

Plot.TK(mod, comp=c(1,2),bloc=1, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(3,4),bloc=1, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(5,6),bloc=1, ClassY=ClassY, ellipse=TRUE)

Plot.TK(mod, comp=c(1,2),bloc=2, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(3,4),bloc=2, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(5,6),bloc=2, ClassY=ClassY, ellipse=TRUE)

Plot.TK(mod, comp=c(1,2),bloc=3, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(3,4),bloc=3, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(5,6),bloc=3, ClassY=ClassY, ellipse=TRUE)

Plot.TK(mod, comp=c(1,2),bloc=4, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(3,4),bloc=4, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(5,6),bloc=4, ClassY=ClassY, ellipse=TRUE)



png(file=paste0(pathI,"/RGCCA_T12.png"), width = 777, height = 600)
Plot.TK(mod, comp=c(1,2),bloc=2, ClassY=ClassY, ellipse=TRUE)
dev.off()
png(file=paste0(pathI,"/RGCCA_T45.png"), width = 777, height = 600)
Plot.TK(mod, comp=c(4,5),bloc="G", ClassY=ClassY, ellipse=TRUE)
dev.off()
png(file=paste0(pathI,"/RGCCA_T14.png"), width = 777, height = 600)
Plot.TK(mod, comp=c(1,4),bloc="G", ClassY=ClassY, ellipse=TRUE)
dev.off()


# 2.3 Lien entre variables et composantes -------------------------------------

CORR.VAR <- list()

for(i in 1:4){
  nv <- ncol(NCI[[i]])
  CORR.VAR[[i]] <- matrix(NA,nrow=nv,ncol=6)
  rownames(CORR.VAR[[i]]) <- colnames(NCI[[i]])
  for(j in 1:nv){
    for(k in 1:6){
      CORR.VAR[[i]][j,k] <- cor(NCI[[i]][,j],mod$variates[[i]][,k])
    }  }
}

nv <- ncol(X.NR)
CORR.VAR[[5]] <- matrix(NA,nrow=nv,ncol=6)
rownames(CORR.VAR[[5]]) <- colnames(X.NR)
SuperTK <- matrix(0,nrow=n,ncol=6)
for(i in 1:6){
  for(b in 1:4){
    SuperTK[,i] <- SuperTK[,i] + mod$variates[[b]][,i]
    
  }
}
for(j in 1:nv){
  for(k in 1:6){
    CORR.VAR[[5]][j,k] <- cor(X.NR[,j],SuperTK[,k])
  }  }


for(i in 1:4){
  MAT <- NULL
  for(k in 1:6){
    xx <- CORR.VAR[[i]]
    yy <- xx[order(abs(xx[,k]),decreasing=TRUE),k]
    ss <- ifelse(yy>0,"+","-")
    MAT <- cbind(MAT,names(yy[1:10]),ss[1:10])}
  colnames(MAT) <- c("Comp. 1"," ","Comp. 2"," ","Comp. 3"," ","Comp. 4", " ","Composante 5"," ","Comp. 6"," ")
  print(xtable(MAT[1:5,],caption=paste("Les 5 variables explicatives les plus corrélées avec la composante - bloc $X_",i,"$")), include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after =0,table.placement = getOption("xtable.table.placement", "H"),floating=FALSE,file=paste0(pathR,"/NCI_DIABLO_VARCOV_",i,".tex"))
  
}

i = 5
MAT <- NULL
for(k in 1:6){
  xx <- CORR.VAR[[i]]
  yy <- xx[order(abs(xx[,k]),decreasing=TRUE),k]
  ss <- ifelse(yy>0,"+","-")
  MAT <- cbind(MAT,names(yy[1:10]),ss[1:10])}
colnames(MAT) <- c("Comp. 1"," ","Comp. 2"," ","Comp. 3"," ","Comp. 4", " ","Composante 5"," ","Comp. 6"," ")
print(xtable(MAT[1:5,],caption=paste("Les 5 variables explicatives les plus corrélées avec le vecteur de scores global pour chaque composante")), include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after =0,table.placement = getOption("xtable.table.placement", "H"),floating=FALSE,file=paste0(pathR,"/NCI_DIABLO_VARCOV_",i,".tex"))

# 2.4 Importance des variables en prédiction ------------------------------------
# 2.4.1 VIP --------------------------------------------------------------------- 

R2 <- matrix(NA,ncol=7,nrow=6)

for(k in 1:6){
  for(i in 1:7){
    PRESS <-  sum((Class.mat[,i] - Pred$WeightedPredict[,i,k])^2)
    DEN <- sum((Class.mat - colMeans(Class.mat))^2)
    R2[k,i] <- 1 - PRESS / DEN
  }
}


B2 <- matrix(NA, nrow=q, ncol=6)
for(k in 1:6){
  for(b in 1:nb){
    B2[N0[b]:NT[b],k] <- mod$loadings[[b]][,k]
  }
  B2[,k] <- B2[,k] / sqrt( t(B2[,k]) %*% B2[,k])
}

rownames(B2) <- colnames(X.R)

VIP <- numeric(q)
for(j in 1:q){
  NUM <- sum(rowSums(R2) * B2[j,]^2)
  DEN <- sum(rowSums(R2))
  VIP[j] <- sqrt(q*NUM/DEN )
}

sum(VIP > mean(VIP))

for(i in 1:4) {print(mean(VIP[N0[i]:NT[i]]>mean(VIP)))}
names(VIP) <- rownames(B2)
VIP_1 <- head(VIP[order(VIP,decreasing=T)], 10)
VIP_m <- matrix(VIP_1,ncol=1)
rownames(VIP_m) <- names(VIP_1)
colnames(VIP_m) <- "VIP"

x <- data.frame(matrix(NA,ncol=4,nrow=5))
x[,1] <- names(VIP_m[1:5,])
x[,3] <- names(VIP_m[6:10,])
x[,2] <- format(VIP_m[1:5,],digits=3)
x[,4] <- format(VIP_m[6:10,],digits=3)

print(xtable(t(x),caption="Les 10 variables et leur coefficient ayant une valeur du VIP les plus importantes",digits=2,scientific=T), floating = getOption("xtable.floating", TRUE),include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames=getOption("xtable.include.colnames",FALSE),sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"), hline.after = c(1,3),file=paste0(pathR,"/NCI_RGCCA_VIP.tex"), sanitize.text.function = bold.somerows)

for(i in 1:10)  {
  png(file=paste0(pathI,"/NCI_DIABLO_BOX_V",i,".png"), width = 500, height = 400)
  boxplot(X.NR[,rownames(VIP_m)[i]]~ClassY, main=paste("Boxplot de la variable",rownames(VIP_m)[i]),ylab=paste0("Expression de ",rownames(VIP_m)[i]), cex.main=2, cex.lab=1.5)
  dev.off()
}


# 2.4.2 Coefficient de régression -------------------------------

xx <- rbind(Pred$B.hat[[1]][,,6],Pred$B.hat[[2]][,,6],Pred$B.hat[[3]][,,6],Pred$B.hat[[4]][,,6])
for(i in 1:7){assign(paste0("B",i), xx[order(abs(xx[,i]),decreasing=TRUE),i])}



BETA <- matrix(NA,ncol=14,nrow=5)
for(i in 1:7){
  BETA[,(2*i-1)] <- names(get(paste0("B",i)))[1:5]
  BETA[,(2*i)] <- format(get(paste0("B",i))[1:5],digits=2)
}

colnames(BETA) <- rep(c("CN","CO","LC","LE","ME","OV","RE"),each=2)

print(xtable(BETA[,1:8],caption="Les 5 variables avec les coefficients les plus importants en valeurs absolues et la valeur du coefficient de régression pour chaque variable réponse",label="tab:RGCCA-Beta1"),   include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after = 0,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_RGCCA_Beta1.tex"))

print(xtable(BETA[,9:14],caption="(suite) Les 5 variables avec les coefficients les plus importants en valeurs absolues et la valeur du coefficient de régression pour chaque variable réponse",label="tab:RGCCA-Beta2"), include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after = 0,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_RGCCA_Beta2.tex"))



for(i in c(1,3,5,7,9,11,13))  {
  for(k in 1:5) {
    CC <- colnames(BETA)[i]
    png(file=paste0(pathI,"/NCI_DIABLO_BOX_",CC,k,".png"), width = 500, height = 400)
    boxplot(X.NR[,BETA[k,i]]~ClassY, main=paste("Boxplot de la variable",BETA[k,i]),ylab=paste0("Expression de ",BETA[k,i]), cex.main=2, cex.lab=1.5)
    dev.off()
  }
}

