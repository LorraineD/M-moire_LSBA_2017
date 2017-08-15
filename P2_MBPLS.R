# Analyse multibloc PLS sur les données NCI-60
# LDA - Mémoire septembre 2017

# Chargement données & packages -----------------------------------

setwd("C:/Users/Lorraine/Desktop/Code_propre")

require(pls)
source("Settings_Data.R")
source("Fonctions/Fct_MBPLS.txt")
source("Fonctions/Graphe_MBPLS.R")

options("encoding" = "UTF-8")

##  1. Déterminer le nombre de composantes ----------------------------------------------

# 1.1 RMSE et taux de bonnes classification ------------

RMSE_K <- Prop.BC <- Q2_K <- matrix(NA,ncol=4,nrow=10)
EY <- CP <- list()
for(i in 1:2) {EY[[i]] <- list()
for(a in 1:10) EY[[i]][[a]] <- matrix(NA,ncol=7,nrow=n)}

for(i in 1:n){
  cat("i = ",i, "\n")
  XX <- list(NCI[[1]][-i,],NCI[[2]][-i,],NCI[[3]][-i,],NCI[[4]][-i,])
  YY <- Class.mat[-i,]
  
  for(a in 1:10){
    mod1 <- MBPLS(XX, YY, A=a, deflY =TRUE, ScaleX = FALSE)
    EY[[1]][[a]][i,] <- X.NR[i,] %*% mod1$CoeffBeta
    
    mod2 <- MBPLS(XX, YY, A=a, deflY =TRUE, ScaleX = TRUE)
    EY[[2]][[a]][i,] <- X.R[i,] %*% mod2$CoeffBeta
    
  }}


Press_k <- RMSE_k <- BC <- Q2_k <- matrix(NA,nrow=10, ncol=2)
for(i in 1:2){
  for(a in 1:10){
    Q2_k[a,i] <- Q2(EY[[i]][[a]],Class.C)[[1]]
    RMSE_k[a,i] <- RMSE(EY[[i]][[a]],Class.C,k=a)
    Press_k[a,i] <- PRESS_GOMCIA(Class.C,EY[[i]][[a]])
    BC[a,i] <- mean(apply(EY[[i]][[a]],1,which.max)==ClassY.Num)*100
  }
}

RMSE_k <- data.frame(RMSE_k)
BC <- data.frame(BC)
Q2_k <- data.frame(Q2_k)

colnames(RMSE_k) <- colnames(BC) <- c("Données Non-réduites","Données Réduites")

# Enregistrement des graphiques et tables 
# RMSE 
png(file=paste0(pathI,"/NCI_MBPLS_RMSE.png"), width = 777, height = 463)
ggplot(RMSE_k) + geom_point(aes(seq(1,10),RMSE_k[,1],color="Non-réduites")) + geom_line(aes(seq(1,10),RMSE_k[,1],color="Non-réduites",lty="Non-réduites")) + ylim(0.4,1)  + geom_point(aes(seq(1,10),RMSE_k[,2],color="Réduites")) + geom_line(aes(seq(1,10),RMSE_k[,2],color="Réduites",lty="Réduites")) +  labs(title="RMSE selon le nombre de composantes",x="Nombre de composantes",y="RMSE")  + theme_bw() + guides(color=guide_legend("Modèle :"),lty=guide_legend("Modèle :"))+ scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()
# Taux de bonnes classifications
png(file=paste0(pathI,"/NCI_MBPLS_BC.png"), width = 777, height = 463)
ggplot(BC) + geom_point(aes(seq(1,10),BC[,1],color="Non-réduites")) + geom_line(aes(seq(1,10),BC[,1],color="Non-réduites",lty="Non-réduites")) + ylim(0,100)  + geom_point(aes(seq(1,10),BC[,2],color="Réduites")) + geom_line(aes(seq(1,10),BC[,2],color="Réduites",lty="Réduites"))  +  labs(title="Proportion de bonnes classifications selon le nombre de composantes",x="Nombre de composantes",y="% de bonnes classifications")  + theme_bw() + guides(color=guide_legend("Modèle :"),lty=guide_legend("Modèle :"))+ scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()
# Table RMSE
print(xtable(t(RMSE_k),caption="Valeur du RMSE obtenu selon le nombre de composantes et selon le modèle",digits=2), floating = getOption("xtable.floating", TRUE), table.placement = getOption("xtable.table.placement", "H"),sanitize.rownames.function = bold,sanitize.colnames.function =bold,file=paste0(pathR,"/NCI_MBPLS_RMSE.tex"))
# Table bonnes classifications
print(xtable(t(BC),caption="Proportion de bonnes classifications (en pourcentage) obtenues selon le nombre de composantes et selon le modèle",digits=2), floating = getOption("xtable.floating", TRUE),sanitize.rownames.function = bold,sanitize.colnames.function =bold, table.placement = getOption("xtable.table.placement", "H"),scalebox="0.9",file=paste0(pathR,"/NCI_MBPLS_BC.tex"))


# Nombre de composantes choisi : 
Ncomp = c(8,6)

# 1.2 Capacité de prédiction des modèles ---------------------

# Qualité de prédiction
Comp.mod <- matrix(NA,ncol=2,nrow=6) 
colnames(Comp.mod) <- colnames(BC) 
rownames(Comp.mod) <- c("PRESS-CV","PRESS","RMSE-CV","RMSE","Bonne classification - CV","Bonne classification")

Comp.mod[1,1] <- Press_k[8,1] ; Comp.mod[1,2] <- Press_k[6,2]
Comp.mod[3,1] <- unlist(RMSE_k[8,1]) ; Comp.mod[3,2] <- unlist(RMSE_k[8,2])
Comp.mod[5,1] <- unlist(BC[8,1]) ; Comp.mod[5,2] <- unlist(BC[8,2])

mod1 <- MBPLS(NCI, Class.mat, A=8, deflY =TRUE, ScaleX = FALSE)
mod2 <- MBPLS(NCI, Class.mat, A=6, deflY =TRUE, ScaleX = TRUE)

Comp.mod[2,] <- c(PRESS_GOMCIA(Class.C,mod1$prediction),PRESS_GOMCIA(Class.C,mod2$prediction))
Comp.mod[4,] <- c(RMSE(Y=Class.C,Ychap=mod1$prediction,8),RMSE(Y=Class.C,Ychap=mod2$prediction,6))
Comp.mod[6,] <- c(mean(apply(mod1$prediction,1,which.max)==ClassY.Num),mean(apply(mod2$prediction,1,which.max)==ClassY.Num))*100

print(xtable(Comp.mod,caption="Capacité de prédiction pour les deux modèles",digits=2,scientific=T,label="NCI:MBPLS_COMP"), floating = getOption("xtable.floating", TRUE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_MBPLS_ModComp.tex"))


# Inertie
Inert <- matrix(NA, ncol=2, nrow=5)
Inert[1,] <- c(sum(mod1$InertieG),sum(mod2$InertieG))*100
Inert[2:5,1] <- colSums(mod1$InertieX)*100
Inert[2:5,2] <- colSums(mod2$InertieX)*100
colnames(Inert) <- colnames(BC)
rownames(Inert) <- c("Inertie de Y","Inertie de X1","Inertie de X2","Inertie de X3","Inertie de X4")
print(xtable(Inert,caption="Inertie expliquée pour les différentes matrices pour les deux modèles ",digits=2,scientific=T,label="NCI:MBPLS_IN"), floating = getOption("xtable.floating", TRUE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_MBPLS_Iner.tex"))


# 2. Modèle MBPLS - 6 composantes - données réduites ----------------------------
a = 6
mod <-  MBPLS(NCI, Class.mat, A=a, deflY =TRUE, ScaleX = TRUE)


# 2.1 Classe prédite vs classe observée ----------------------------------------

YPred <- mod$prediction
for(i in 1:n){
  XX <- list(NCI[[1]][-i,],NCI[[2]][-i,],NCI[[3]][-i,],NCI[[4]][-i,])
  YY <- Class.mat[-i,]
  
  modCV <- MBPLS(XX, YY, A=6, deflY =TRUE, ScaleX = TRUE)
  YPred[i,] <- X.R[i,] %*% modCV$CoeffBeta
}

table(apply(YPred,1,which.max),ClassY.Num)
colnames(Class.mat)

# 2.2 Projection des observations ------------------------------------------

Plot.TK(mod, comp=c(1,2), bloc="G", ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(3,4), bloc="G", ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(5,6), bloc="G", ClassY=ClassY, ellipse=TRUE)

Plot.TK(mod, comp=c(1,2), bloc=1, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(3,4), bloc=1, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(5,6), bloc=1, ClassY=ClassY, ellipse=TRUE)

Plot.TK(mod, comp=c(1,2), bloc=2, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(3,4), bloc=2, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(5,6), bloc=2, ClassY=ClassY, ellipse=TRUE)

Plot.TK(mod, comp=c(1,2), bloc=3, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(3,4), bloc=3, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(5,6), bloc=3, ClassY=ClassY, ellipse=TRUE)

Plot.TK(mod, comp=c(1,2), bloc=4, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(3,4), bloc=4, ClassY=ClassY, ellipse=TRUE)
Plot.TK(mod, comp=c(5,6), bloc=4, ClassY=ClassY, ellipse=TRUE)


png(file=paste0(pathI,"/NCI_MBPLS_T12.png"), width = 777, height = 600)
Plot.TK(mod, comp=c(1,2), bloc="G", ClassY=ClassY, ellipse=TRUE)
dev.off()
png(file=paste0(pathI,"/NCI_MBPLS_T46.png"), width = 777, height = 600)
Plot.TK(mod, comp=c(4,6), bloc="G", ClassY=ClassY, ellipse=TRUE)
dev.off()
png(file=paste0(pathI,"/NCI_MBPLS_T45.png"), width = 777, height = 600)
Plot.TK(mod, comp=c(4,5), bloc="G", ClassY=ClassY, ellipse=TRUE)
dev.off()

# 2.3  Importance des variables --------------------------------------------
# 
# 2.3.1 VIP ----------------------------------------------------------------

R2 <- matrix(NA,ncol=7,nrow=6)

for(k in 1:6){
  res <- MBPLS(NCI,Class.mat,A = k, deflY =TRUE, ScaleX = TRUE)
  for(i in 1:7){
    PRESS <-  sum((Class.C[,i] - res$prediction[,i])^2)
    DEN <- sum((Class.C - colMeans(Class.C))^2)
    R2[k,i] <- 1 - PRESS / DEN
  }
}

B2 <- mod$globalweight
for(k in 1:6){
  B2[,k] <- mod$globalweight[,k] / sqrt( t(mod$globalweight[,k]) %*% mod$globalweight[,k])
}

VIP <- numeric(q)
for(j in 1:q){
  NUM <- sum(rowSums(R2) * B2[j,]^2)
  DEN <- sum(rowSums(R2))
  VIP[j] <- sqrt(q*NUM/DEN )
}

# Proportion et nombre de variables significatives
mean(VIP > mean(VIP))
sum(VIP > mean(VIP))

for(i in 1:4) {print(mean(VIP[N0[i]:NT[i]]>mean(VIP)))}


names(VIP) <- rownames(mod$CoeffBeta)
VIP_1 <- head(VIP[order(VIP,decreasing=T)], 10)
VIP_m <- matrix(VIP_1,ncol=1)
rownames(VIP_m) <- names(VIP_1)
colnames(VIP_m) <- "VIP"

# Boxplots des variables 
for(i in 1:10)  {
  png(file=paste0(pathI,"/NCI_MBPLS_BOX_V",i,".png"), width = 500, height = 400)
  boxplot(X.NR[,rownames(VIP_m)[i]]~ClassY, main=paste("Boxplot de la variable",rownames(VIP_m)[i]),ylab=paste0("Expression de ",rownames(VIP_m)[i]), cex.main=2, cex.lab=1.5)
  dev.off()
}

x <- data.frame(matrix(NA,ncol=4,nrow=5))
x[,1] <- names(VIP_m[1:5,])
x[,3] <- names(VIP_m[6:10,])
x[,2] <- format(VIP_m[1:5,],digits=3)
x[,4] <- format(VIP_m[6:10,],digits=3)

print(xtable(t(x),caption="Les 10 variables et leur coefficient ayant une valeur du VIP les plus importantes",digits=2,scientific=T), floating = getOption("xtable.floating", TRUE),include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames=getOption("xtable.include.colnames",FALSE),hline.after = c(1,3),file=paste0(pathR,"/NCI_MBPLS_VIP.tex"),table.placement = getOption("xtable.table.placement", "H"), sanitize.text.function = bold.somerows)


# 2.3.2 Coefficients de régression ------------------------------------------

xx <- mod$CoeffBeta
for(i in 1:7){assign(paste0("B",i), xx[order(abs(xx[,i]),decreasing=TRUE),i])}



BETA <- matrix(NA,ncol=14,nrow=5)
for(i in 1:7){
  BETA[,(2*i-1)] <- names(get(paste0("B",i)))[1:5]
  BETA[,(2*i)] <- format(get(paste0("B",i))[1:5],digits=2)
}

colnames(BETA) <- rep(colnames(Class.mat),each=2)


print(xtable(BETA[,1:8],caption="Les 5 variables avec les coefficients les plus importants en valeurs absolues et la valeur du coefficient de régression pour chaque variable réponse",label="tab:NCI:MBPLS-Beta1"),   include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after = 0,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_MBPLS-CoeffBeta1.tex"),floating=FALSE)

print(xtable(BETA[,9:14],caption="(suite) Les 5 variables avec les coefficients les plus importants en valeurs absolues et la valeur du coefficient de régression pour chaque variable réponse",label="tab:NCI:MBPLS-Beta2"), include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after =0,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_MBPLS-CoeffBeta2.tex"),floating=FALSE)

for(i in c(1,3,5,7,9,11,13))  {
  for(k in 1:5) {
    CC <- colnames(BETA)[i]
    png(file=paste0(pathI,"/NCI_MBPLS_BOX_",CC,k,".png"), width = 500, height = 400)
    boxplot(X.NR[,BETA[k,i]]~ClassY, main=paste("Boxplot de la variable",BETA[k,i]),ylab=paste0("Expression de ",BETA[k,i]), cex.main=2, cex.lab=1.5)
    dev.off()
  }
}


# 2.4 Liens entre variables et composantes ---------------------------------


CORR.VAR <- list()

for(i in 1:4){
  nv <- ncol(NCI[[i]])
  CORR.VAR[[i]] <- matrix(NA,nrow=nv,ncol=6)
  rownames(CORR.VAR[[i]]) <- colnames(NCI[[i]])
  for(j in 1:nv){
    for(k in 1:6){
      CORR.VAR[[i]][j,k] <- cor(NCI[[i]][,j],mod$blocscore[[i]][,k])
    }  }
}

nv <- ncol(X.NR)
CORR.VAR[[5]] <- matrix(NA,nrow=nv,ncol=6)
rownames(CORR.VAR[[5]]) <- colnames(X.NR)
for(j in 1:nv){
  for(k in 1:6){
    CORR.VAR[[5]][j,k] <- cor(X.NR[,j],mod$scoreX[,k])
  }  }


for(i in 1:4){
  MAT <- NULL
  for(k in 1:6){
    xx <- CORR.VAR[[i]]
    yy <- xx[order(abs(xx[,k]),decreasing=TRUE),k]
    names(yy) <- gsub("/.*", '',  names(yy))
    ss <- ifelse(yy>0,"+","-")
    MAT <- cbind(MAT,names(yy[1:10]),ss[1:10])}
  colnames(MAT) <- c("Comp. 1"," ","Comp. 2"," ","Comp. 3"," ","Comp. 4", " ","Composante 5"," ","Comp. 6"," ")
  print(xtable(MAT[1:5,],caption=paste("Les 5 variables explicatives les plus corrélées avec la composante - bloc $X_",i,"$")), include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after =0,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_MBPLS-VARCOV_",i,".tex"),floating=FALSE)
  
}

i = 5
MAT <- NULL
for(k in 1:6){
  xx <- CORR.VAR[[i]]
  yy <- xx[order(abs(xx[,k]),decreasing=TRUE),k]
  names(yy) <- gsub("/.*", '',  names(yy))
  ss <- ifelse(yy>0,"+","-")
  MAT <- cbind(MAT,names(yy[1:10]),ss[1:10])}
colnames(MAT) <- c("Comp. 1"," ","Comp. 2"," ","Comp. 3"," ","Comp. 4", " ","Composante 5"," ","Comp. 6"," ")
print(xtable(MAT[1:5,],caption=paste("Les 5 variables explicatives les plus corrélées avec le vecteur de scores global pour chaque composante")), include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after =0,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_MBPLS-VARCOV_",i,".tex"),floating=FALSE)
