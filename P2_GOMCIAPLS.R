# Analyse GOMCIA-PLS sur les données NCI-60
# LDA - Mémoire septembre 2017

# Chargement données & packages -----------------------------------

setwd("C:/Users/Lorraine/Desktop/Code_propre")

require(pls)
source("Settings_Data.R")

source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/acimogen.txt")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/acimogencv.txt")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/acimogen.graph")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/acimogen0.graph")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/acimogen0.txt")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/Dcentred.txt")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/Diagobloc.txt")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/RdtV.txt")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/star.graph3")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/Graphe_GOMCIA.R")

options("encoding" = "UTF-8")

##  1. Déterminer le nombre de composantes ----------------------------------------------

# 1.1 RMSE et taux de bonnes classification ------------


blocY <-list(Class.mat)

Type <- c(1,1,3,3)

Pred_1 <- list() ; Pred_2 <- list() ; Pred_3 <- list(); Pred_4 <- list()
for(i in 1:10) { 
  Pred_1[[i]] <- matrix(NA,ncol=7,nrow=n)
  Pred_2[[i]] <- matrix(NA,ncol=7,nrow=n)
  Pred_3[[i]] <- matrix(NA,ncol=7,nrow=n)
  Pred_4[[i]] <- matrix(NA,ncol=7,nrow=n)}


for(i in 1:n){
  cat("i = ",i)
  XX <- list(NCI[[1]][-i,],NCI[[2]][-i,],NCI[[3]][-i,],NCI[[4]][-i,])
  YY <- list(blocY[[1]][-i,])
  
  res1 <- acimogen(XX,YY,1,centX=T,centY=T,corX=F,corY=F,A=10,impres=F)
  res2 <- acimogen(XX,YY,1,centX=T,centY=T,corX=T,corY=F,A=10,impres=F)
  res3 <- acimogen(XX,YY,3,centX=T,centY=T,corX=F,corY=F,A=10,impres=F)
  res4 <- acimogen(XX,YY,3,centX=T,centY=T,corX=T,corY=F,A=10,impres=F)
  
  for(a in 1:10){
    Pred_1[[a]][i,] <- t(X.NR[i,]) %*% res1$betaAchapconj[[a]]
    Pred_2[[a]][i,] <- t(X.R[i,])  %*% res2$betaAchapconj[[a]]
    Pred_3[[a]][i,] <- t(X.NR[i,]) %*% res3$betaAchapconj[[a]]
    Pred_4[[a]][i,] <- t(X.R[i,])  %*% res4$betaAchapconj[[a]]
  }}



R2_K <- Prop.Class <-  RMSE_K <- Q2_K <- PRESS_K <- matrix(NA,ncol=4,nrow=10)
for(j in 1:4){
  for(a in 1:10){
    x <- get(paste0("Pred_",j))[[a]]
    RMSE_K[a,j] <- RMSE(x,Class.C,k=a)
    PRESS_K[a,j] <- PRESS_GOMCIA(x,Class.C)
    Q2_K[a,j] <- Q2(x,Class.C)[[1]]
    Prop.Class[a,j] <- mean(apply(x,1,which.max) == ClassY.Num) *100
  }}


res1 <- acimogen(NCI,blocY,1,centX=T,centY=T,corX=F,corY=F,A=10,impres=F)
res2 <- acimogen(NCI,blocY,1,centX=T,centY=T,corX=T,corY=F,A=10,impres=F)
res3 <- acimogen(NCI,blocY,3,centX=T,centY=T,corX=F,corY=F,A=10,impres=F)
res4 <- acimogen(NCI,blocY,3,centX=T,centY=T,corX=T,corY=F,A=10,impres=F)

for(j in 1:4){
  for(a in 1:10){
    res <- get(paste0("res",j))
    R2_K[a,j] <- R2(res,A = a,gen=Type[j],bloc="Y") 
  }
}


PRESS_K <- data.frame(PRESS_K) ; R2_K <- data.frame(R2_K) ; Prop.Class <- data.frame(Prop.Class)
colnames(R2_K) <- colnames(PRESS_K) <- colnames(Prop.Class)<- c("1 - NR","1 - R","3 - NR","3 - R")

# Graphiques
png(file=paste0(pathI,"/NCI_GOMCIA_PRESS.png"), width = 777, height = 463)
ggplot(PRESS_K) + geom_point(aes(seq(1,10),PRESS_K[,1],color="1 - NR")) + geom_line(aes(seq(1,10),PRESS_K[,1],color="1 - NR")) + ylim(0,2)  + geom_point(aes(seq(1,10),PRESS_K[,2],color="1 - R")) + geom_line(aes(seq(1,10),PRESS_K[,2],color="1 - R"))  + geom_point(aes(seq(1,10),PRESS_K[,3],color="3 - NR")) + geom_line(aes(seq(1,10),PRESS_K[,3],color="3 - NR"))  + geom_point(aes(seq(1,10),PRESS_K[,4],color="3 - R")) + geom_line(aes(seq(1,10),PRESS_K[,4],color="3 - R")) +  labs(title="PRESS selon le nombre de composantes",x="Nombre de composantes",y="PRESS")  + theme_bw() + guides(color=guide_legend("Modèle :"))   + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()

png(file=paste0(pathI,"/NCI_GOMCIA_R2.png"), width = 777, height = 463)
ggplot(R2_K) + geom_point(aes(seq(1,10),R2_K[,1],color="1 - NR")) + geom_line(aes(seq(1,10),R2_K[,1],color="1 - NR")) + ylim(0,1)  + geom_point(aes(seq(1,10),R2_K[,2],color="1 - R")) + geom_line(aes(seq(1,10),R2_K[,2],color="1 - R"))  + geom_point(aes(seq(1,10),R2_K[,3],color="3 - NR")) + geom_line(aes(seq(1,10),R2_K[,3],color="3 - NR"))  + geom_point(aes(seq(1,10),R2_K[,4],color="3 - R")) + geom_line(aes(seq(1,10),R2_K[,4],color="3 - R")) +  labs(title="R2 selon le nombre de composantes",x="Nombre de composantes",y="R2")  + theme_bw() + guides(color=guide_legend("Modèle :"))  + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10))+ theme(text = element_text(size=18),legend.position="bottom")
dev.off()

png(file=paste0(pathI,"/NCI_GOMCIA_BC.png"), width = 777, height = 463)
ggplot(Prop.Class) + geom_point(aes(seq(1,10),Prop.Class[,1],color="1 - NR")) + geom_line(aes(seq(1,10),Prop.Class[,1],color="1 - NR")) + ylim(0,100)  + geom_point(aes(seq(1,10),Prop.Class[,2],color="1 - R")) + geom_line(aes(seq(1,10),Prop.Class[,2],color="1 - R"))  + geom_point(aes(seq(1,10),Prop.Class[,3],color="3 - NR")) + geom_line(aes(seq(1,10),Prop.Class[,3],color="3 - NR"))  + geom_point(aes(seq(1,10),Prop.Class[,4],color="3 - R")) + geom_line(aes(seq(1,10),Prop.Class[,4],color="3 - R")) +  labs(title="Taux de bonnes classifications selon le nombre de composantes",x="Nombre de composantes",y="% de bonnes classifications")  + theme_bw() + guides(color=guide_legend("Modèle :"))  + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()

print(xtable(t(PRESS_K),caption="Valeur du PRESS des différents modèles obtenus par cross-validation - NR : données non-réduites - R: données réduites",digits=2,scientific=F), floating = getOption("xtable.floating", TRUE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_GOMCIA_PRESS.tex"))

print(xtable(t(R2_K),caption="Valeur du $R^2$ des différents modèles obtenus par cross-validation - NR : données non-réduites - R: données réduites",digits=2,scientific=F), floating = getOption("xtable.floating", TRUE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_GOMCIA_R2.tex"))

print(xtable(t(Prop.Class),caption="Pourcentage d'observations bien classifiées pour les différents modèles obtenus par cross-validation - NR : données non-réduites - R: données réduites",digits=2,scientific=F), floating = getOption("xtable.floating", TRUE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_GOMCIA_BC.tex"))



# Nombres de composantes sélectionnées
A <- c(6,6,6,6)

# 1.2 Comparaison des modèles -------------------------------


Comparaison.Modele.1 <- matrix(NA,ncol=4,nrow=6)
rownames(Comparaison.Modele.1) <- c("PRESS-CV","PRESS","RMSE-CV","RMSE","Bonne classification - CV", "Bonne classificiation")
colnames(Comparaison.Modele.1) <- colnames(PRESS_K)

# Résultats par cross-validation
for(j in 1:4){
  Comparaison.Modele.1[1,j] <- PRESS_K[A[j],j]
  Comparaison.Modele.1[3,j] <- RMSE_K[A[j],j]
  Comparaison.Modele.1[5,j] <- Prop.Class[A[j],j]
}

# Les quatres modèles 
res1 <- acimogen(NCI,blocY,gen=1,centX=T,centY=T,corX=F,corY=F,A=6,impres=F)
res2 <- acimogen(NCI,blocY,gen=1,centX=T,centY=T,corX=T,corY=F,A=6,impres=F)
res3 <- acimogen(NCI,blocY,gen=3,centX=T,centY=T,corX=F,corY=F,A=6,impres=F)
res4 <- acimogen(NCI,blocY,gen=3,centX=T,centY=T,corX=T,corY=F,A=6,impres=F)

# Résultats sans cross-validation
for(j in 1:4){
  res <- get(paste0("res",j))
  Comparaison.Modele.1[2,j] <- PRESS_GOMCIA(Class.C,res$Ychap[[1]][[A[j]]])
  Comparaison.Modele.1[4,j] <- RMSE(Ychap = res$Ychap[[1]][[A[j]]],Y=Class.C,k=A[j])
  Comparaison.Modele.1[6,j] <- mean(apply(res$Ychap[[1]][[A[j]]],1,which.max) == ClassY.Num) *100
}


print(xtable(Comparaison.Modele.1,caption="Capacité de prédiction des différents modèles - NR : données non-réduites - R: données réduites",digits=2,scientific=T,label="NCI:GOMCIA_COMP"), floating = getOption("xtable.floating", TRUE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_GOMCIA_Comp1.tex"))


# Critère R2 pour les différents modèles

Comparaison.Modele.2 <- matrix(NA,nrow=5,ncol=4)
colnames(Comparaison.Modele.2) <- colnames(Comparaison.Modele.1)
rownames(Comparaison.Modele.2) <- c("R2 - Y", "R2 - X1", "R2 - X2","R2 - X3","R2 - X4")


for(j in 1:4){
  res <- get(paste0("res",j))
  Comparaison.Modele.2[1,j] <- R2(res,A=6,gen=Type[j],"Y") * 100
  Comparaison.Modele.2[2:5,j] <- R2(res,A=6,gen=Type[j],"X") *100}

print(xtable(Comparaison.Modele.2,caption="Coefficient de corrélation multiple des différents modèles - NR : données non-réduites - R: données réduites",digits=3,scientific=T,label="NCI:GOMCIA_COMP2"), floating = getOption("xtable.floating", TRUE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_GOMCIA_Comp2.tex"))

# 2. Modèle GOMCIA avec données réduites ---------------------------------------------

mod <- res4
A <- 6

# 2.1 Valeurs prédites vs valeurs observées -----------------------------------

table(apply(Pred_4[[6]],1,which.max),ClassY.Num)

# 2.2 Proximité bloc et composantes  --------------------------
png(file=paste0(pathI,"/GOMCIA_NCI-MU1.png"),width=446,height = 398)
Plot.MU(mod,xlab="Première composante",ylab="Deuxième composante", main="Coefficient µ - Composantes 1 & 2",label=c(1,2,3,4))
dev.off()

png(file=paste0(pathI,"/GOMCIA_NCI-MU2.png"),width=446,height = 398)
Plot.MU(mod,comp1=3,comp2=4,label=c(1,2,3,4),main="Coefficient µ - Composantes 3 & 4")
dev.off()

png(file=paste0(pathI,"/GOMCIA_NCI-MU3.png"),width=446,height = 398)
Plot.MU(mod,comp1=5,comp2=6,label=c(1,2,3,4),main="Coefficient µ - Composantes 5 & 6")
dev.off()

##  2.3 Projection des observations --------------------------



## Global - ALL 
Plot.TK(mod, ClassY = ClassY, comp=c(1,2), bloc = "G", ellipse=TRUE)
Plot.TK(mod, ClassY = ClassY, comp=c(3,4), bloc = "G", ellipse=TRUE)
Plot.TK(mod, ClassY = ClassY, comp=c(5,6), bloc = "G", ellipse=TRUE)

##  Bloc 1 - ALL 
Plot.TK(mod, ClassY = ClassY, comp=c(1,2), bloc = 1, ellipse=TRUE)
Plot.TK(mod, ClassY = ClassY, comp=c(3,4), bloc = 1, ellipse=TRUE)
Plot.TK(mod, ClassY = ClassY, comp=c(5,6), bloc = 1, ellipse=TRUE)

## Bloc 2 - ALL 
Plot.TK(mod, ClassY = ClassY, comp=c(1,2), bloc = 2, ellipse=TRUE)
Plot.TK(mod, ClassY = ClassY, comp=c(3,4), bloc = 2, ellipse=TRUE)
Plot.TK(mod, ClassY = ClassY, comp=c(5,6), bloc = 2, ellipse=TRUE)

## Bloc 3 - ALL
Plot.TK(mod, ClassY = ClassY, comp=c(1,2), bloc = 3, ellipse=TRUE)
Plot.TK(mod, ClassY = ClassY, comp=c(3,4), bloc = 3, ellipse=TRUE)
Plot.TK(mod, ClassY = ClassY, comp=c(5,6), bloc = 3, ellipse=TRUE)

## Bloc 4 - ALL
Plot.TK(mod, ClassY = ClassY, comp=c(1,2), bloc = 4, ellipse=TRUE)
Plot.TK(mod, ClassY = ClassY, comp=c(3,4), bloc = 4, ellipse=TRUE)
Plot.TK(mod, ClassY = ClassY, comp=c(5,6), bloc = 4, ellipse=TRUE)


# figure 
png(file=paste0(pathI,"/GOMCIAPLS_T12.png"),width=777,height = 600)
Plot.TK(mod, ClassY = ClassY, comp=c(1,2), bloc = "G", ellipse=TRUE)
dev.off()
png(file=paste0(pathI,"/GOMCIAPLS_T46.png"),width=777,height = 600)
Plot.TK(mod, ClassY = ClassY, comp=c(4,6), bloc = "G", ellipse=TRUE)
dev.off()
png(file=paste0(pathI,"/GOMCIAPLS_T45.png"),width=777,height = 600)
Plot.TK(mod, ClassY = ClassY, comp=c(4,5), bloc = "G", ellipse=TRUE)
dev.off()

## 2.4 Les variables ----------------------------------------

CORR.VAR <- list()

for(i in 1:4){
  nv <- ncol(NCI[[i]])
  CORR.VAR[[i]] <- matrix(NA,nrow=nv,ncol=6)
  rownames(CORR.VAR[[i]]) <- colnames(NCI[[i]])
  for(j in 1:nv){
    for(k in 1:6){
      CORR.VAR[[i]][j,k] <- cor(NCI[[i]][,j],mod$Tk[[k]][,i])
    }  }
}

nv <- ncol(X.NR)
CORR.VAR[[5]] <- matrix(NA,nrow=nv,ncol=6)
rownames(CORR.VAR[[5]]) <- colnames(X.NR)
SuperTK <- matrix(NA,nrow=n,ncol=6)
for(i in 1:6){
  SuperTK[,i] <- mod$Tk[[i]] %*% mod$MU[,i]
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
  print(xtable(MAT[1:5,],caption=paste("Les 5 variables explicatives les plus corrélées avec la composante - bloc $X_",i,"$")), include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after =0,table.placement = getOption("xtable.table.placement", "H"),floating=FALSE,file=paste0(pathR,"/NCI_GOMCIA-VARCOV_",i,".tex"))
  
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
print(xtable(MAT[1:5,],caption=paste("Les 5 variables explicatives les plus corrélées avec le vecteur de scores global pour chaque composante")), include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after =0,table.placement = getOption("xtable.table.placement", "H"),floating=FALSE,file=paste0(pathR,"/NCI_GOMCIA-VARCOV_",i,".tex"))



# 2.5 BIP -----------------------------------
D <- mod$D

BIP <- numeric(4)
TT <- NULL 
for(j in 1:6) {
  T_A <- mod$MU[1,j] * mod$Tk[[1]][,1] + mod$MU[2,j] * mod$Tk[[1]][,2] +mod$MU[3,j] * mod$Tk[[1]][,3] + mod$MU[4,j] * mod$Tk[[1]][,4]
  T_A <- matrix(T_A,ncol=1)
  B_A <- T_A %*% ginv( t(T_A) %*% D %*% T_A) %*% t(T_A) %*% D %*% blocY[[1]]
  assign(paste0("T",j),T_A)
  assign(paste0("A.",j),B_A)
  TT <- cbind(TT,T_A)
}
BTT  <- TT %*% ginv( t(TT) %*% D %*% TT) %*% t(TT) %*% D %*% blocY[[1]]
den <- 0
for(j in 1:6) {
  for(l in 1:4){
    den <- den + mod$MU[l,j] *normD(BTT,D)^2 
    assign(paste0("D",j),den)
  }}


for(k in 1:4){
  K <- 4  
  fr <- 0
  pBIP <- 0
  for(j in 1:6){
    fr <- mod$MU[k,j] * normD(get(paste0("A.",j)),D)^2 / (get(paste0("D",j)))
    pBIP <- pBIP + fr
    
  }
  BIP[k] <- sqrt(4 * pBIP)}




# 2.6 Importance des variables en prédiction -----------------------------------
# 2.6.1 Coefficient de régression ---------------------------------------------

xx <- mod$betaAchapconj[[6]]
rownames(xx) <- c(colnames(NCI[[1]]),colnames(NCI[[2]]),colnames(NCI[[3]]),colnames(NCI[[4]]))

for(i in 1:7){assign(paste0("B",i), xx[order(abs(xx[,i]),decreasing=TRUE),i])}


BETA <- matrix(NA,ncol=14,nrow=10)
for(i in 1:7){
  BETA[,(2*i-1)] <- names(get(paste0("B",i)))[1:10]
  BETA[,(2*i)] <- format(get(paste0("B",i))[1:10],digits=2)
}

colnames(BETA) <- rep(colnames(Class.mat),each=2)

print(xtable(BETA[1:5,1:8],caption="Les 5 variables avec les coefficients les plus importants en valeurs absolues et la valeur du coefficient de régression pour chaque variable réponse",label="tab:GOMCIA-Beta"),   include.rownames = getOption("xtable.include.rownames", FALSE),table.placement = getOption("xtable.table.placement", "H"),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after =0,file=paste0(pathR,"/NCI_GOMCIA-CoeffBeta1.tex"),floating=FALSE)

print(xtable(BETA[1:5,9:14],caption="(suite) Les 5 variables avec les coefficients les plus importants en valeurs absolues et la valeur du coefficient de régression pour chaque variable réponse",label="tab:GOMCIA-Beta2"), include.rownames = getOption("xtable.include.rownames", FALSE),table.placement = getOption("xtable.table.placement", "H"),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after = 0,file=paste0(pathR,"/NCI_GOMCIA-CoeffBeta2.tex"),floating=FALSE)


Name <- c("SNC","Colon","Poumon","Leucémique","Mélanome","Ovarien","Rénal")

for(i in c(1,3,5,7,9,11,13))  {
  for(k in 1:5) {
    CC <- colnames(BETA)[i]
    png(file=paste0(pathI,"/NCI_GOMCIA_BOX_",CC,k,".png"), width = 500, height = 400)
    boxplot(X.NR[,BETA[k,i]]~ClassY, main=paste("Boxplot de la variable",BETA[k,i]),ylab=paste0("Expression de ",BETA[k,i]), cex.main=2, cex.lab=1.5)
    dev.off()
  }
}

# 2.6.2 ---------------------------

R2 <- matrix(NA,ncol=7,nrow=6)
res <- acimogen(NCI,blocY,gen=3,centX=T,centY=T,corX=T,corY=F,A=6,impres=F)
for(k in 1:6){
  
  for(i in 1:7){
    PRESS <-  sum((Class.C[,i] - res$Ychap[[1]][[k]][,i])^2)
    DEN <- sum((Class.C - colMeans(Class.mat))^2)
    R2[k,i] <- 1 - PRESS / DEN
  }
}
globalweight <- rbind(res$Ak[[1]],res$Ak[[2]],res$Ak[[3]],res$Ak[[4]])
B2 <- globalweight
for(k in 1:6) B2[,k] <- globalweight[,k] / sqrt( t(globalweight[,k]) %*% globalweight[,k])


VIP <- numeric(q)
for(j in 1:q){
  NUM <- sum(rowSums(R2) * B2[j,]^2)
  DEN <- sum(rowSums(R2))
  VIP[j] <- sqrt(q*NUM/DEN )
}

mean(VIP > mean(VIP)) ;sum(VIP > mean(VIP)) 

for(i in 1:4) {print(mean(VIP[N0[i]:NT[i]]>mean(VIP)))}
names(VIP) <- c(rownames(res$Ak[[1]]),rownames(res$Ak[[2]]),rownames(res$Ak[[3]]),rownames(res$Ak[[4]]))
VIP_1 <- head(VIP[order(VIP,decreasing=T)], 10)
VIP_m <- matrix(VIP_1,ncol=1)
rownames(VIP_m) <- names(VIP_1)
colnames(VIP_m) <- "VIP"

x <- data.frame(matrix(NA,ncol=4,nrow=5))
x[,1] <- names(VIP_m[1:5,])
x[,3] <- names(VIP_m[6:10,])
x[,2] <- format(VIP_m[1:5,],digits=3)
x[,4] <- format(VIP_m[6:10,],digits=3)

print(xtable(t(x),caption="Les 10 variables et leur coefficient ayant une valeur du VIP les plus importantes",digits=2,scientific=T), floating = getOption("xtable.floating", TRUE),include.rownames = getOption("xtable.include.rownames", FALSE),sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"), hline.after = c(1,3),file=paste0(pathR,"/NCI_GOMCIA_VIP.tex"), sanitize.text.function = bold.somerows,include.colnames=getOption("xtable.include.colnames",FALSE))

for(i in 1:10)  {
  png(file=paste0(pathI,"/NCI_GOMCIA_BOX_V",i,".png"), width = 500, height = 400)
  boxplot(X.NR[,rownames(VIP_m)[i]]~ClassY, main=paste("Boxplot de la variable",rownames(VIP_m)[i]),ylab=paste0("Expression de ",rownames(VIP_m)[i]), cex.main=2, cex.lab=1.5)
  dev.off()
}

