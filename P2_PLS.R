# Analyse PLS sur les données NCI-60
# LDA - Mémoire septembre 2017

# Chargement données & packages -----------------------------------

setwd("C:/Users/Lorraine/Desktop/Code_propre")

require(pls)
source("Settings_Data.R")
source("Fonctions/Graphe_MBPLS.R")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/Graphe_PLS.R")

options("encoding" = "UTF-8")

##  1. Déterminer le nombre de composantes ----------------------------------------------

# 1.1 RMSE et taux de bonnes classification ------------

EY <- CP <- list()

for(i in 1:2) {
  EY[[i]]<- list()
  for(a in 1:10){
    EY[[i]][[a]] <- matrix(NA,ncol=7,nrow=n)
  }
}

for(i in 1:n){
  cat("i = ",i, "\n")
  
  mod1 <- mvr(Class.mat[-i,] ~ X.R[-i,],10,method="oscorespls",scale=FALSE)
  mod2 <- mvr(Class.mat[-i,] ~ X.NR[-i,],10,method="oscorespls",scale=FALSE)
  for(a in 1:10){
    EY[[1]][[a]][i,] <- predict(mod1,comps=1:a,newdata=X.R)[i,]
    EY[[2]][[a]][i,] <- predict(mod2,comps=1:a,newdata=X.NR)[i,]
  }}

Press_k <- RMSE_k <- BC <- Q2_k <- matrix(NA,ncol=2,nrow=10)
colnames(BC) <- colnames(RMSE_k) <- c("Données réduites","Données non-réduites")
for(i in 1:2){
  for(a in 1:10){
    Q2_k[a,i] <- Q2(EY[[i]][[a]],Class.C)[[1]]
    RMSE_k[a,i] <- RMSE(EY[[i]][[a]],Class.C,k=a)
    Press_k[a,i] <- PRESS_GOMCIA(Class.C,EY[[i]][[a]])
    BC[a,i] <- mean(apply(EY[[i]][[a]],1,which.max)==ClassY.Num)*100
  }}


RMSE_k <- data.frame(RMSE_k)
BC <- data.frame(BC)

png(file=paste0(pathI,"/NCI_PLS_RMSE.png"), width = 777, height = 463)
ggplot(RMSE_k) + geom_point(aes(seq(1,10),RMSE_k[,1],col="Réduites")) + geom_line(aes(seq(1,10),RMSE_k[,1],col="Réduites",lty="Réduites")) + geom_point(aes(seq(1,10),RMSE_k[,2],col="Non-réduites")) + geom_line(aes(seq(1,10),RMSE_k[,2],col="Non-réduites",lty="Non-réduites")) + labs(title="RMSE selon le nombre de composantes",x="Nombre de composantes",y="RMSE")  + theme_bw() + guides(color=guide_legend("Modèle :"),lty=guide_legend("Modèle :")) + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()

png(file=paste0(pathI,"/NCI_PLS_BC.png"), width = 777, height = 463)
ggplot(BC) + geom_point(aes(seq(1,10),BC[,1],col="Réduites")) + geom_line(aes(seq(1,10),BC[,1],col="Réduites",lty="Réduites")) + geom_point(aes(seq(1,10),BC[,2],col="Non-réduites")) + geom_line(aes(seq(1,10),BC[,2],col="Non-réduites",lty="Non-réduites")) + labs(title="% de bonnes classifications selon le nombre de composantes",x="Nombre de composantes",y="% de bonnes classifications")  + theme_bw() + guides(color=guide_legend("Modèle :"),lty=guide_legend("Modèle :")) + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()


print(xtable(t(RMSE_k),caption="Valeur du RMSE obtenu selon le nombre de composantes et selon le modèle",digits=2), floating = getOption("xtable.floating", TRUE), table.placement = getOption("xtable.table.placement", "H"),sanitize.rownames.function = bold,sanitize.colnames.function =bold,file=paste0(pathR,"/NCI_PLS_RMSE.tex"))

print(xtable(t(BC),caption="Proportion de bonnes classifications (en pourcentage) obtenues selon le nombre de composantes et selon le modèle",digits=2), floating = getOption("xtable.floating", TRUE),sanitize.rownames.function = bold,sanitize.colnames.function =bold, table.placement = getOption("xtable.table.placement", "H"),scalebox="0.9",file=paste0(pathR,"/NCI_PLS_BC.tex"))

# 1.2 Qualité de prédiction -----------------------------------------------------

mod <- mvr(Class.mat ~ X.R,6,method="oscorespls",scale=TRUE)

pred10 <- predict(mod,newdata=X.R,1:6,type="response")
pred <- pred10[,,6]

COMP <- matrix(NA,nrow=7,ncol=1)
rownames(COMP) <- c("PRESS-CV","PRESS","RMSE-CV","RMSE","Bonne classification - CV","Bonne classification","Inertie de Y expliquée")
colnames(COMP) <- "PLS"

Tk <- mod$scores
c <- mod$Yloadings
blocY <- Class.mat
I.Y <- numeric(6)
for(i in 1:6) I.Y[i] <- t(Tk[,i]) %*% Tk[,i] %*% t(c[,i]) %*% c[,i] / tr(t(blocY)%*%blocY)

COMP[1,] <- Press_k[6,1]
COMP[2,] <- PRESS_GOMCIA(Ychap=pred,Y=Class.C)
COMP[3,] <- RMSE_k[6,1]
COMP[4,] <- RMSE(pred,Class.C,6)
COMP[5,] <- BC[6,1]
COMP[6,] <-  mean(apply(pred,1,which.max)==ClassY.Num)*100
COMP[7,] <- sum(I.Y)

print(xtable(COMP,caption="Capacité de prédiction pour le modèle PLS avec données réduites et 6 composantes",digits=2,scientific=T,label="NCI:PLS_COMP"), floating = getOption("xtable.floating", TRUE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_PLS_Comp.tex"))

# 2. Modèle PLS avec données réduites ----------------------------------------

# 2.1 Projections des observations --------------------------------------------

xx <- mod$scores
class(xx) <- NULL
xx <- data.frame(xx)

png(file=paste0(pathI,"/NCI_PLS_T12.png"), width = 777, height = 600)
ggplot(xx,aes(xx[,1],xx[,2],col=ClassY)) + geom_point(size=2) + labs(title="Projection des observations sur les deux premières composantes") + ylab("Composante 2") + xlab("Composante 1") + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :"))+ theme(text = element_text(size=18),legend.position="bottom") + stat_ellipse(type = "norm", level = 0.9)
dev.off()

png(file=paste0(pathI,"/NCI_PLS_T34.png"), width = 777, height = 600)
ggplot(xx,aes(xx[,3],xx[,4],col=ClassY)) + geom_point(size=2) + labs(title="Projection des observations sur composantes 3 et 4") + ylab("Composante 4") + xlab("Composante 3") + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :")) + theme(text = element_text(size=18),legend.position="bottom") + stat_ellipse(type = "norm", level = 0.9)
dev.off()

png(file=paste0(pathI,"/NCI_PLS_T56.png"), width = 777, height = 600)
ggplot(xx,aes(xx[,5],xx[,6],col=ClassY)) + geom_point(size=2) + labs(title="Projection des observations sur les composantes 5 et 6") + ylab("Composante 6") + xlab("Composante 5") + geom_vline(xintercept=0,colour="grey") + geom_hline(yintercept=0,colour='grey') + theme_bw() + guides(color=guide_legend("Classe des observations :"))+ theme(text = element_text(size=18),legend.position="bottom") + stat_ellipse(type = "norm", level = 0.9)
dev.off()


# 2.2 importance des variables en prédiction ------------------------

# 2.2.1  Coefficients de régression --------------------------------------

xx <- mod$coefficients[,,6]
for(i in 1:7){assign(paste0("B",i), xx[order(abs(xx[,i]),decreasing=TRUE),i])}

BETA <- matrix(NA,ncol=14,nrow=5)
for(i in 1:7){
  BETA[,(2*i-1)] <- names(get(paste0("B",i)))[1:5]
  BETA[,(2*i)] <- format(get(paste0("B",i))[1:5],digits=2)
}

colnames(BETA) <- rep(colnames(Class.mat),each=2)

print(xtable(BETA[,1:8],caption="Les 5 variables avec les coefficients les plus importants en valeurs absolues et la valeur du coefficient de régression pour chaque variable réponse",label="tab:PLS-Beta1"),   include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after = 0,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_PLS-CoeffBeta1.tex"),floating=FALSE)

print(xtable(BETA[,9:14],caption="(suite) Les 5 variables avec les coefficients les plus importants en valeurs absolues et la valeur du coefficient de régression pour chaque variable réponse",label="tab:PLS-Beta2"), include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames = getOption("xtable.include.colnames", TRUE),sanitize.colnames.function =bold, hline.after =0,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_PLS-CoeffBeta2.tex"),floating=FALSE)

Name <- c("SNC","Colon","Poumon","Leucémique","Mélanome","Ovarien","Rénal")


for(i in c(1,3,5,7,9,11,13))  {
  for(k in 1:5) {
    CC <- colnames(BETA)[i]
    png(file=paste0(pathI,"/NCI_PLS_BOX_",CC,k,".png"), width = 500, height = 400)
    boxplot(X.NR[,BETA[k,i]]~ClassY, main=paste("Boxplot de la variable",BETA[k,i]),ylab=paste0("Expression de ",BETA[k,i]), cex.main=2, cex.lab=1.5)
    dev.off()
  }
}


# 2.2.2 Coefficient VIP -----------------------------------------------------------

R2 <- matrix(NA,ncol=7,nrow=6)

for(k in 1:6){
  for(i in 1:7){
    PRESS <-  sum((Class.mat[,i] - pred10[,i,k])^2)
    DEN <- sum((Class.mat - colMeans(Class.mat))^2)
    R2[k,i] <- 1 - PRESS / DEN
  }
}

B2 <- mod$loadings[,1:6]
for(k in 1:6){
  B2[,k] <- mod$loadings[,k] / sqrt( t(mod$loadings[,k]) %*% mod$loadings[,k])
}

VIP <- numeric(q)
for(j in 1:q){
  NUM <- sum(rowSums(R2) * B2[j,]^2)
  DEN <- sum(rowSums(R2))
  VIP[j] <- sqrt(q*NUM/DEN )
}

mean(VIP > mean(VIP))

names(VIP) <- rownames(mod$loadings)
VIP_1 <- head(VIP[order(VIP,decreasing=T)], 10)
VIP_m <- matrix(VIP_1,ncol=1)
rownames(VIP_m) <- names(VIP_1)
colnames(VIP_m) <- "VIP"

x <- data.frame(matrix(NA,ncol=4,nrow=5))
x[,1] <- names(VIP_m[1:5,])
x[,3] <- names(VIP_m[6:10,])
x[,2] <- format(VIP_m[1:5,],digits=3)
x[,4] <- format(VIP_m[6:10,],digits=3)


print(xtable(t(x),caption="Les 10 variables et leur coefficient ayant une valeur du VIP les plus importantes",digits=2,scientific=T), floating = getOption("xtable.floating", TRUE),include.rownames = getOption("xtable.include.rownames", FALSE),include.colnames=getOption("xtable.include.colnames",FALSE),hline.after = c(1,3),file=paste0(pathR,"/NCI_PLS_VIP.tex"), sanitize.text.function = bold.somerows)

for(k in 1:10) {
  png(file=paste0(pathI,"/NCI_PLS_BOX_V",k,".png"), width = 500, height = 400)
  boxplot(X.NR[,rownames(VIP_m)[k]]~ClassY, main=paste("Boxplot de la variable",rownames(VIP_m)[k]),ylab=paste0("Expression de ",rownames(VIP_m)[k]), cex.main=2, cex.lab=1.5)
  dev.off()
}


