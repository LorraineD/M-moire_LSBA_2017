# Chargement des données et des fonctions --------------------------------------

setwd("C:/Users/Lorraine/Desktop/Code_propre")
source("Settings_Data.R")

# Données
NCI1 <- list("agilent"=NCI[[1]],"hgu95"=NCI[[4]])
NCI1.R <- lapply(NCI1, scale, scale=TRUE)
NCI2 <- list("hgu133"=NCI[[2]],"hgu133p2"=NCI[[3]])
NCI2.R <- lapply(NCI2, scale, scale=TRUE)
nb.x <- length(NCI1)
nb.y <- length(NCI2)
NCI_Name <- c("agilent","hgu95","hgu133","hgu133p2")

# Fonctions 
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/acimogen.txt")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/acimogencv.txt")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/acimogen.graph")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/acimogen0.graph")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/acimogen0.txt")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/Dcentred.txt")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/Diagobloc.txt")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/RdtV.txt")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/GOMCIA/star.graph3")

options("encoding" = "UTF-8")


# 1. Nombre optimal de composantes -------------------------------------

res1 <- acimogen0(NCI1,NCI2,gen=1,corX=T,corY=T,A=10,impres=F,itermax=150)
res2 <- acimogen0(NCI1,NCI2,gen=2,corX=T,corY=T,A=10,impres=F,itermax=150)

# 1.1 RMSE -----------------------------------------------------

RMSE_K <- list()
for(i in 1:2) RMSE_K[[i]] <- matrix(NA,ncol=4,nrow=10,dimnames=list(paste0("Comp",seq(1:10)),NCI_Name))

for(a in 1:10){
  for(i in 1:nb.x) {
    RMSE_K[[1]][a,i] <- RMSE(res1$Xchap[[i]][[a]],NCI1.R[[i]],k=a)
    RMSE_K[[2]][a,i] <- RMSE(res2$Xchap[[i]][[a]],NCI1.R[[i]],k=a)
  }
  for(i in 1:nb.y) {
    RMSE_K[[1]][a,(i+2)] <- RMSE(res1$Ychap[[i]][[a]],NCI2.R[[i]],k=a)
    RMSE_K[[2]][a,(i+2)] <- RMSE(res2$Ychap[[i]][[a]],NCI2.R[[i]],k=a)
  }
  
}

RMSE_K[[1]] <- data.frame(RMSE_K[[1]])
RMSE_K[[2]] <- data.frame(RMSE_K[[2]])


png(file=paste0(pathI,"/NCI_GOMCIAC_RMSE1.png"), width = 777, height = 463)
ggplot() + geom_point(aes(seq(1,10),RMSE_K[[1]][,1],color="agilent")) + geom_line(aes(seq(1,10),RMSE_K[[1]][,1],color="agilent"))  + geom_point(aes(seq(1,10),RMSE_K[[1]][,2],color="hgu95")) + geom_line(aes(seq(1,10),RMSE_K[[1]][,2],color="hgu95"))  + geom_point(aes(seq(1,10),RMSE_K[[1]][,3],color="hgu133")) + geom_line(aes(seq(1,10),RMSE_K[[1]][,3],color="hgu133"))  + geom_point(aes(seq(1,10),RMSE_K[[1]][,4],color="hgu133p2")) + geom_line(aes(seq(1,10),RMSE_K[[1]][,4],color="hgu133p2")) +  labs(title="RMSE selon le nombre de composantes - GOMCIA1 données réduites",x="Nombre de composantes",y="RMSE")  + theme_bw() + guides(color=guide_legend("Bloc :"))+ scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()

png(file=paste0(pathI,"/NCI_GOMCIAC_RMSE2.png"), width = 777, height = 463)
ggplot() + geom_point(aes(seq(1,10),RMSE_K[[2]][,1],color="agilent")) + geom_line(aes(seq(1,10),RMSE_K[[2]][,1],color="agilent"))  + geom_point(aes(seq(1,10),RMSE_K[[2]][,2],color="hgu95")) + geom_line(aes(seq(1,10),RMSE_K[[2]][,2],color="hgu95"))  + geom_point(aes(seq(1,10),RMSE_K[[2]][,3],color="hgu133")) + geom_line(aes(seq(1,10),RMSE_K[[2]][,3],color="hgu133"))  + geom_point(aes(seq(1,10),RMSE_K[[2]][,4],color="hgu133p2")) + geom_line(aes(seq(1,10),RMSE_K[[2]][,4],color="hgu133p2")) +  labs(title="RMSE selon le nombre de composantes - GOMCIA2 données réduites",x="Nombre de composantes",y="RMSE")  + theme_bw() + guides(color=guide_legend("Bloc :")) + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()


png(file=paste0(pathI,"/NCI_GOMCIAC_RMSE5.png"), width = 777, height = 463)
ggplot() + geom_point(aes(seq(1,10),RMSE_K[[1]][,1],color="GOMCIA1")) + geom_line(aes(seq(1,10),RMSE_K[[1]][,1],color="GOMCIA1",lty="GOMCIA1"))   + geom_point(aes(seq(1,10),RMSE_K[[2]][,1],color="GOMCIA2")) + geom_line(aes(seq(1,10),RMSE_K[[2]][,1],color="GOMCIA2",lty="GOMCIA2")) +  labs(title="RMSE selon le nombre de composantes - bloc agilent",x="Nombre de composantes",y="RMSE")  + theme_bw()+ guides(color=guide_legend("Méthode :"),lty=guide_legend("Méthode :")) + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()

png(file=paste0(pathI,"/NCI_GOMCIAC_RMSE6.png"), width = 777, height = 463)
ggplot() + geom_point(aes(seq(1,10),RMSE_K[[1]][,2],color="GOMCIA1")) + geom_line(aes(seq(1,10),RMSE_K[[1]][,2],color="GOMCIA1",lty="GOMCIA1"))   + geom_point(aes(seq(1,10),RMSE_K[[2]][,2],color="GOMCIA2")) + geom_line(aes(seq(1,10),RMSE_K[[2]][,2],color="GOMCIA2",lty="GOMCIA2")) +  labs(title="RMSE selon le nombre de composantes - bloc hgu95",x="Nombre de composantes",y="RMSE")  + theme_bw() + guides(color=guide_legend("Méthode :"),lty=guide_legend("Méthode :")) + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()

png(file=paste0(pathI,"/NCI_GOMCIAC_RMSE7.png"), width = 777, height = 463)
ggplot() + geom_point(aes(seq(1,10),RMSE_K[[1]][,3],color="GOMCIA1")) + geom_line(aes(seq(1,10),RMSE_K[[1]][,3],color="GOMCIA1",lty="GOMCIA1"))   + geom_point(aes(seq(1,10),RMSE_K[[2]][,3],color="GOMCIA2")) + geom_line(aes(seq(1,10),RMSE_K[[2]][,3],color="GOMCIA2",lty="GOMCIA2")) +  labs(title="RMSE selon le nombre de composantes - bloc hgu133",x="Nombre de composantes",y="RMSE")  + theme_bw() + guides(color=guide_legend("Méthode :"),lty=guide_legend("Méthode :")) + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()

png(file=paste0(pathI,"/NCI_GOMCIAC_RMSE8.png"), width = 777, height = 463)
ggplot() + geom_point(aes(seq(1,10),RMSE_K[[1]][,4],color="GOMCIA1")) + geom_line(aes(seq(1,10),RMSE_K[[1]][,4],color="GOMCIA1",lty="GOMCIA1"))   + geom_point(aes(seq(1,10),RMSE_K[[2]][,4],color="GOMCIA2")) + geom_line(aes(seq(1,10),RMSE_K[[2]][,4],color="GOMCIA2",lty="GOMCIA2")) +  labs(title="RMSE selon le nombre de composantes - bloc hgu133p2",x="Nombre de composantes",y="RMSE")  + theme_bw() + guides(color=guide_legend("Méthode :"),lty=guide_legend("Méthode :")) + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()

print(xtable(t(RMSE_K[[1]]),caption="Valeur du RMSE pour les différents blocs de variables explicatives selon le nombre de composantes - Modèle GOMCIA1 avec données réduites",digits=3,scientific=F), floating = getOption("xtable.floating", FALSE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_GOMCIAC_RMSE1.tex"))



print(xtable(t(RMSE_K[[2]]),caption="Valeur du RMSE pour les différents blocs de variables explicatives selon le nombre de composantes - Modèle GOMCIA2 avec données réduites",digits=3,scientific=F), floating = getOption("xtable.floating", FALSE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_GOMCIAC_RMSE2.tex"))


## 1.2 Inertie expliquée ----------------------------------------
Inertie <- list()
for(i in 1:2){
  res <- get(paste0("res",i))
  Inertie[[i]] <- data.frame(res$inerXexppc[1,],res$inerXexppc[2,],res$inerYexppc[1,],res$inerYexppc[2,])
  colnames(Inertie[[i]]) <- NCI_Name
}


png(file=paste0(pathI,"/NCI_GOMCIAC_Inert1.png"), width = 777, height = 463)
ggplot() + geom_point(aes(seq(1,10),Inertie[[1]][,1],color="agilent")) + geom_line(aes(seq(1,10),Inertie[[1]][,1],color="agilent")) + geom_point(aes(seq(1,10),Inertie[[1]][,2],color="hgu95")) + geom_line(aes(seq(1,10),Inertie[[1]][,2],color="hgu95")) +  geom_point(aes(seq(1,10),Inertie[[1]][,3],color="hgu133")) + geom_line(aes(seq(1,10),Inertie[[1]][,3],color="hgu133")) +  geom_point(aes(seq(1,10),Inertie[[1]][,4],color="hgu133p2")) + geom_line(aes(seq(1,10),Inertie[[1]][,4],color="hgu133p2")) +  labs(title="% d'inertie expliquée selon le nombre de composantes - GOMCIA1",x="Nombre de composantes",y="Inertie expliquée (en %)")  + theme_bw() + guides(color=guide_legend("Bloc :"))  + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()

png(file=paste0(pathI,"/NCI_GOMCIAC_Inert2.png"), width = 777, height = 463)
ggplot() + geom_point(aes(seq(1,10),Inertie[[2]][,1],color="agilent")) + geom_line(aes(seq(1,10),Inertie[[2]][,1],color="agilent")) + geom_point(aes(seq(1,10),Inertie[[2]][,2],color="hgu95")) + geom_line(aes(seq(1,10),Inertie[[2]][,2],color="hgu95")) +  geom_point(aes(seq(1,10),Inertie[[2]][,3],color="hgu133")) + geom_line(aes(seq(1,10),Inertie[[2]][,3],color="hgu133")) +  geom_point(aes(seq(1,10),Inertie[[2]][,4],color="hgu133p2")) + geom_line(aes(seq(1,10),Inertie[[2]][,4],color="hgu133p2")) +  labs(title="% d'inertie expliquée selon le nombre de composantes - GOMCIA2",x="Nombre de composantes",y="Inertie expliquée (en %)")  + theme_bw() + guides(color=guide_legend("Bloc :")) + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()


print(xtable(t(Inertie[[1]]),caption="Inertie expliquée en pourcentage pour les différents blocs de variables explicatives selon le nombre de composantes - Modèle GOMCIA1 avec données réduites",digits=4,scientific=F), floating = getOption("xtable.floating", FALSE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_GOMCIAC_Inert1.tex"))

print(xtable(t(Inertie[[2]]),caption="Inertie expliquée en pourcentage  pour les différents blocs de variables explicatives selon le nombre de composantes - Modèle GOMCIA2 avec données réduites",digits=4,scientific=F), floating = getOption("xtable.floating", FALSE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_GOMCIAC_Inert2.tex"))


# 2. Modèle GOMCIA avec 3 composantes ------------------------------------------

A = 3
mod <-  acimogen0(NCI1,NCI2,gen=1,corX=T,corY=T,A=A,impres=F,itermax=150)



## 2.1 Projection des observations ---------------------------------------

png(file=paste0(pathI,"/NCI_GOMCIAC_TK12.png"), width = 777, height = 463)
Plot.TK(mod)
dev.off()
png(file=paste0(pathI,"/NCI_GOMCIAC_TK34.png"), width = 777, height = 463)
Plot.TK(mod,comp=c(1,3))
dev.off()


## 2.2 Proximité entre les blocs  ---------------------------------------
Tk <- list()
for(i in 1:A) Tk[[i]] <- cbind(mod$Tk[[i]],mod$Ul[[i]])
DD <- matrix(0,ncol=6,nrow=53) 
Pos <- combn(seq(1,4),2)
CC <- numeric(6)
for(i in 1:A){
  for(j in 1:6){
    x <- Pos[,j]
    DD[,j] <- DD[,j] + (Tk[[i]][,x[1]] - Tk[[i]][,x[2]])^2
    CC[j] <- CC[j] + cor(Tk[[i]][,x[1]],Tk[[i]][,x[2]])
    
  }}


table(apply(sqrt(DD),1,which.min))

table(apply(sqrt(DD),1,which.max))

colMeans(DD)
CC/4


## 2.3 Proximité entre les variables des différents blocs  ----------------------------                                              
# 2.3.1 Heatmap des poids des variables
require(gplots) ; require(dendextend)
Ak <- NULL
for(i in 1:2){
  Ak.a <- NULL
  for(j in 1:A) {
    ord1 <- order(abs(mod$Ak[[i]][,j]),decreasing = TRUE)
    data1 <- mod$Ak[[i]][ord1,]
    ord2 <- order(abs(mod$Bl[[i]][,j]),decreasing = TRUE)
    data2 <- mod$Bl[[i]][ord2,]
    Ak.a <- rbind(Ak.a,data1[1:10,],data2[1:10,])
  }
  Ak <- rbind(Ak,Ak.a)
}
Ak <- unique(Ak)
ID <- substr(rownames(Ak),1,1)
colnames(Ak) <- paste("Composante",1:A)

dAK <- dist(Ak)
hAk <- hclust(dAK)
labels(hAk) <- rep(" ",120)

cols_4 <- c("#CE4EB0","#2F5FF2","#4ECE59","#F29B2F")
col_car_type <- cols_4[factor(substr(rownames(Ak),1,1))]
dend <- as.dendrogram(hAk)
dend <- set(dend,"leaves_pch", 19)
dend <- set(dend,"leaves_col", col_car_type[order.dendrogram(dend)]) 
dend <- set(dend,"leaves_cex", 2)
png(file=paste0(pathI,"/GOMCIAC_ClusterAk.png"), width = 1500, height = 750)
plot(dend,main="Clustering hiérarchique sur base du poids des variables",horiz = FALSE,cex.lab=2, cex.axis=2, cex.main=4, cex.sub=2)
abline(h=0.25,lty=2)
legend("topright",legend=NCI_Name,col=cols_4,pch=19,cex=2.5)
dev.off()

groupe <- cutree(dend,k=4)
bloc <- as.factor(substr(rownames(Ak),1,1))
for(i in 1:4) print(table(as.factor(substr(rownames(Ak),1,1))[groupe==i])/table(groupe)[i])
table(as.factor(substr(rownames(Ak),1,1)),groupe)

distAk <- as.matrix(dAK)
distAk[upper.tri(distAk,diag=TRUE)] <- 1

for(i in 1:5) {
  coord <- Nmin(distAk,i)
  print(paste(colnames(distAk)[coord[2]],rownames(distAk)[coord[1]],distAk[coord]))
}

# 2.3.2 Clustering variables
require(made4)
COR <- Var.Covar(mod,comp=seq(1,3),type="cor")
rownames(COR) <- c(rownames(mod$Ak[[1]]),rownames(mod$Ak[[2]]),rownames(mod$Bl[[1]]),rownames(mod$Bl[[2]]))
ID <- as.factor(COR[,A+1])
COR <- COR[,1:A]
xx <- apply(COR,1,function(x) all(abs(x)>0.20))
distM <- as.matrix(as.dist(1-cor(t(COR[xx,]))))
distM[upper.tri(distM,diag=TRUE)] <- 2

colnames(COR) <- paste('Comp',1:A)
png(file=paste0(pathI,"/GOMCIAC_HeatCorr.png"), width = 1500, height = 750)
par(mar = c(5, 4, 1.4, 0.2))
heatplot(COR[xx,],dend='row',dualScale = T,classvec=ID[xx],classvecCol=cols_4,labRow=" ",labCol=paste('Composante',1:4),adjCol=c(0.5,1),srtCol=360)
add_legend("topright",legend=NCI_Name,col=cols_4,pch=19,cex=2.5,bty='n')
dev.off()
distM <- as.matrix(as.dist(1-cor(t(COR[xx,]))))
distM[upper.tri(distM,diag=TRUE)] <- 2

for(i in 1:5) {
  coord <- Nmin(distM,i)
  print(paste(colnames(distM)[coord[2]],rownames(distM)[coord[1]],distM[coord]))
}

hCor <- hclust(as.dist(1-cor(t(COR[xx,]))))
groupeC <- cutree(hCor,k=6)
for(i in 1:5) print(table(as.factor(substr(rownames(COR[xx,]),1,1))[groupeC==i])/table(groupeC)[i])


table(groupeC,as.factor(substr(rownames(COR[xx,]),1,1)))
