# Chargement des données et des fonctions --------------------------------------

setwd("C:/Users/Lorraine/Desktop/Code_propre")
source("Settings_Data.R")

NCI_b <-  lapply(NCI,t)

require(mogsa)

options("encoding" = "UTF-8")


# 1. Nombre de composantes ----------------------------------------


modR <-mbpca(NCI_b, ncomp=10,method="globalScore",scale=TRUE,moa=FALSE,verbose=FALSE)

## 1.1 RMSE ----------------------------------------

RMSE_K <- X1  <- list()
RMSE_K[[1]] <- matrix(NA,ncol=4,nrow=10,dimnames=list(paste0("Comp",seq(1:10)),NCI_Name))

for(a in 1:10){
  for(i in 1:nb) {
    T.A <- modR$t[,1:a]
    X1[[i]] <- T.A %*% ginv(t(T.A)%*%T.A)%*%t(T.A) %*% NCI.R[[i]]
    
    RMSE_K[[1]][a,i] <- RMSE(X1[[i]],NCI.R[[i]],k=a)
  }}

png(file=paste0(pathI,"/NCI_MBPCA_RMSE1.png"), width = 777, height = 463)
ggplot() + geom_point(aes(seq(1,10),RMSE_K[[1]][,1],color="agilent")) + geom_line(aes(seq(1,10),RMSE_K[[1]][,1],color="agilent"))  + geom_point(aes(seq(1,10),RMSE_K[[1]][,2],color="hgu133")) + geom_line(aes(seq(1,10),RMSE_K[[1]][,2],color="hgu133"))  + geom_point(aes(seq(1,10),RMSE_K[[1]][,3],color="hgu133p2")) + geom_line(aes(seq(1,10),RMSE_K[[1]][,3],color="hgu133p2"))  + geom_point(aes(seq(1,10),RMSE_K[[1]][,4],color="hgu95")) + geom_line(aes(seq(1,10),RMSE_K[[1]][,4],color="hgu95")) +  labs(title="RMSE selon le nombre de composantes",x="Nombre de composantes",y="RMSE")  + theme_bw() + guides(color=guide_legend("Bloc :"))+ scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()

print(xtable(t(RMSE_K[[1]]),caption="Valeur du RMSE pour les différents blocs de variables explicatives selon le nombre de composantes - Modèle consensus PCA avec données réduites",digits=3,scientific=F), floating = getOption("xtable.floating", FALSE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_MBPCA_RMSE1.tex"))


## 1.2 Inertie expliquée ----------------------------------------

Inertie <- list()

Inertie[[1]] <- t(VE_PCA(NCI,modR$t,10)) * 100
colnames(Inertie[[1]]) <- names(NCI)

png(file=paste0(pathI,"/NCI_MBPCA_Inert1.png"), width = 777, height = 463)
ggplot() + geom_point(aes(seq(1,10),Inertie[[1]][,1],color="agilent")) + geom_line(aes(seq(1,10),Inertie[[1]][,1],color="agilent")) + geom_point(aes(seq(1,10),Inertie[[1]][,2],color="hgu133")) + geom_line(aes(seq(1,10),Inertie[[1]][,2],color="hgu133")) +  geom_point(aes(seq(1,10),Inertie[[1]][,3],color="hgu133p2")) + geom_line(aes(seq(1,10),Inertie[[1]][,3],color="hgu133p2")) +  geom_point(aes(seq(1,10),Inertie[[1]][,4],color="hgu95")) + geom_line(aes(seq(1,10),Inertie[[1]][,4],color="hgu95")) +  labs(title="Pourcentage de variance expliquée selon le nombre de composantes",x="Nombre de composantes",y="Inertie expliquée (en %)")  + theme_bw() + guides(color=guide_legend("Bloc :")) + scale_x_continuous(breaks=seq(1,10),labels=seq(1,10)) + theme(text = element_text(size=18),legend.position="bottom")
dev.off()


print(xtable(t(Inertie[[1]]),caption="Variance expliquée en pourcentage pour les différents blocs de variables explicatives selon le nombre de composantes",digits=4,scientific=F), floating = getOption("xtable.floating", FALSE), sanitize.rownames.function = bold,sanitize.colnames.function =bold,table.placement = getOption("xtable.table.placement", "H"),file=paste0(pathR,"/NCI_MBPCA_Inert1.tex"))

# 2. Modèle multibloc PCA ------------------------------------------

A = 3
mod <- mbpca(NCI_b, ncomp=A,method="globalScore",scale=TRUE,moa=FALSE,verbose=FALSE)


## 2.1 Projection des observations ---------------------------------------

png(file=paste0(pathI,"/NCI_MBPCA_TK12.png"), width = 777, height = 463)
Plot.TK(mod)
dev.off()
png(file=paste0(pathI,"/NCI_MBPCA_TK34.png"), width = 777, height = 463)
Plot.TK(mod,comp=c(1,3))
dev.off()


## 2.2 Proximité entre les blocs  ---------------------------------------

Tk.Super <- mod$t[,1:A]
colnames(Tk.Super) <- paste("Composante",seq(1,A))
rownames(Tk.Super) <- rownames(NCI[[1]])

DD <- matrix(0,ncol=6,nrow=53) 
Pos <- combn(seq(1,4),2)
CC <- numeric(6)
for(i in 1:A){
  for(j in 1:6){
    x <- Pos[,j]
    DD[,j] <- DD[,j] + (mod$tb[[x[1]]][,i] - mod$tb[[x[2]]][,i])^2
    CC[j] <- CC[j] + cor(mod$tb[[x[1]]][,i],mod$tb[[x[2]]][,i])
  }}


table(apply(sqrt(DD),1,which.min))

table(apply(sqrt(DD),1,which.max))

colMeans(DD)
CC/A


mod$X <- NCI
## 2.3 Proximité entre les variables des différents blocs  ----------------------------                                              
# 2.3.1 Heatmap des poids des variables
require(gplots) ; require(dendextend)
Ak <- NULL
for(i in 1:4){
  Ak.a <- NULL
  for(j in 1:A) {
    ord1 <- order(abs(mod$pb[[i]][,j]),decreasing = TRUE)
    data1 <- mod$pb[[i]][ord1,]
    Ak.a <- rbind(Ak.a,data1[1:10,])
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
png(file=paste0(pathI,"/MBPCA_ClusterAk.png"), width = 1500, height = 750)
plot(dend,main="Clustering hiérarchique sur base du poids des variables",horiz = FALSE,cex.lab=2, cex.axis=2, cex.main=4, cex.sub=2)
abline(h=0.25,lty=2)
legend("topright",legend=NCI_Name,col=cols_4,pch=19,cex=2.5)
dev.off()

groupe <- cutree(dend,k=4)
bloc <- as.factor(substr(rownames(Ak),1,1))
for(i in 1:4) print(table(as.factor(substr(rownames(Ak),1,1))[groupe==i])/table(groupe)[i])
table(groupe)

distAk <- as.matrix(dAK)
distAk[upper.tri(distAk,diag=TRUE)] <- 1

for(i in 1:5) {
  coord <- Nmin(distAk,i)
  print(paste(colnames(distAk)[coord[2]],rownames(distAk)[coord[1]],distAk[coord]))
}

# 2.3.2 Clustering variables
require(made4)
COR <- Var.Covar(mod,comp=A,type="cor")
rownames(COR) <- unlist(lapply(NCI,colnames))
ID <- as.factor(COR[,A+1])
COR <- COR[,1:A]
xx <- apply(COR,1,function(x) all(abs(x)>0.20))
distM <- as.matrix(as.dist(1-cor(t(COR[xx,]))))
distM[upper.tri(distM,diag=TRUE)] <- 2

colnames(COR) <- paste('Comp',1:A)
png(file=paste0(pathI,"/MBPCA_HeatCorr.png"), width = 1500, height = 750)
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