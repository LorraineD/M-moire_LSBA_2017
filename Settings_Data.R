# Chargement des données ----------------------------------

load("Data/Class_mat.Rdata")
load("Data/ClassY.Rdata")
load("Data/NCI.Rdata")

for(i in 1:4){
  colnames(NCI[[i]]) <- gsub("/.*", '',  colnames(NCI[[i]]))
  colnames(NCI[[i]]) <- paste0(i,"-",colnames(NCI[[i]]))
}

ClassY.Num <- numeric(length(ClassY))

for(i in 1:length(ClassY)){
  if(ClassY[i]=="CN") ClassY.Num[i] <- 1 
  if(ClassY[i]=="CO") ClassY.Num[i] <- 2
  if(ClassY[i]=="LE") ClassY.Num[i] <- 3
  if(ClassY[i]=="ME") ClassY.Num[i] <- 4
  if(ClassY[i]=="LC") ClassY.Num[i] <- 5
  if(ClassY[i]=="OV") ClassY.Num[i] <- 6
  if(ClassY[i]=="RE") ClassY.Num[i] <- 7}

Class.C <- scale(Class.mat,scale=FALSE)

n <- nrow(Class.mat)
q <- ncol(NCI[[1]]) + ncol(NCI[[2]]) + ncol(NCI[[3]]) +ncol(NCI[[4]])
nb <- length(NCI)

X.NR <- cbind(NCI[[1]],NCI[[2]],NCI[[3]],NCI[[4]])
NCI.R <- lapply(NCI,scale,scale=T)
X.R <-  cbind(NCI.R[[1]],NCI.R[[2]],NCI.R[[3]],NCI.R[[4]])

NCI_Name <- names(NCI)


N0 <- c(1,ncol(NCI[[1]])+1,ncol(NCI[[1]]) + ncol(NCI[[2]])+1,ncol(NCI[[1]]) + ncol(NCI[[2]])+ncol(NCI[[3]])+1)
NT <- c(ncol(NCI[[1]]),ncol(NCI[[1]]) + ncol(NCI[[2]]),ncol(NCI[[1]]) + ncol(NCI[[2]])+ncol(NCI[[3]]),ncol(NCI[[1]]) + ncol(NCI[[2]])+ncol(NCI[[3]])+ncol(NCI[[4]]))

# Packages et options ---------------------------------------------------

require(ggplot2) ; require(xtable) ; require(ggrepel)
options(xtable.size="\\fontsize{8pt}{9pt}\\selectfont")
source("Fonctions/Fonction_Qualite.R")
source("C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Code/Graphe_Comp.R")

# Chemin d'accès -  table 
pathR <- "C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Mémoire_Ecriture/Table"

# Chemin d'accès - image
pathI <- "C:/Users/Lorraine/Documents/BSTA2/2ème Master/Mémoire/Mémoire_Ecriture/Image"

bold <- function(x){paste0('{\\bfseries ', x, '}')}
bold.somerows <-  function(x) gsub('BOLD(.*)',paste('\\\\textbf{\\1','}'),x)