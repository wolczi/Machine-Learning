####### SEED ##########

set.seed(20210615)

#######################

sciezka_pliki <- "C:/Users/Przemek/Documents/R/praca dyplomowa"
setwd(sciezka_pliki)      
getwd()

daneNr1 <- read.delim( file = "E-GEOD-72056.sdrf.txt")
daneNr2 <- read.delim( file = "GSE72056_melanoma_single_cell.txt")
View(daneNr1[1:15, 1:15])
View(daneNr2[1:15, 1:15])

daneNr2_transpozycja <- t(daneNr2) # kolumny jako wiersze, wiersze jako kolumny
daneNr2_transpozycja <- as.data.frame(daneNr2_transpozycja)

View(daneNr2_transpozycja[1:15, 1:15])

# pobieranie danych z wszystkich kolumn pierwszego wiersza i robienie z nich nazw kolumn
# potem usuwanie wiersza pierwszego (tego z ktorego korzystalismy)
names(daneNr2_transpozycja) <- as.character(daneNr2_transpozycja[1,])
daneNr2_transpozycja <- daneNr2_transpozycja[-1,]

library(tidyverse)

daneNr2_transpozycja$`malignant(1=no,2=yes,0=unresolved)` <- as.integer(daneNr2_transpozycja$`malignant(1=no,2=yes,0=unresolved)`)

# summary(as.factor(daneNr1$Characteristics..classification.based.on.inferred.cnvs.)) # podsumowanie (liczebnosc poszczegolnych klasyfikacji)
summary(as.factor(daneNr2_transpozycja$`malignant(1=no,2=yes,0=unresolved)`)) 
# wybór tylko tych danych, ktÃ³re sklasyfikowane jako zlosliwe lub niezlosliwe
dane_ZlosliwyNiezlosliwy <- daneNr2_transpozycja[(daneNr2_transpozycja$`malignant(1=no,2=yes,0=unresolved)`== 1 |  daneNr2_transpozycja$`malignant(1=no,2=yes,0=unresolved)`== 2), ]
View(dane_ZlosliwyNiezlosliwy[1:15, 1:15])

dane_ZlosliwyNiezlosliwy <- dane_ZlosliwyNiezlosliwy[,-1] #usuwanie pierwszej kolumny
dane_ZlosliwyNiezlosliwy <- dane_ZlosliwyNiezlosliwy[,-2] #usuwanie drugiej kolumny

names(dane_ZlosliwyNiezlosliwy)[1] <- "malignant(1=no,2=yes)"
names(dane_ZlosliwyNiezlosliwy)[1]

# data.frame bez nazw wierszy
row.names(dane_ZlosliwyNiezlosliwy) <- NULL

View(dane_ZlosliwyNiezlosliwy[1:15, 1:15])

apply(dane_ZlosliwyNiezlosliwy[1,1:5], 2, class)
table(dane_ZlosliwyNiezlosliwy$`malignant(1=no,2=yes)`)
# 1    2 
# 3256 1257 

dane_ZlosliwyNiezlosliwy2 <- dane_ZlosliwyNiezlosliwy
dane_ZlosliwyNiezlosliwy2 <- as.data.frame(lapply(dane_ZlosliwyNiezlosliwy2, as.numeric))
View(dane_ZlosliwyNiezlosliwy2[1:15, 1:15])
apply(dane_ZlosliwyNiezlosliwy2[1,1:5], 2, class)

nieZlosliwy <- dane_ZlosliwyNiezlosliwy2 %>% filter(malignant.1.no.2.yes. == 1)
Zlosliwy <- dane_ZlosliwyNiezlosliwy2 %>% filter(malignant.1.no.2.yes. == 2)


zbilansowane.Dane <- NA
zbilansowane.Dane <- rbind(nieZlosliwy[sample(nrow(nieZlosliwy))[1:1000],], Zlosliwy[sample(nrow(Zlosliwy))[1:1000],])
table(zbilansowane.Dane$malignant.1.no.2.yes.)
# 1    2 
# 1000 1000

View(zbilansowane.Dane[1:15,1:15])

################################
# Ttest - pomniejszanie zbioru do 1000-ca cech

sciezka_funkcje <- "C:/Users/Przemek/Documents/R/stare pliki z R/bioinfa/Laboratorium 12/modules"
source(paste0(sciezka_funkcje,'/varImpTtest.R'))

Ttest.iter <- list()

for(i in 1:10){
  rows.dane <- sample.int(nrow(zbilansowane.Dane), size = round(nrow(zbilansowane.Dane)/3), replace = F)
  dane.train <- zbilansowane.Dane[-rows.dane,]
  dane.test <- zbilansowane.Dane[rows.dane,]
  Ttest.iter[[i]] <- varImpTtest(data = dane.train, p.adjust.methods = "bonferroni")
}

Ttest.1000=as.data.frame(matrix(nrow=1000,ncol=10))

for(i in 1:10){                      
  Ttest.1000[,i] <- Ttest.iter[[i]]$gene[1:1000]
  colnames(Ttest.1000)[i] <- paste0("ttest", i)
}

library(base)

rank <- data.frame("gene" = c(Ttest.1000[1,1]),
                   "points" = c(0))

rank$gene <- as.character(rank$gene)
rank$points <- as.integer(rank$points)

for(i in 1:1000) {
  for(j in 1:10){
    if(Ttest.1000[i,j] %in% rank$gene){
      rank[rank$gene==Ttest.1000[i,j],2] = (as.integer(rank[rank$gene==Ttest.1000[i,j],2]) + as.integer((1001-i)))
    }
    else{
      rank[nrow(rank) + 1,] = c(Ttest.1000[i,j],as.integer(1001-i))
    }
  }
}

rank$points <- as.integer(rank$points)
order.rank <- rank[order(rank$points, decreasing=TRUE),]


dane.1000.Cech <- zbilansowane.Dane[,c("malignant.1.no.2.yes.", order.rank$gene[1:1000])]
dane.1000.Cech$malignant.1.no.2.yes.[dane.1000.Cech$malignant.1.no.2.yes. == 1] <- 0
dane.1000.Cech$malignant.1.no.2.yes.[dane.1000.Cech$malignant.1.no.2.yes. == 2] <- 1
dane.1000.Cech$malignant.1.no.2.yes. <- as.character(dane.1000.Cech$malignant.1.no.2.yes.)
dane.1000.Cech$malignant.1.no.2.yes.<- as.factor(dane.1000.Cech$malignant.1.no.2.yes.)


ttest.auc.mcc.error <- data.frame("ttest.auc" = numeric(),
                                  "ttest.mcc" = numeric(),
                                  "ttest.error" = numeric())

mdfs.auc.mcc.error <- data.frame("mdfs.auc" = numeric(),
                                 "mdfs.mcc" = numeric(),
                                 "mdfs.error" = numeric())

rforest.auc.mcc.error <- data.frame("rforest.auc" = numeric(),
                                    "rforest.mcc" = numeric(),
                                    "rforest.error" = numeric())

boruta.auc.mcc.error <- data.frame("boruta.auc" = numeric(),
                                   "boruta.mcc" = numeric(),
                                   "boruta.error" = numeric())

mrmr.auc.mcc.error <- data.frame("mrmr.auc" = numeric(),
                                 "mrmr.mcc" = numeric(),
                                 "mrmr.error" = numeric())

relief.auc.mcc.error <- data.frame("relief.auc" = numeric(),
                                   "relief.mcc" = numeric(),
                                   "relief.error" = numeric())

################################
# wykonywanie metod na zbiorze 1000-ca cech 
# ze zbilansowanymi obserwacjami
#
#
# Ttest na zbiorze 1000-ca cech

library(pROC)
library(randomForest)
library(mltools)

Ttest.iter.1000 <- list()
ttest.rf.prob.list <- list()
ttest.rf.result.list <- list()
ttest.rf.correct.list <- list()

suma_auc = 0
suma_mcc = 0 
auc_tmp = 0 
mcc_tmp = 0 
suma = 0

for(i in 1:50){
  # selekcja cech
  rows.dane <- sample.int(nrow(dane.1000.Cech), size = round(nrow(dane.1000.Cech)/3), replace = F)
  dane.train.1000 <- dane.1000.Cech[-rows.dane,]
  dane.test.1000 <- dane.1000.Cech[rows.dane,]
  Ttest.iter.1000[[i]] <- varImpTtest(data = dane.train.1000, p.adjust.methods = "bonferroni")
  
  ttest.best100 <- Ttest.iter.1000[[i]]$gene[1:100]
  ttest.best100[101] <- "malignant.1.no.2.yes."
  
  # klasyfikacja
  dane.train.1000 <- dane.train.1000[, ttest.best100]
  dane.test.1000 <- dane.test.1000[, ttest.best100]
  
  model.rf <- randomForest(x=dane.train.1000[,-101], y = dane.train.1000[,101],
                           ntree = 500, do.trace = 100)
  rf.result <- predict(model.rf, newdata = dane.test.1000[,-101])
  rf.prob <- predict(model.rf, newdata = dane.test.1000[,-101], type = 'prob')[,2]
  
  # auc
  auc_tmp <- auc_roc(as.numeric(rf.prob)-1, as.numeric(dane.test.1000[,101])-1)
  suma_auc <- suma_auc + auc_tmp
  
  # mcc
  mcc_tmp = mcc(as.numeric(rf.result)-1, as.numeric(dane.test.1000[,101])-1)
  suma_mcc <- suma_mcc + mcc_tmp
  
  # error
  error.rf1 <- sum(dane.test.1000[,101] != rf.result)/length(rf.result)
  suma <- suma + error.rf1
  
  # auc, mcc, error ka¿dej iteracji
  ttest.auc.mcc.error[i,1] <- c(auc_tmp)
  ttest.auc.mcc.error[i,2] <- c(mcc_tmp)
  ttest.auc.mcc.error[i,3] <- c(error.rf1)
  
  ttest.rf.prob.list[[i]] <- rf.prob
  ttest.rf.result.list[[i]] <- rf.result
  ttest.rf.correct.list[[i]] <- dane.test.1000[,101]
}


auc.ttest <- suma_auc/i
mcc.ttest <- suma_mcc/i
error.ttest <- suma/i


###### rankingowanie cech ######


Ttest.100=as.data.frame(matrix(nrow=100,ncol=10))

for(i in 1:10){                      
  Ttest.100[,i] <- Ttest.iter.1000[[i]]$gene[1:100]
  colnames(Ttest.100)[i] <- paste0("ttest", i)
}

library(base)

rank.Ttest.100 <- data.frame("gene" = c(Ttest.100[1,1]),
                             "points" = c(0))

rank.Ttest.100$gene <- as.character(rank.Ttest.100$gene)
rank.Ttest.100$points <- as.integer(rank.Ttest.100$points)

for(i in 1:100) {
  for(j in 1:10){
    if(Ttest.100[i,j] %in% rank.Ttest.100$gene){
      rank.Ttest.100[rank.Ttest.100$gene==Ttest.100[i,j],2] = (as.integer(rank.Ttest.100[rank.Ttest.100$gene==Ttest.100[i,j],2]) + as.integer((101-i)))
    }
    else{
      rank.Ttest.100[nrow(rank.Ttest.100) + 1,] = c(Ttest.100[i,j],as.integer(101-i))
    }
  }
}

rank.Ttest.100$points <- as.integer(rank.Ttest.100$points)
order.rank.Ttest.100 <- rank.Ttest.100[order(rank.Ttest.100$points, decreasing=TRUE)[1:100],]


################################
# MDFS 


sciezka_funkcje <- "C:/Users/Przemek/Documents/R/stare pliki z R/bioinfa/Laboratorium 12/modules"
source(paste0(sciezka_funkcje,'/varImpMDFS1D.R'))

library(MDFS)

MDFS.iter.1000 <- list()
mdfs.rf.prob.list <- list()
mdfs.rf.result.list <- list()
mdfs.rf.correct.list <- list()

suma_auc = 0
suma_mcc = 0 
auc_tmp = 0 
mcc_tmp = 0 
suma = 0

for(i in 1:50){
  # selekcja cech
  rows.dane <- sample.int(nrow(dane.1000.Cech), size = round(nrow(dane.1000.Cech)/3), replace = F)
  dane.train.1000 <- dane.1000.Cech[-rows.dane,]
  dane.test.1000 <- dane.1000.Cech[rows.dane,]
  MDFS.iter.1000[[i]] <- varImpMDFS1D(data = dane.train.1000, p.adjust.methods = "bonferroni")

  mdfs.best100 <-  MDFS.iter.1000[[i]]$gene[1:100]
  mdfs.best100[101] <- "malignant.1.no.2.yes."
  
  # klasyfikacja
  dane.train.1000 <- dane.train.1000[, mdfs.best100]
  dane.test.1000 <- dane.test.1000[, mdfs.best100]
  
  model.rf <- randomForest(x=dane.train.1000[,-101], y = dane.train.1000[,101],
                           ntree = 500, do.trace = 100)
  rf.result <- predict(model.rf, newdata = dane.test.1000[,-101])
  rf.prob <- predict(model.rf, newdata = dane.test.1000[,-101], type = 'prob')[,2]
  
  # auc
  auc_tmp <- auc_roc(as.numeric(rf.prob)-1, as.numeric(dane.test.1000[,101])-1)
  suma_auc <- suma_auc + auc_tmp
  
  # mcc
  mcc_tmp = mcc(as.numeric(rf.result)-1, as.numeric(dane.test.1000[,101])-1)
  suma_mcc <- suma_mcc + mcc_tmp
  
  # error
  error.rf1 <- sum(dane.test.1000[,101] != rf.result)/length(rf.result)
  suma <- suma + error.rf1
  
  # auc, mcc, error ka¿dej iteracji
  mdfs.auc.mcc.error[i,1] <- c(auc_tmp)
  mdfs.auc.mcc.error[i,2] <- c(mcc_tmp)
  mdfs.auc.mcc.error[i,3] <- c(error.rf1)
  
  mdfs.rf.prob.list[[i]] <- rf.prob
  mdfs.rf.result.list[[i]] <- rf.result
  mdfs.rf.correct.list[[i]] <- dane.test.1000[,101]
  
}

auc.mdfs <- suma_auc/i
mcc.mdfs <- suma_mcc/i
error.mdfs <- suma/i

###### rankingowanie cech ######


MDFS.100=as.data.frame(matrix(nrow=100,ncol=10))

for(i in 1:10){                      
  MDFS.100[,i] <- MDFS.iter.1000[[i]]$gene[1:100]
  colnames(MDFS.100)[i] <- paste0("mdfs", i)
}

library(base)

rank.MDFS.100 <- data.frame("gene" = c(MDFS.100[1,1]),
                            "points" = c(0))

rank.MDFS.100$gene <- as.character(rank.MDFS.100$gene)
rank.MDFS.100$points <- as.integer(rank.MDFS.100$points)

for(i in 1:100) {
  for(j in 1:10){
    if(MDFS.100[i,j] %in% rank.MDFS.100$gene){
      rank.MDFS.100[rank.MDFS.100$gene==MDFS.100[i,j],2] = (as.integer(rank.MDFS.100[rank.MDFS.100$gene==MDFS.100[i,j],2]) + as.integer((101-i)))
    }
    else{
      rank.MDFS.100[nrow(rank.MDFS.100) + 1,] = c(MDFS.100[i,j],as.integer(101-i))
    }
  }
}

rank.MDFS.100$points <- as.integer(rank.MDFS.100$points)
order.rank.MDFS.100 <- rank.MDFS.100[order(rank.MDFS.100$points, decreasing=TRUE)[1:100],]

################################
# Boruta

library(Boruta)
library(dplyr)

Boruta.iter.1000 <- list()
boruta.rf.prob.list <- list()
boruta.rf.result.list <- list()
boruta.rf.correct.list <- list()

suma_auc = 0
suma_mcc = 0 
auc_tmp = 0 
mcc_tmp = 0 
suma = 0

for(i in 1:50){
  # selekcja cech
  rows.dane <- sample.int(nrow(dane.1000.Cech), size = round(nrow(dane.1000.Cech)/3), replace = F)
  dane.train.1000 <- dane.1000.Cech[-rows.dane,]
  dane.test.1000 <- dane.1000.Cech[rows.dane,]
  
  result.Boruta <- Boruta(malignant.1.no.2.yes.~., data = dane.train.1000, doTrace = 2)
  
  stats.Boruta <- attStats(result.Boruta)
  
  imp.Boruta <- stats.Boruta %>% filter(decision == "Confirmed")
  imp.Boruta <- imp.Boruta[order(-imp.Boruta$medianImp),]
  imp.Boruta <- imp.Boruta[-c(1,3,4,5,6)]
  imp.Boruta <- as.data.frame(cbind("gene" = rownames(imp.Boruta), "medianImp" = imp.Boruta$medianImp))
  
  Boruta.iter.1000[[i]] <- imp.Boruta
  
  boruta.best100 <- Boruta.iter.1000[[i]]$gene[1:100]
  boruta.best100[101] <- "malignant.1.no.2.yes."
  
  # klasyfikacja
  dane.train.1000 <- dane.train.1000[, boruta.best100]
  dane.test.1000 <- dane.test.1000[, boruta.best100]
  
  model.rf <- randomForest(x=dane.train.1000[,-101], y = dane.train.1000[,101],
                           ntree = 500, do.trace = 100)
  rf.result <- predict(model.rf, newdata = dane.test.1000[,-101])
  rf.prob <- predict(model.rf, newdata = dane.test.1000[,-101], type = 'prob')[,2]
  
  # auc
  auc_tmp <- auc_roc(as.numeric(rf.prob)-1, as.numeric(dane.test.1000[,101])-1)
  suma_auc <- suma_auc + auc_tmp
  
  # mcc
  mcc_tmp = mcc(as.numeric(rf.result)-1, as.numeric(dane.test.1000[,101])-1)
  suma_mcc <- suma_mcc + mcc_tmp
  
  # error
  error.rf1 <- sum(dane.test.1000[,101] != rf.result)/length(rf.result)
  suma <- suma + error.rf1
  
  # auc, mcc, error ka¿dej iteracji
  boruta.auc.mcc.error[i,1] <- c(auc_tmp)
  boruta.auc.mcc.error[i,2] <- c(mcc_tmp)
  boruta.auc.mcc.error[i,3] <- c(error.rf1)
  
  boruta.rf.prob.list[[i]] <- rf.prob
  boruta.rf.result.list[[i]] <- rf.result
  boruta.rf.correct.list[[i]] <- dane.test.1000[,101]
}


auc.boruta <- suma_auc/i
mcc.boruta <- suma_mcc/i
error.boruta <- suma/i


###### rankingowanie cech ######


Boruta.100=as.data.frame(matrix(nrow=100,ncol=10))

for(i in 1:10){                      
  Boruta.100[,i] <- Boruta.iter.1000[[i]]$gene[1:100]
  colnames(Boruta.100)[i] <- paste0("boruta", i)
}

library(base)

rank.Boruta.100 <- data.frame("gene" = c(Boruta.100[1,1]),
                              "points" = c(0))

rank.Boruta.100$gene <- as.character(rank.Boruta.100$gene)
rank.Boruta.100$points <- as.integer(rank.Boruta.100$points)

for(i in 1:100) {
  for(j in 1:10){
    if(Boruta.100[i,j] %in% rank.Boruta.100$gene){
      rank.Boruta.100[rank.Boruta.100$gene==Boruta.100[i,j],2] = (as.integer(rank.Boruta.100[rank.Boruta.100$gene==Boruta.100[i,j],2]) + as.integer((101-i)))
    }
    else{
      rank.Boruta.100[nrow(rank.Boruta.100) + 1,] = c(Boruta.100[i,j],as.integer(101-i))
    }
  }
}

rank.Boruta.100$points <- as.integer(rank.Boruta.100$points)
order.rank.Boruta.100 <- rank.Boruta.100[order(rank.Boruta.100$points, decreasing=TRUE)[1:100],]

################################
# mRmR
library(mRMRe)

mRmR.iter.1000 <- list()
mrmr.rf.prob.list <- list()
mrmr.rf.result.list <- list()
mrmr.rf.correct.list <- list()

suma_auc = 0
suma_mcc = 0 
auc_tmp = 0 
mcc_tmp = 0 
suma = 0

for(i in 1:50){
  # selekcja cech
  rows.dane <- sample.int(nrow(dane.1000.Cech), size = round(nrow(dane.1000.Cech)/3), replace = F)
  dane.train.1000 <- dane.1000.Cech[-rows.dane,]
  dane.test.1000 <- dane.1000.Cech[rows.dane,]
  
  data.For.mRmR <- mRMR.data(data = dane.train.1000[,-1])
  result.mRmR <- mRMR.classic(data = data.For.mRmR, target_indices = c(ncol(dane.train.1000[,-1])))
  imp.mRmR <-data.frame('importance'=result.mRmR@mi_matrix[nrow(result.mRmR@mi_matrix),])
  imp.mRmR$feature <- rownames(imp.mRmR)
  row.names(imp.mRmR) <- NULL
  imp.mRmR <- na.omit(imp.mRmR)
  imp.mRmR <- imp.mRmR[order(imp.mRmR$importance, decreasing=TRUE),]
  colnames(imp.mRmR)[2] <- "gene"
  
  mRmR.iter.1000[[i]] <- imp.mRmR
  
  mrmr.best100 <-  mRmR.iter.1000[[i]]$gene[1:100]
  mrmr.best100[101] <- "malignant.1.no.2.yes."
  
  # klasyfikacja
  dane.train.1000 <- dane.train.1000[, mrmr.best100]
  dane.test.1000 <- dane.test.1000[, mrmr.best100]
  
  model.rf <- randomForest(x=dane.train.1000[,-101], y = dane.train.1000[,101],
                           ntree = 500, do.trace = 100)
  rf.result <- predict(model.rf, newdata = dane.test.1000[,-101])
  rf.prob <- predict(model.rf, newdata = dane.test.1000[,-101], type = 'prob')[,2]
  
  # auc
  auc_tmp <- auc_roc(as.numeric(rf.prob)-1, as.numeric(dane.test.1000[,101])-1)
  suma_auc <- suma_auc + auc_tmp
  
  # mcc
  mcc_tmp = mcc(as.numeric(rf.result)-1, as.numeric(dane.test.1000[,101])-1)
  suma_mcc <- suma_mcc + mcc_tmp
  
  # error
  error.rf1 <- sum(dane.test.1000[,101] != rf.result)/length(rf.result)
  suma <- suma + error.rf1
  
  # auc, mcc, error ka¿dej iteracji
  mrmr.auc.mcc.error[i,1] <- c(auc_tmp)
  mrmr.auc.mcc.error[i,2] <- c(mcc_tmp)
  mrmr.auc.mcc.error[i,3] <- c(error.rf1)
  
  mrmr.rf.prob.list[[i]] <- rf.prob
  mrmr.rf.result.list[[i]] <- rf.result
  mrmr.rf.correct.list[[i]] <- dane.test.1000[,101]
  
}

auc.mrmr <- suma_auc/i
mcc.mrmr <- suma_mcc/i
error.mrmr <- suma/i


###### rankingowanie cech ######


mRmR.100=as.data.frame(matrix(nrow=100,ncol=10))

for(i in 1:10){                      
  mRmR.100[,i] <- mRmR.iter.1000[[i]]$gene[1:100]
  colnames(mRmR.100)[i] <- paste0("mrmr", i)
}

library(base)

rank.mRmR.100 <- data.frame("gene" = c(mRmR.100[1,1]),
                            "points" = c(0))

rank.mRmR.100$gene <- as.character(rank.mRmR.100$gene)
rank.mRmR.100$points <- as.integer(rank.mRmR.100$points)

for(i in 1:100) {
  for(j in 1:10){
    if(mRmR.100[i,j] %in% rank.mRmR.100$gene){
      rank.mRmR.100[rank.mRmR.100$gene==mRmR.100[i,j],2] = (as.integer(rank.mRmR.100[rank.mRmR.100$gene==mRmR.100[i,j],2]) + as.integer((101-i)))
    }
    else{
      rank.mRmR.100[nrow(rank.mRmR.100) + 1,] = c(mRmR.100[i,j],as.integer(101-i))
    }
  }
}

rank.mRmR.100$points <- as.integer(rank.mRmR.100$points)
order.rank.mRmR.100 <- rank.mRmR.100[order(rank.mRmR.100$points, decreasing=TRUE)[1:100],]

################################
# ReliefF

library(FSelector)

Relief.iter.1000 <- list()
relief.rf.prob.list <- list()
relief.rf.result.list <- list()
relief.rf.correct.list <- list()

suma_auc = 0
suma_mcc = 0 
auc_tmp = 0 
mcc_tmp = 0 
suma = 0

for(i in 1:50){
  # selekcja cech
  rows.dane <- sample.int(nrow(dane.1000.Cech), size = round(nrow(dane.1000.Cech)/3), replace = F)
  dane.train.1000 <- dane.1000.Cech[-rows.dane,]
  dane.test.1000 <- dane.1000.Cech[rows.dane,]
  
  imp.Relief <- relief(malignant.1.no.2.yes.~., dane.train.1000, neighbours.count = 3, sample.size = 5)
  
  subset <- cutoff.k(imp.Relief, 100)
  imp.Relief <- as.data.frame(subset)
  colnames(imp.Relief)[1] <- "gene"
  
  Relief.iter.1000[[i]] <- imp.Relief
  
  relief.best100 <-  Relief.iter.1000[[i]]$gene[1:100]
  relief.best100[101] <- "malignant.1.no.2.yes."
  
  # klasyfikacja 
  dane.train.1000 <- dane.train.1000[, relief.best100]
  dane.test.1000 <- dane.test.1000[, relief.best100]
  
  model.rf <- randomForest(x=dane.train.1000[,-101], y = dane.train.1000[,101],
                           ntree = 500, do.trace = 100)
  rf.result <- predict(model.rf, newdata = dane.test.1000[,-101])
  rf.prob <- predict(model.rf, newdata = dane.test.1000[,-101], type = 'prob')[,2]
  
  # auc
  auc_tmp <- auc_roc(as.numeric(rf.prob)-1, as.numeric(dane.test.1000[,101])-1)
  suma_auc <- suma_auc + auc_tmp
  
  # mcc
  mcc_tmp = mcc(as.numeric(rf.result)-1, as.numeric(dane.test.1000[,101])-1)
  suma_mcc <- suma_mcc + mcc_tmp
  
  # error
  error.rf1 <- sum(dane.test.1000[,101] != rf.result)/length(rf.result)
  suma <- suma + error.rf1
  
  # auc, mcc, error ka¿dej iteracji
  relief.auc.mcc.error[i,1] <- c(auc_tmp)
  relief.auc.mcc.error[i,2] <- c(mcc_tmp)
  relief.auc.mcc.error[i,3] <- c(error.rf1)
  
  relief.rf.prob.list[[i]] <- rf.prob
  relief.rf.result.list[[i]] <- rf.result
  relief.rf.correct.list[[i]] <- dane.test.1000[,101]
}

auc.relief <- suma_auc/i
mcc.relief <- suma_mcc/i
error.relief <- suma/i

###### rankingowanie cech ######


Relief.100=as.data.frame(matrix(nrow=100,ncol=10))

for(i in 1:10){                      
  Relief.100[,i] <- Relief.iter.1000[[i]]$gene[1:100]
  colnames(Relief.100)[i] <- paste0("relief", i)
}

library(base)

rank.Relief.100 <- data.frame("gene" = c(Relief.100[1,1]),
                            "points" = c(0))

rank.Relief.100$gene <- as.character(rank.Relief.100$gene)
rank.Relief.100$points <- as.integer(rank.Relief.100$points)

for(i in 1:100) {
  for(j in 1:10){
    if(Relief.100[i,j] %in% rank.Relief.100$gene){
      rank.Relief.100[rank.Relief.100$gene==Relief.100[i,j],2] = (as.integer(rank.Relief.100[rank.Relief.100$gene==Relief.100[i,j],2]) + as.integer((101-i)))
    }
    else{
      rank.Relief.100[nrow(rank.Relief.100) + 1,] = c(Relief.100[i,j],as.integer(101-i))
    }
  }
}

rank.Relief.100$points <- as.integer(rank.Relief.100$points)
order.rank.Relief.100 <- rank.Relief.100[order(rank.Relief.100$points, decreasing=TRUE)[1:100],]

################################
# Random Forest Importance

library(FSelector)

randomForest.iter.1000 <- list()
rforest.rf.prob.list <- list()
rforest.rf.result.list <- list()
rforest.rf.correct.list <- list()

suma_auc = 0
suma_mcc = 0 
auc_tmp = 0 
mcc_tmp = 0 
suma = 0

for(i in 1:50){
  # selekcja cech
  rows.dane <- sample.int(nrow(dane.1000.Cech), size = round(nrow(dane.1000.Cech)/3), replace = F)
  dane.train.1000 <- dane.1000.Cech[-rows.dane,]
  dane.test.1000 <- dane.1000.Cech[rows.dane,]
  
  imp.RF <- random.forest.importance(malignant.1.no.2.yes.~., dane.train.1000, importance.type = 1)
  imp.RandomForest <- cutoff.k(imp.RF, 100)
  imp.RandomForest  <- as.data.frame(imp.RandomForest)
  colnames(imp.RandomForest)[1] <- "gene"
  
  randomForest.iter.1000[[i]] <- imp.RandomForest
  
  randomforest.best100 <-  randomForest.iter.1000[[i]]$gene[1:100]
  randomforest.best100[101] <- "malignant.1.no.2.yes."
  
  # klasyfikacja
  dane.train.1000 <- dane.train.1000[, randomforest.best100]
  dane.test.1000 <- dane.test.1000[, randomforest.best100]
  
  model.rf <- randomForest(x=dane.train.1000[,-101], y = dane.train.1000[,101],
                           ntree = 500, do.trace = 100)
  rf.result <- predict(model.rf, newdata = dane.test.1000[,-101])
  rf.prob <- predict(model.rf, newdata = dane.test.1000[,-101], type = 'prob')[,2]
  
  # auc
  auc_tmp <- auc_roc(as.numeric(rf.prob)-1, as.numeric(dane.test.1000[,101])-1)
  suma_auc <- suma_auc + auc_tmp
  
  # mcc
  mcc_tmp = mcc(as.numeric(rf.result)-1, as.numeric(dane.test.1000[,101])-1)
  suma_mcc <- suma_mcc + mcc_tmp
  
  # error
  error.rf1 <- sum(dane.test.1000[,101] != rf.result)/length(rf.result)
  suma <- suma + error.rf1
  
  # auc, mcc, error ka¿dej iteracji
  rforest.auc.mcc.error[i,1] <- c(auc_tmp)
  rforest.auc.mcc.error[i,2] <- c(mcc_tmp)
  rforest.auc.mcc.error[i,3] <- c(error.rf1)
  
  rforest.rf.prob.list[[i]] <- rf.prob
  rforest.rf.result.list[[i]] <- rf.result
  rforest.rf.correct.list[[i]] <- dane.test.1000[,101]
}

auc.randomforest <- suma_auc/i
mcc.randomforest <- suma_mcc/i
error.randomforest <- suma/i

###### rankingowanie cech ######


randomForest.100=as.data.frame(matrix(nrow=100,ncol=10))

for(i in 1:10){                      
  randomForest.100[,i] <- randomForest.iter.1000[[i]]$gene[1:100]
  colnames(randomForest.100)[i] <- paste0("rforest", i)
}

library(base)

rank.randomForest.100 <- data.frame("gene" = c(randomForest.100[1,1]),
                                    "points" = c(0))

rank.randomForest.100$gene <- as.character(rank.randomForest.100$gene)
rank.randomForest.100$points <- as.integer(rank.randomForest.100$points)

for(i in 1:100) {
  for(j in 1:10){
    if(randomForest.100[i,j] %in% rank.randomForest.100$gene){
      rank.randomForest.100[rank.randomForest.100$gene==randomForest.100[i,j],2] = (as.integer(rank.randomForest.100[rank.randomForest.100$gene==randomForest.100[i,j],2]) + as.integer((101-i)))
    }
    else{
      rank.randomForest.100[nrow(rank.randomForest.100) + 1,] = c(randomForest.100[i,j],as.integer(101-i))
    }
  }
}

rank.randomForest.100$points <- as.integer(rank.randomForest.100$points)
order.rank.randomForest.100 <- rank.randomForest.100[order(rank.randomForest.100$points, decreasing=TRUE)[1:100],]

#############################################
# porównanie klasyfikacji wszystkich metod

classi.all <- data.frame("ttest" = c("auc" = auc.ttest, 
                                     "mcc" = mcc.ttest, 
                                     "error" = error.ttest),
                         "mdfs" = c("auc" = auc.mdfs, 
                                     "mcc" = mcc.mdfs, 
                                     "error" = error.mdfs),
                         "rforest" = c("auc" = auc.randomforest, 
                                    "mcc" = mcc.randomforest, 
                                    "error" = error.randomforest),
                         "boruta" = c("auc" = auc.boruta, 
                                    "mcc" = mcc.boruta, 
                                    "error" = error.boruta),
                         "mrmr" = c("auc" = auc.mrmr, 
                                    "mcc" = mcc.mrmr, 
                                    "error" = error.mrmr),
                         "relief" = c("auc" = auc.relief, 
                                    "mcc" = mcc.relief, 
                                    "error" = error.relief)
                         )

#############################################
# wykresy cz. 1


# ttest
box.auc.ttest <- ttest.auc.mcc.error$ttest.auc
box.mcc.ttest <- ttest.auc.mcc.error$ttest.mcc
box.error.ttest <- ttest.auc.mcc.error$ttest.error
colMeans(ttest.auc.mcc.error)

# mdfs
box.auc.mdfs <- mdfs.auc.mcc.error$mdfs.auc
box.mcc.mdfs <- mdfs.auc.mcc.error$mdfs.mcc
box.error.mdfs <- mdfs.auc.mcc.error$mdfs.error
colMeans(mdfs.auc.mcc.error)


# boruta
box.auc.boruta <- boruta.auc.mcc.error$boruta.auc
box.mcc.boruta <- boruta.auc.mcc.error$boruta.mcc
box.error.boruta <- boruta.auc.mcc.error$boruta.error
colMeans(boruta.auc.mcc.error)

# mrmr
box.auc.mrmr <- mrmr.auc.mcc.error$mrmr.auc
box.mcc.mrmr <- mrmr.auc.mcc.error$mrmr.mcc
box.error.mrmr <- mrmr.auc.mcc.error$mrmr.error
colMeans(mrmr.auc.mcc.error)

# relief
box.auc.relief <- relief.auc.mcc.error$relief.auc
box.mcc.relief <- relief.auc.mcc.error$relief.mcc
box.error.relief <- relief.auc.mcc.error$relief.error
colMeans(relief.auc.mcc.error)

# rforest
box.auc.rforest <- rforest.auc.mcc.error$rforest.auc
box.mcc.rforest <- rforest.auc.mcc.error$rforest.mcc
box.error.rforest <- rforest.auc.mcc.error$rforest.error
colMeans(rforest.auc.mcc.error)

boxplot(box.auc.ttest, box.auc.mdfs, box.auc.boruta, box.auc.mrmr, box.auc.relief, box.auc.rforest,
        col=c("green","blue","orange","grey", "red", "purple"), 
        names =c("t test","mdfs", "boruta", "mrmr", "relief", "rforest"), ylab ="AUC", notch = T,
        cex.main=2.0, cex.lab=1.5, cex.axis=1.5)

boxplot(box.mcc.ttest, box.mcc.mdfs, box.mcc.boruta, box.mcc.mrmr, box.mcc.relief, box.mcc.rforest,
        col=c("green","blue","orange","grey", "red", "purple"),
        names =c("t test","mdfs", "boruta", "mrmr", "relief", "rforest"), ylab ="MCC", notch = T,
        cex.main=2.0, cex.lab=1.5, cex.axis=1.5)

boxplot(box.error.ttest, box.error.mdfs, box.error.boruta, box.error.mrmr, box.error.relief, box.error.rforest,
        col=c("green","blue","orange","grey", "red", "purple"),
        names =c("t test","mdfs", "boruta", "mrmr", "relief", "rforest"), ylab ="Error", notch = T,
        cex.main=2.0, cex.lab=1.5, cex.axis=1.5)

#############################################
# wykresy cz. 2

plotroc<-function(y, x=(1:length(y)), type='l', xlim=c(1,0), xlab='Specificity', ylab='Sensitivity', identity=T, ...) {
  yo<-y[order(x)]
  plot(cumsum(1-yo)/sum(1-yo),1-cumsum(yo)/sum(yo),xlim=xlim,type=type,xlab=xlab,ylab=ylab,...)
  if (identity) lines(c(0,1),c(1,0),lty=3,col='gray')
}

pointsroc<-function(y, x=(1:length(y)), type='l', xlim=c(1,0),...) {
  yo<-y[order(x)]
  points(cumsum(1-yo)/sum(1-yo),1-cumsum(yo)/sum(yo),type=type,...)
}


# ttest
probs <-vector()
correct <-vector()
for (i in 1:50)
{
  probs <-c(probs, as.numeric(ttest.rf.prob.list[[i]]))
  correct <-c(correct, ttest.rf.prob.list[[i]]) 
}

plotroc(y= as.numeric(correct), x=probs, col='red', cex.main=2, cex.lab=1.5, cex.axis=2)


# mdfs
probs <-vector()
correct <-vector()
for (i in 1:50)
{
  probs <-c(probs, as.numeric(mdfs.rf.prob.list[[i]]))
  correct <-c(correct, mdfs.rf.prob.list[[i]])
}

pointsroc(y= as.numeric(correct), x=probs, col='green')


# boruta
probs <-vector()
correct <-vector()
for (i in 1:50)
{
  probs <-c(probs, as.numeric(boruta.rf.prob.list[[i]]))
  correct <-c(correct, boruta.rf.prob.list[[i]])
}

pointsroc(y= as.numeric(correct), x=probs, col='blue')


# mrmr
probs <-vector()
correct <-vector()
for (i in 1:50)
{
  probs <-c(probs, as.numeric(mrmr.rf.prob.list[[i]]))
  correct <-c(correct, mrmr.rf.prob.list[[i]])
}


pointsroc(y= as.numeric(correct), x=probs, col='magenta')


# relief
probs <-vector()
correct <-vector()
for (i in 1:50)
{
  probs <-c(probs, as.numeric(relief.rf.prob.list[[i]]))
  correct <-c(correct, relief.rf.prob.list[[i]])
}

pointsroc(y= as.numeric(correct), x=probs, col='gray')


# rforest
probs <-vector()
correct <-vector()
for (i in 1:50)
{
  probs <-c(probs, as.numeric(rforest.rf.prob.list[[i]]))
  correct <-c(correct, rforest.rf.prob.list[[i]])
}

pointsroc(y= as.numeric(correct), x=probs, col='dark green')


# legenda
legend( 0.25, 0.65, legend = c("ttest","mdfs", "boruta", "mrmr", "relief", "rforest"), pch = '___',
        col = c("red", "green", "blue", "magenta", "gray", "dark green"), cex=2 )


#############################################
# porównanie rankingów selekcji cech wszystkich metod

all.methods.rank.100 <- as.data.frame(cbind("Ttest" = order.rank.Ttest.100$gene, 
                               "MDFS" = order.rank.MDFS.100$gene, 
                               "Boruta" = order.rank.Boruta.100$gene, 
                               "mRmR" = order.rank.mRmR.100$gene, 
                               "ReliefF" = order.rank.Relief.100$gene, 
                               "randomForest" = order.rank.randomForest.100$gene))

#################################
# PUNKTACJA CECH Z WSZYSTKICH METOD
# RANKING

library(base)

rank <- data.frame("gene" = c(all.methods.rank.100[1,1]),
                   "points" = c(0))

rank$gene <- as.character(rank$gene)
rank$points <- as.integer(rank$points)

for(i in 1:100) {
  for(j in 1:6){
    if(all.methods.rank.100[i,j] %in% rank$gene){
      rank[rank$gene==all.methods.rank.100[i,j],2] = (as.integer(rank[rank$gene==all.methods.rank.100[i,j],2]) + as.integer((101-i)))
    }
    else{
      rank[nrow(rank) + 1,] = c(all.methods.rank.100[i,j],as.integer(101-i))
    }
  }
}

rank$points <- as.integer(rank$points)
order.rank <- rank[order(rank$points, decreasing=TRUE),]
