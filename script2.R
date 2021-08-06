library(kernlab)
data(spam)

####### SEED ##########
set.seed(20210615)

#######################
spam$type <- as.character(spam$type)
spam$type[spam$type == "spam"] <- 1
spam$type[spam$type == "nonspam"] <- 0
spam$type <- as.factor(spam$type)


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
# T-Test

sciezka_funkcje <- "C:/Users/Przemek/Documents/R/stare pliki z R/bioinfa/Laboratorium 12/modules"
source(paste0(sciezka_funkcje,'/varImpTtestSPAM.R'))

library(pROC)
library(randomForest)
library(mltools)

Ttest.iter <- list()
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
  rows.dane <- sample.int(nrow(spam), size = round(nrow(spam)/3), replace = F)
  dane.train <- spam[-rows.dane,]
  dane.test <- spam[rows.dane,]
  Ttest.iter[[i]] <- varImpTtestSPAM(data = dane.train, p.adjust.methods = "bonferroni")
  
  ttest.best10 <- Ttest.iter[[i]]$word[1:10]
  ttest.best10[11] <- "type"
  
  # klasyfikacja
  dane.train <- dane.train[, ttest.best10]
  dane.test <- dane.test[, ttest.best10]
  
  model.rf <- randomForest(x=dane.train[,-11], y = dane.train[,11],
                           ntree = 500, do.trace = 100)
  rf.result <- predict(model.rf, newdata = dane.test[,-11])
  rf.prob <- predict(model.rf, newdata = dane.test[,-11], type = 'prob')[,2]
  
  # auc
  auc_tmp <- auc_roc(as.numeric(rf.prob)-1, as.numeric(dane.test[,11])-1)
  suma_auc <- suma_auc + auc_tmp
  
  # mcc
  mcc_tmp = mcc(as.numeric(rf.result)-1, as.numeric(dane.test[,11])-1)
  suma_mcc <- suma_mcc + mcc_tmp
  
  # error
  error.rf1 <- sum(dane.test[,11] != rf.result)/length(rf.result)
  suma <- suma + error.rf1

  # auc, mcc, error ka¿dej iteracji
  ttest.auc.mcc.error[i,1] <- c(auc_tmp)
  ttest.auc.mcc.error[i,2] <- c(mcc_tmp)
  ttest.auc.mcc.error[i,3] <- c(error.rf1)
  
  ttest.rf.prob.list[[i]] <- rf.prob
  ttest.rf.result.list[[i]] <- rf.result
  ttest.rf.correct.list[[i]] <- dane.test[,11]
}


auc.ttest <- suma_auc/i
mcc.ttest <- suma_mcc/i
error.ttest <- suma/i

###### rankingowanie cech ######


Ttest.10=as.data.frame(matrix(nrow=10,ncol=10))

for(i in 1:50){                      
  Ttest.10[,i] <- Ttest.iter[[i]]$word[1:10]
  colnames(Ttest.10)[i] <- paste0("ttest", i)
}

library(base)

rank.Ttest.10 <- data.frame("word" = c(Ttest.10[1,1]),
                             "points" = c(0))

rank.Ttest.10$word <- as.character(rank.Ttest.10$word)
rank.Ttest.10$points <- as.integer(rank.Ttest.10$points)

for(i in 1:10) {
  for(j in 1:50){
    if(Ttest.10[i,j] %in% rank.Ttest.10$word){
      rank.Ttest.10[rank.Ttest.10$word==Ttest.10[i,j],2] = (as.integer(rank.Ttest.10[rank.Ttest.10$word==Ttest.10[i,j],2]) + as.integer((11-i)))
    }
    else{
      rank.Ttest.10[nrow(rank.Ttest.10) + 1,] = c(Ttest.10[i,j],as.integer(11-i))
    }
  }
}

rank.Ttest.10$points <- as.integer(rank.Ttest.10$points)
order.rank.Ttest.10 <- rank.Ttest.10[order(rank.Ttest.10$points, decreasing=TRUE)[1:10],]

################################
# MDFS 


sciezka_funkcje <- "C:/Users/Przemek/Documents/R/stare pliki z R/bioinfa/Laboratorium 12/modules"
source(paste0(sciezka_funkcje,'/varImpMDFS1Dspam.R'))

library(MDFS)

MDFS.iter <- list()
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
  rows.dane <- sample.int(nrow(spam), size = round(nrow(spam)/3), replace = F)
  dane.train <- spam[-rows.dane,]
  dane.test <- spam[rows.dane,]
  MDFS.iter[[i]] <- varImpMDFS1Dspam(data = dane.train, p.adjust.methods = "bonferroni")
  
  mdfs.best10 <-  MDFS.iter[[i]]$word[1:10]
  mdfs.best10[11] <- "type"
  
  # klasyfikacja
  dane.train <- dane.train[, mdfs.best10]
  dane.test <- dane.test[, mdfs.best10]
  
  model.rf <- randomForest(x=dane.train[,-11], y = dane.train[,11],
                           ntree = 500, do.trace = 100)
  rf.result <- predict(model.rf, newdata = dane.test[,-11])
  rf.prob <- predict(model.rf, newdata = dane.test[,-11], type = 'prob')[,2]
  
  # auc
  auc_tmp <- auc_roc(as.numeric(rf.prob)-1, as.numeric(dane.test[,11])-1)
  suma_auc <- suma_auc + auc_tmp
 
  # mcc
  mcc_tmp = mcc(as.numeric(rf.result)-1, as.numeric(dane.test[,11])-1)
  suma_mcc <- suma_mcc + mcc_tmp
   
  # error
  error.rf1 <- sum(dane.test[,11] != rf.result)/length(rf.result)
  suma <- suma + error.rf1

  # auc, mcc, error ka¿dej iteracji
  mdfs.auc.mcc.error[i,1] <- c(auc_tmp)
  mdfs.auc.mcc.error[i,2] <- c(mcc_tmp)
  mdfs.auc.mcc.error[i,3] <- c(error.rf1)
  
  mdfs.rf.prob.list[[i]] <- rf.prob
  mdfs.rf.result.list[[i]] <- rf.result
  mdfs.rf.correct.list[[i]] <- dane.test[,11]
}

auc.mdfs <- suma_auc/i
mcc.mdfs <- suma_mcc/i
error.mdfs <- suma/i

###### rankingowanie cech ######


MDFS.10=as.data.frame(matrix(nrow=10,ncol=10))

for(i in 1:50){                      
  MDFS.10[,i] <- MDFS.iter[[i]]$word[1:10]
  colnames(MDFS.10)[i] <- paste0("mdfs", i)
}

library(base)

rank.MDFS.10 <- data.frame("word" = c(MDFS.10[1,1]),
                            "points" = c(0))

rank.MDFS.10$word <- as.character(rank.MDFS.10$word)
rank.MDFS.10$points <- as.integer(rank.MDFS.10$points)

for(i in 1:10) {
  for(j in 1:50){
    if(MDFS.10[i,j] %in% rank.MDFS.10$word){
      rank.MDFS.10[rank.MDFS.10$word==MDFS.10[i,j],2] = (as.integer(rank.MDFS.10[rank.MDFS.10$word==MDFS.10[i,j],2]) + as.integer((11-i)))
    }
    else{
      rank.MDFS.10[nrow(rank.MDFS.10) + 1,] = c(MDFS.10[i,j],as.integer(11-i))
    }
  }
}

rank.MDFS.10$points <- as.integer(rank.MDFS.10$points)
order.rank.MDFS.10 <- rank.MDFS.10[order(rank.MDFS.10$points, decreasing=TRUE)[1:10],]

################################
# Boruta

library(Boruta)
library(dplyr)

Boruta.iter <- list()
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
  rows.dane <- sample.int(nrow(spam), size = round(nrow(spam)/3), replace = F)
  dane.train <- spam[-rows.dane,]
  dane.test <- spam[rows.dane,]
  
  result.Boruta <- Boruta(type~., data = dane.train, doTrace = 2)
  
  stats.Boruta <- attStats(result.Boruta)
  
  imp.Boruta <- stats.Boruta %>% filter(decision == "Confirmed")
  imp.Boruta <- imp.Boruta[order(-imp.Boruta$medianImp),]
  imp.Boruta <- imp.Boruta[-c(1,3,4,5,6)]
  imp.Boruta <- as.data.frame(cbind("word" = rownames(imp.Boruta), "medianImp" = imp.Boruta$medianImp))
  
  Boruta.iter[[i]] <- imp.Boruta
  
  boruta.best10 <- Boruta.iter[[i]]$word[1:10]
  boruta.best10[11] <- "type"
  
  # klasyfikacja
  dane.train <- dane.train[, boruta.best10]
  dane.test <- dane.test[, boruta.best10]
  
  model.rf <- randomForest(x=dane.train[,-11], y = dane.train[,11],
                           ntree = 500, do.trace = 100)
  rf.result <- predict(model.rf, newdata = dane.test[,-11])
  rf.prob <- predict(model.rf, newdata = dane.test[,-11], type = 'prob')[,2]
  
  # auc
  auc_tmp <- auc_roc(as.numeric(rf.prob)-1, as.numeric(dane.test[,11])-1)
  suma_auc <- suma_auc + auc_tmp
  
  # mcc
  mcc_tmp = mcc(as.numeric(rf.result)-1, as.numeric(dane.test[,11])-1)
  suma_mcc <- suma_mcc + mcc_tmp

  # error
  error.rf1 <- sum(dane.test[,11] != rf.result)/length(rf.result)
  suma <- suma + error.rf1
  
  # auc, mcc, error ka¿dej iteracji
  boruta.auc.mcc.error[i,1] <- c(auc_tmp)
  boruta.auc.mcc.error[i,2] <- c(mcc_tmp)
  boruta.auc.mcc.error[i,3] <- c(error.rf1)
  
  boruta.rf.prob.list[[i]] <- rf.prob
  boruta.rf.result.list[[i]] <- rf.result
  boruta.rf.correct.list[[i]] <- dane.test[,11]
}


auc.boruta <- suma_auc/i
mcc.boruta <- suma_mcc/i
error.boruta <- suma/i

###### rankingowanie cech ######


Boruta.10=as.data.frame(matrix(nrow=10,ncol=10))

for(i in 1:50){                      
  Boruta.10[,i] <- Boruta.iter[[i]]$word[1:10]
  colnames(Boruta.10)[i] <- paste0("boruta", i)
}

library(base)

rank.Boruta.10 <- data.frame("word" = c(Boruta.10[1,1]),
                              "points" = c(0))

rank.Boruta.10$word <- as.character(rank.Boruta.10$word)
rank.Boruta.10$points <- as.integer(rank.Boruta.10$points)

for(i in 1:10) {
  for(j in 1:50){
    if(Boruta.10[i,j] %in% rank.Boruta.10$word){
      rank.Boruta.10[rank.Boruta.10$word==Boruta.10[i,j],2] = (as.integer(rank.Boruta.10[rank.Boruta.10$word==Boruta.10[i,j],2]) + as.integer((11-i)))
    }
    else{
      rank.Boruta.10[nrow(rank.Boruta.10) + 1,] = c(Boruta.10[i,j],as.integer(11-i))
    }
  }
}

rank.Boruta.10$points <- as.integer(rank.Boruta.10$points)
order.rank.Boruta.10 <- rank.Boruta.10[order(rank.Boruta.10$points, decreasing=TRUE)[1:10],]

################################
# mRmR
library(mRMRe)

mRmR.iter <- list()
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
  rows.dane <- sample.int(nrow(spam), size = round(nrow(spam)/3), replace = F)
  dane.train <- spam[-rows.dane,]
  dane.test <- spam[rows.dane,]
  
  data.For.mRmR <- mRMR.data(data = dane.train[,-58])
  result.mRmR <- mRMR.classic(data = data.For.mRmR, target_indices = c(ncol(dane.train[,-58])))
  imp.mRmR <-data.frame('importance'=result.mRmR@mi_matrix[nrow(result.mRmR@mi_matrix),])
  imp.mRmR$feature <- rownames(imp.mRmR)
  row.names(imp.mRmR) <- NULL
  imp.mRmR <- na.omit(imp.mRmR)
  imp.mRmR <- imp.mRmR[order(imp.mRmR$importance, decreasing=TRUE),]
  colnames(imp.mRmR)[2] <- "word"
  
  mRmR.iter[[i]] <- imp.mRmR
  
  mrmr.best10 <-  mRmR.iter[[i]]$word[1:10]
  mrmr.best10[11] <- "type"
  
  # klasyfikacja
  dane.train <- dane.train[, mrmr.best10]
  dane.test <- dane.test[, mrmr.best10]
  
  model.rf <- randomForest(x=dane.train[,-11], y = dane.train[,11],
                           ntree = 500, do.trace = 100)
  rf.result <- predict(model.rf, newdata = dane.test[,-11])
  rf.prob <- predict(model.rf, newdata = dane.test[,-11], type = 'prob')[,2]
  
  # auc
  auc_tmp <- auc_roc(as.numeric(rf.prob)-1, as.numeric(dane.test[,11])-1)
  suma_auc <- suma_auc + auc_tmp
  
  # mcc
  mcc_tmp = mcc(as.numeric(rf.result)-1, as.numeric(dane.test[,11])-1)
  suma_mcc <- suma_mcc + mcc_tmp
  
  # error
  error.rf1 <- sum(dane.test[,11] != rf.result)/length(rf.result)
  suma <- suma + error.rf1
  
  # auc, mcc, error ka¿dej iteracji
  mrmr.auc.mcc.error[i,1] <- c(auc_tmp)
  mrmr.auc.mcc.error[i,2] <- c(mcc_tmp)
  mrmr.auc.mcc.error[i,3] <- c(error.rf1)
  
  mrmr.rf.prob.list[[i]] <- rf.prob
  mrmr.rf.result.list[[i]] <- rf.result
  mrmr.rf.correct.list[[i]] <- dane.test[,11]
}

auc.mrmr <- suma_auc/i
mcc.mrmr <- suma_mcc/i
error.mrmr <- suma/i

###### rankingowanie cech ######


mRmR.10=as.data.frame(matrix(nrow=10,ncol=10))

for(i in 1:50){                      
  mRmR.10[,i] <- mRmR.iter[[i]]$word[1:10]
  colnames(mRmR.10)[i] <- paste0("mrmr", i)
}

library(base)

rank.mRmR.10 <- data.frame("word" = c(mRmR.10[1,1]),
                            "points" = c(0))

rank.mRmR.10$word <- as.character(rank.mRmR.10$word)
rank.mRmR.10$points <- as.integer(rank.mRmR.10$points)

for(i in 1:10) {
  for(j in 1:50){
    if(mRmR.10[i,j] %in% rank.mRmR.10$word){
      rank.mRmR.10[rank.mRmR.10$word==mRmR.10[i,j],2] = (as.integer(rank.mRmR.10[rank.mRmR.10$word==mRmR.10[i,j],2]) + as.integer((11-i)))
    }
    else{
      rank.mRmR.10[nrow(rank.mRmR.10) + 1,] = c(mRmR.10[i,j],as.integer(11-i))
    }
  }
}

rank.mRmR.10$points <- as.integer(rank.mRmR.10$points)
order.rank.mRmR.10 <- rank.mRmR.10[order(rank.mRmR.10$points, decreasing=TRUE)[1:10],]

################################
# ReliefF

library(FSelector)

Relief.iter <- list()
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
  rows.dane <- sample.int(nrow(spam), size = round(nrow(spam)/3), replace = F)
  dane.train <- spam[-rows.dane,]
  dane.test <- spam[rows.dane,]
  
  imp.Relief <- relief(type~., dane.train, neighbours.count = 3, sample.size = 5)
  
  subset <- cutoff.k(imp.Relief, 57)
  imp.Relief <- as.data.frame(subset)
  colnames(imp.Relief)[1] <- "word"
  
  Relief.iter[[i]] <- imp.Relief
  
  relief.best10 <-  Relief.iter[[i]]$word[1:10]
  relief.best10[11] <- "type"
  
  # klasyfikacja
  dane.train <- dane.train[, relief.best10]
  dane.test <- dane.test[, relief.best10]
  
  model.rf <- randomForest(x=dane.train[,-11], y = dane.train[,11],
                           ntree = 500, do.trace = 100)
  rf.result <- predict(model.rf, newdata = dane.test[,-11])
  rf.prob <- predict(model.rf, newdata = dane.test[,-11], type = 'prob')[,2]
  
  # auc
  auc_tmp <- auc_roc(as.numeric(rf.prob)-1, as.numeric(dane.test[,11])-1)
  suma_auc <- suma_auc + auc_tmp
  
  # mcc
  mcc_tmp = mcc(as.numeric(rf.result)-1, as.numeric(dane.test[,11])-1)
  suma_mcc <- suma_mcc + mcc_tmp
  
  # error
  error.rf1 <- sum(dane.test[,11] != rf.result)/length(rf.result)
  suma <- suma + error.rf1
  
  # auc, mcc, error ka¿dej iteracji
  relief.auc.mcc.error[i,1] <- c(auc_tmp)
  relief.auc.mcc.error[i,2] <- c(mcc_tmp)
  relief.auc.mcc.error[i,3] <- c(error.rf1)
  
  relief.rf.prob.list[[i]] <- rf.prob
  relief.rf.result.list[[i]] <- rf.result
  relief.rf.correct.list[[i]] <- dane.test[,11]
}

auc.relief <- suma_auc/i
mcc.relief <- suma_mcc/i
error.relief <- suma/i

###### rankingowanie cech ######


Relief.10=as.data.frame(matrix(nrow=10,ncol=10))

for(i in 1:50){                      
  Relief.10[,i] <- Relief.iter[[i]]$word[1:10]
  colnames(Relief.10)[i] <- paste0("relief", i)
}

library(base)

rank.Relief.10 <- data.frame("word" = c(Relief.10[1,1]),
                            "points" = c(0))

rank.Relief.10$word <- as.character(rank.Relief.10$word)
rank.Relief.10$points <- as.integer(rank.Relief.10$points)

for(i in 1:10) {
  for(j in 1:50){
    if(Relief.10[i,j] %in% rank.Relief.10$word){
      rank.Relief.10[rank.Relief.10$word==Relief.10[i,j],2] = (as.integer(rank.Relief.10[rank.Relief.10$word==Relief.10[i,j],2]) + as.integer((11-i)))
    }
    else{
      rank.Relief.10[nrow(rank.Relief.10) + 1,] = c(Relief.10[i,j],as.integer(11-i))
    }
  }
}

rank.Relief.10$points <- as.integer(rank.Relief.10$points)
order.rank.Relief.10 <- rank.Relief.10[order(rank.Relief.10$points, decreasing=TRUE)[1:10],]

################################
# Random Forest Importance

library(FSelector)

randomForest.iter <- list()
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
  rows.dane <- sample.int(nrow(spam), size = round(nrow(spam)/3), replace = F)
  dane.train <- spam[-rows.dane,]
  dane.test <- spam[rows.dane,]
  
  imp.RF <- random.forest.importance(type~., dane.train, importance.type = 1)
  imp.RandomForest <- cutoff.k(imp.RF, 57)
  imp.RandomForest  <- as.data.frame(imp.RandomForest)
  colnames(imp.RandomForest)[1] <- "word"
  
  randomForest.iter[[i]] <- imp.RandomForest
  
  randomforest.best10 <-  randomForest.iter[[i]]$word[1:10]
  randomforest.best10[11] <- "type"
  
  # klasyfikacja
  dane.train <- dane.train[, randomforest.best10]
  dane.test <- dane.test[, randomforest.best10]
  
  model.rf <- randomForest(x=dane.train[,-11], y = dane.train[,11],
                           ntree = 500, do.trace = 100)
  rf.result <- predict(model.rf, newdata = dane.test[,-11])
  rf.prob <- predict(model.rf, newdata = dane.test[,-11], type = 'prob')[,2]
  
  # auc
  auc_tmp <- auc_roc(as.numeric(rf.prob)-1, as.numeric(dane.test[,11])-1)
  suma_auc <- suma_auc + auc_tmp
  
  # mcc
  mcc_tmp = mcc(as.numeric(rf.result)-1, as.numeric(dane.test[,11])-1)
  suma_mcc <- suma_mcc + mcc_tmp
  
  # error
  error.rf1 <- sum(dane.test[,11] != rf.result)/length(rf.result)
  suma <- suma + error.rf1
  
  # auc, mcc, error ka¿dej iteracji
  rforest.auc.mcc.error[i,1] <- c(auc_tmp)
  rforest.auc.mcc.error[i,2] <- c(mcc_tmp)
  rforest.auc.mcc.error[i,3] <- c(error.rf1)
  
  rforest.rf.prob.list[[i]] <- rf.prob
  rforest.rf.result.list[[i]] <- rf.result
  rforest.rf.correct.list[[i]] <- dane.test[,11]
}

auc.randomforest <- suma_auc/i
mcc.randomforest <- suma_mcc/i
error.randomforest <- suma/i

###### rankingowanie cech ######


randomForest.10=as.data.frame(matrix(nrow=10,ncol=10))

for(i in 1:50){                      
  randomForest.10[,i] <- randomForest.iter[[i]]$word[1:10]
  colnames(randomForest.10)[i] <- paste0("rforest", i)
}

library(base)

rank.randomForest.10 <- data.frame("word" = c(randomForest.10[1,1]),
                                    "points" = c(0))

rank.randomForest.10$word <- as.character(rank.randomForest.10$word)
rank.randomForest.10$points <- as.integer(rank.randomForest.10$points)

for(i in 1:10) {
  for(j in 1:50){
    if(randomForest.10[i,j] %in% rank.randomForest.10$word){
      rank.randomForest.10[rank.randomForest.10$word==randomForest.10[i,j],2] = (as.integer(rank.randomForest.10[rank.randomForest.10$word==randomForest.10[i,j],2]) + as.integer((11-i)))
    }
    else{
      rank.randomForest.10[nrow(rank.randomForest.10) + 1,] = c(randomForest.10[i,j],as.integer(11-i))
    }
  }
}

rank.randomForest.10$points <- as.integer(rank.randomForest.10$points)
order.rank.randomForest.10 <- rank.randomForest.10[order(rank.randomForest.10$points, decreasing=TRUE)[1:10],]

#############################################
# porównanie klasyfikacji dla wszystkich metod

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
                                      "error" = error.relief))


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
        col = c("red", "green", "blue", "magenta", "gray", "dark green"), cex=2)


#############################################
# porównanie rankingów selekcji cech wszystkich metod

all.methods.rank.10 <- as.data.frame(cbind("Ttest" = order.rank.Ttest.10$word, 
                               "MDFS" = order.rank.MDFS.10$word, 
                               "Boruta" = order.rank.Boruta.10$word, 
                               "mRmR" = order.rank.mRmR.10$word, 
                               "ReliefF" = order.rank.Relief.10$word, 
                               "randomForest" = order.rank.randomForest.10$word))

#################################
# PUNKTACJA CECH Z WSZYSTKICH METOD
# RANKING

library(base)

rank <- data.frame("word" = c(all.methods.rank.10[1,1]),
                   "points" = c(0))

rank$word <- as.character(rank$word)
rank$points <- as.integer(rank$points)

for(i in 1:10) {
  for(j in 1:6){
    if(all.methods.rank.10[i,j] %in% rank$word){
      rank[rank$word==all.methods.rank.10[i,j],2] = (as.integer(rank[rank$word==all.methods.rank.10[i,j],2]) + as.integer((11-i)))
    }
    else{
      rank[nrow(rank) + 1,] = c(all.methods.rank.10[i,j],as.integer(11-i))
    }
  }
}

rank$points <- as.integer(rank$points)
order.rank <- rank[order(rank$points, decreasing=TRUE),]
