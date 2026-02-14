##### ML Practical A - Classification with ML #########

##### k-nn and Decision Trees #######


rm(list=ls())

# Close any open graphics devices
graphics.off()

library(class)
library(gmodels)
library(caret)
library(rpart)
library(mixOmics)
library(rpart.plot)

#Load the data
enose   <- read.table("EnoseAllSamples.csv", sep=",", header=TRUE, row.names=1)
sensory <- read.table("SensoryAllSamples.csv ", header = TRUE, sep =",",row.names=1)
# Match rows from enose to rows from sensory
merged <- merge(enose, sensory, by="row.names")
#Define rownames
rownames(merged) = merged[,1]
#remove row names column
AllData<-as.data.frame(merged[,-1])

# PCA
pca.enose <- pca(AllData[,-ncol(AllData)], ncomp=4, scale=TRUE)
x11();plotIndiv(pca.enose, ind.names=as.character(AllData[,ncol(AllData)]), group=as.factor(AllData[,ncol(AllData)]), style="lattice")
# no differentiation

# Define classes as factor
AllData$sensory<- factor(AllData$sensory, levels = c("1", "2", "3"))
#Check proportion of each class
round(prop.table(table(AllData$sensory)) * 100, digits = 1)

#Partition data into training and test set
##Randomly select 70% of sample row numbers
set.seed(8)
trainIndex <- createDataPartition(AllData$sensory, p = .7, 
                                  list = FALSE, 
                                  times = 1)

#Create training set and test set
trainSet <- AllData[trainIndex,]
testSet <- AllData[-trainIndex,]

#separate training and test classes
trainCl <- trainSet[,ncol(trainSet)]
testCl <- testSet[,ncol(testSet)]


##Run knn model with k=3
trainSet.knn <-trainSet[, -ncol(trainSet)]
testSet.knn <- testSet[, -ncol(testSet)]
model.k3 <- knn(trainSet.knn, testSet.knn, trainCl, k=3)

summary(model.k3) #Outputs how many samples from the training set were attributed to each class. No model is produced

cross.table <- CrossTable(testCl, model.k3, prop.chisq=FALSE, prop.t=FALSE, prop.c=FALSE, prop.r=FALSE)
confusion.matrix <- confusionMatrix(model.k3, testCl, positive="3")
confusion.matrix

##Try with different k(1:20)
array.of.ks <- c(1:20)
k.results<-function(n, trS, tstS, trCl, tstCl){
  k.accuracy<-c()
  for (k in 1:n) {
    model.k<-knn(trS, tstS, trCl, k)
    confusion.matrix <- confusionMatrix(model.k, tstCl, positive="3")
    k.accuracy[k]<-confusion.matrix$overall[1]
  }
  return(k.accuracy)
}

#Apply function
k.accuracy<-k.results(20, trainSet.knn, testSet.knn, trainCl, testCl)
plot(array.of.ks, k.accuracy, type='l')

###Scale data and try again with k=3
preProcValues <- preProcess(trainSet.knn, method = c("center", "scale"))

trainTransformed <- predict(preProcValues, trainSet.knn)
testTransformed <- predict(preProcValues, testSet.knn)

#Create model with scaled data
sc.model.k3<-knn(trainTransformed, testTransformed, trainCl, 3)

summary(sc.model.k3)

cross.table <- CrossTable(testCl, sc.model.k3, prop.chisq=FALSE, prop.t=FALSE, prop.c=FALSE, prop.r=FALSE)
confusion.matrix <- confusionMatrix(sc.model.k3, testCl, positive="3")
confusion.matrix


##Try with different k(1:20)
array.of.ks <- c(1:20)
k.accuracy.sc<-k.results(20, trainTransformed, testTransformed, trainCl, testCl)
plot(array.of.ks, k.accuracy.sc, type='l')


##Decision trees
model.tree <- rpart(sensory ~., data=trainSet)
#visualise the tree
x11()
rpart.plot(model.tree, box.palette = "BuBn", type = 5)
summary(model.tree)

#Confusion matrix
predicted <- predict(model.tree, testSet, type="class")
confusion.matrix.tree <- confusionMatrix(predicted, testCl, positive="3")
confusion.matrix.tree
cat('Tree accuracy: ', confusion.matrix.tree$overall[1])

#Prune tree 
pruned_tree <-prune(model.tree, cp = 0.02)

x11()
rpart.plot(pruned_tree, box.palette = "BuBn", type = 5)
summary(pruned_tree)

#performance of pruned tree
predicted.pr <- predict(pruned_tree, testSet, type="class")
confusion.matrix.pruned <- confusionMatrix(predicted.pr, testCl, positive="3")
confusion.matrix.pruned
cat('Pruned Tree accuracy: ', confusion.matrix.pruned$overall[1])

#create function to store array of accuracy for different best values
accuracy.cp<-function(model.tr, tstS, tstCl) {
  cp.accuracy<-c()
  array.of.cp <- seq(from = 0.01, to = 0.1, by = 0.01)
  for (cp in array.of.cp) {
    pruned.tree <- prune(model.tr, cp = cp)
    predicted.pr <- predict(pruned.tree, tstS, type="class")
    confusion.matrix <- confusionMatrix(predicted.pr, tstCl, positive="3")
    current.acc<-confusion.matrix$overall[1]
    cp.accuracy<-c(cp.accuracy, current.acc)
  }
  return(cp.accuracy)
}

accur<-accuracy.cp(model.tree, testSet, testCl) 
array.of.cp <- seq(from = 0.01, to = 0.1, by = 0.01)

x11()
plot(array.of.cp, accur, type='l')

#The accuracy is highest when all the features are included (cp=0.01). Removing features reduces the accuracy of this model. 
#However, for other datasets it may improve the accuracy.

