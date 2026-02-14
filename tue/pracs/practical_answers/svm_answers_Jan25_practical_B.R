##### ML Practical B Answers - Classification with ML #########

##### SVM #######


# Clear workspace
rm(list=ls())

# Close any open graphics devices
graphics.off()


library(caret)
library(kernlab)
library(gmodels)

enose   <- read.table("EnoseAllSamples.csv", sep=",", header=TRUE, row.names=1)
sensory <- read.table("SensoryAllSamples.csv ", header = TRUE, sep =",",row.names=1)
# Match rows from enose to rows from sensory
merged <- merge(enose, sensory, by="row.names")
rownames(merged) = merged[,1]
#remove row names column
AllData<-as.data.frame(merged[,-1])

# Define classes as factor
AllData$sensory<- factor(AllData$sensory, levels = c("1", "2", "3"))
#Check proportion of each class
round(prop.table(table(AllData$sensory)) * 100, digits = 1)

#Partition data into training and test set
##Randomly select 70% of sample row numbers
set.seed(18)
trainIndex <- createDataPartition(AllData$sensory, p = .7, 
                                  list = FALSE, 
                                  times = 1)

#Create training set and test set
trainSet <- AllData[trainIndex,]
testSet <- AllData[-trainIndex,]

#separate training and test classes
trainCl <- trainSet[,ncol(trainSet)]
testCl <- testSet[,ncol(testSet)]


#svm model
model.svm <- train(sensory ~ ., data = trainSet, method="svmLinear2")
model.svm

## Predict test set
predicted <- predict(model.svm, testSet)

#Produce cross table and confusion matrix
cross.table <- CrossTable(testCl, predicted, prop.chisq=FALSE, prop.t=FALSE, prop.c=FALSE, prop.r=FALSE) ##Note here we give the predicted values as argument unlike knn where we give the model
confusion.matrix <- confusionMatrix(predicted, testCl, positive="3")
confusion.matrix
cat('SVM accuracy: ', confusion.matrix$overall[1])


#write a loop to calculate prediction accuracies for 100 data partitions
  accuracies<-c()
  trainIndex1 <- createDataPartition(AllData$sensory, p = .7, 
                                      list = FALSE, 
                                      times = 100)
    for (i in 1:100){
    trainSet <- AllData[trainIndex1[,i],]
    testSet <- AllData[-trainIndex1[,i],]
    testClass <-testSet[, ncol(testSet)]
    model_svm <- train(sensory ~ ., data=trainSet, method="svmLinear2")
    kernel.predicted <- predict(model_svm, testSet)
    kernel.confusion.matrix <- confusionMatrix(kernel.predicted, testClass, positive="3")
    current.accuracy <- kernel.confusion.matrix$overall[1]
    accuracies<-c(accuracies, current.accuracy)
}

iteration <- as.array(c(1:100))
plot(iteration, accuracies, type='l', ylim=c(0.2,1))


##Plot a histogram of the accuracy measurement, which will allow you to assess the accuracy distribution
hist(accuracies, xlab="model accuracy", main="SVM classifier prediction accuracy")

###Model improvement 
###Try different kernels

svm.methods <- c('lssvmRadial', 'lssvmPoly', 'svmRadial', 'svmLinear3')

#function to calculate accuracy for 100 iterations of the dataset for different kernel types

accuracy.m <-function(data, svm.method.name){
  accuracies.m <-c()
  for(j in 1:100){
    classes<-data[,ncol(data)]
    trainIndex <- createDataPartition(as.factor(classes), p = .7, 
                                      list = FALSE, 
                                      times = 1)
    trainIndex<-as.matrix(trainIndex)
    trainSet <- data[trainIndex,]
    testSet <- data[-trainIndex,]
    testC <-testSet[, ncol(testSet)]
    model_svm <- caret::train(sensory ~ ., data=trainSet, kernel=svm.method.name)
    kernel.predicted <- predict(model_svm, testSet)
    kernel.confusion.matrix <- confusionMatrix(kernel.predicted, testC, positive="3")
    current.accuracy <- kernel.confusion.matrix$overall[1]
    accuracies.m<-c(accuracies.m, current.accuracy)
      }
  return(accuracies.m)
  
}


#Test function with one of the kernels
a<-accuracy.m(AllData, "svmRadial")

plot (1:100, a, type="l")

#run a for loop to generate a vector with the mean accuracies for each kernel
mean.accuracies<-c()
for(i in svm.methods){
  b<-accuracy.m(AllData, svm.method.name =i)
  m.mean.accuracy<-mean(b)
  mean.accuracies<-c(mean.accuracies, m.mean.accuracy)
}

barplot(mean.accuracies, names.arg=svm.methods, ylab="mean accuracy", ylim=c(0.0, 0.9))

## All SVM kernel have a mean accuracy around 0.76-0.77. 

## You can try to run the analysis using different svm types available in the caret library:
## the https://topepo.github.io/caret/train-models-by-tag.html#Support_Vector_Machines
## Alternatively you can use a library dedicated to svm such as kernlab. 






