###### Practical Ensemble #########

#####Random Forests##############

## Maria Anastasiadi

# Clear workspace
rm(list=ls())

# Close any open graphics devices
graphics.off()


#install.packages("mlr3")
#install.packages("mlr3verse")
#install.packages("mlbench")
#install.packages("parallelMap")
#install.packages("xgboost")
#install.packages("ranger")
#install.packages("devtools")

library("mlr3verse")
library(mlr3viz)
library(mlr3learners)
library(ranger)
library(mlbench)
library(tidyverse)
library(paradox)
library(bbotk)

data(Zoo)

#zooD <- as.data.frame(Zoo)
head(Zoo)
zooD <- mutate_all(Zoo, as.factor)
head(zooD)


#Check proportion of each class
round(prop.table(table(zooD$type)) * 100, digits = 1)

##remove samples with number of legs 5,8 which are underrepresented 
a<-which(zooD$legs=="8")
b<-which(zooD$legs=="5")
zoo1<-zooD[-c(a,b),]

# create learning task
task_zoo = as_task_classif(type ~ ., data = zoo1)

task_zoo
task_zoo$data()
# default plot: class frequencies
autoplot(task_zoo)


#Check available algorithms for ML 
mlr_learners

#Select classif.ranger

# load learner and set hyperparameter
learner = lrn("classif.ranger")
learner$param_set$ids()
learner$param_set

# train/test split
set.seed(4)
split = mlr3::partition(task_zoo, ratio = 0.67)

# train the model
learner$train(task_zoo, split$train)
learner$model

# predict data
prediction = learner$predict(task_zoo, split$test)

# calculate performance
#classification accuracy
measure = msr("classif.acc")
prediction$score(measure)

#Check all performance measures available
mlr_measures

#confusion matrix
prediction$confusion


library("mlr3viz")
autoplot(prediction)

### Hyperparameter optimisation

## check hyperparameters to tune
learner$param_set$ids()

# set tuning parameters
learner1 = lrn("classif.ranger",
              num.trees = to_tune(200, 500),
              mtry  = to_tune(2, 12),
              min.node.size = 2,
              max.depth = 20
)

resampling = rsmp("cv", folds = 3)

measure = msr("classif.acc")

terminator = trm("evals", n_evals = 20)
terminator

#Now we put everything together into a TuningInstanceSingleCrit with the ti() function.
train_Set<-zoo1[split$train,]
task_train<-as_task_classif(type ~ ., data = train_Set)
instance = ti(
  task = task_train,
  learner = learner1,
  resampling = resampling,
  measures = measure,
  terminator = terminator
)

instance

#Set up the tuner
tuner = tnr("grid_search", resolution = 5, batch_size = 4)
tuner

#Trigerring tuner 
tuner$optimize(instance)

head(as.data.table(instance$archive))

#Final model
learner1$param_set$values = instance$result_learner_param_vals

learner1$train(task_train)
learner1$model

##Test the new model
# predict data
#test_Set <-zoo1[split$test,]
prediction = learner1$predict(task_zoo, split$test)

# calculate performance
#classification accuracy
measure = msr("classif.acc")
prediction$score(measure)

#confusion matrix
prediction$confusion

autoplot(prediction)

##NOTE: Tuning the RF model did not achieve an improvement in the performance, 
# but the performance of the original model was already high

###################################################################

########## Pima Indian Diabetes ##################################

data("PimaIndiansDiabetes2", package = "mlbench")
PimaIndiansDiabetes2 <- as.data.frame(PimaIndiansDiabetes2)
# Omit NAs
PimaIndiansDiabetes2 <- na.omit(PimaIndiansDiabetes2)
# create learning task
task_d = as_task_classif(diabetes ~ ., data = PimaIndiansDiabetes2)

task_d

# default plot: class frequencies
autoplot(task_d)


#Check available algorithms for ML 
mlr_learners

#Select classif.randomForest

# load learner and set hyperparameter
learner = lrn("classif.ranger")
learner$param_set$ids()
learner$param_set

# train/test split
set.seed(4)
split = partition(task_d, ratio = 0.67)

# train the model
learner$train(task_d, split$train)
learner$model

# predict data
prediction = learner$predict(task_d, split$test)

# calculate performance
#classification accuracy
measure = msr("classif.acc")
prediction$score(measure)

#Check all performance measures available
mlr_measures

#confusion matrix
prediction$confusion


learner1 = lrn("classif.ranger",
               num.trees = to_tune(200, 500),
               mtry  = to_tune(2, 8),
               min.node.size = 2,
               max.depth = 20)

resampling = rsmp("cv", folds = 3)

measure = msr("classif.acc")

terminator = trm("evals", n_evals = 50)
terminator

#Now we put everything together into a TuningInstanceSingleCrit with the ti() function.
train_Set<-PimaIndiansDiabetes2[split$train,]
task_train<-as_task_classif(diabetes ~ ., data = train_Set)
instance = ti(
  task = task_train,
  learner = learner1,
  resampling = resampling,
  measures = measure,
  terminator = terminator
)

instance

#Set up the tuner
tuner = tnr("grid_search", resolution = 5, batch_size = 4)
tuner

#Trigerring tuner 
tuner$optimize(instance)

head(as.data.table(instance$archive))

#Final model
learner1$param_set$values = instance$result_learner_param_vals

learner1$train(task_d)
learner1$model

##Test the new model
# predict data
#test_Set <-zoo1[split$test,]
prediction = learner1$predict(task_d, split$test)

# calculate performance
#classification accuracy
measure = msr("classif.acc")
prediction$score(measure)

measures = msrs(c("classif.tpr", "classif.tnr"))
prediction$score(measures)

#confusion matrix
prediction$confusion

autoplot(prediction)

## NOTE: In the case of the PIMA Indian model you should see a substantial improvement from ~0.7 to 1. 

