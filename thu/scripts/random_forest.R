# load libraries
library(mlr3verse)
library(ranger)
library(tidyverse)
library(mlbench)
library(mlr3viz)
library(mlr3learners)
library(paradox)
library(bbotk)

####################################################
################ Data Loading ######################
####################################################

# load the zoo data
data("Zoo", package = "mlbench")

# convert to data.table
zoo_dt <- as.data.table(Zoo)

zoo_dt <- mutate_all(Zoo, as.factor)

round(prop.table(table(zoo_dt$type)) * 100, digits = 1)

a <- which(zoo_dt$legs == "8")
b <- which(zoo_dt$legs == "5")

zoo_dt_balanced <- zoo_dt[-c(a, b), ]
round(prop.table(table(zoo_dt_balanced$type)) * 100, digits = 1)

####################################################
################ Model Building ####################
####################################################
task_zoo <- as_task_classif(type ~ ., data = zoo_dt_balanced)

autoplot(task_zoo)

learner <- lrn("classif.ranger")

####################################################
################ Model Training ####################
####################################################
set.seed(4)
split <- partition(task_zoo, ratio = 0.67)

learner$train(task_zoo, split$train)

learner$model

####################################################
################ Model Testing #####################
####################################################

prediction <- learner$predict(task_zoo, split$test)

# classification accuracy
measure <- msr("classif.acc")
prediction$score(measure)
prediction$confusion

autoplot(prediction)

####################################################
################ Model Tuning ######################
####################################################
learner$param_set$ids()

learner1 <- lrn("classif.ranger",
    num.trees = 500, mtry = to_tune(2, 12), min.node.size = 2, max.depth = 20
)

resampling <- rsmp("cv", folds = 3)

measure <- msr("classif.acc")

terminator <- trm("evals", n_evals = 20)

# deinfe a new task called task_train
train_set <- zoo_dt_balanced[split$train, ]
task_train <- as_task_classif(type ~ ., data = train_set)

# now put everything together in a TuningInstance instance
instance <- ti(
    task = task_train,
    learner = learner1,
    resampling = resampling,
    measure = measure,
    terminator = terminator
)

instance

tuner <- tnr("grid_search", resolution = 5, batch_size = 4)
tuner$optimize(instance)

# build the final model
learner1$param_set$values <- instance$result_learner_param_vals
learner1$train(task_train)
learner1$model

####################################################
################ Model Testing #####################
####################################################
prediction <- learner1$predict(task_zoo, split$test)
# classification accuracy
measure <- msr("classif.acc")
prediction$score(measure)

prediction$confusion
autoplot(prediction)
