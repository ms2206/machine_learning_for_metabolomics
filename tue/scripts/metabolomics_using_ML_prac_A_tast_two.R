# Model training — Decision Tree
library("rpart")
library("rpart.plot")


####################################################
################# Data Loading #####################
####################################################

# read data
data.file.path <- "/Users/mspriggs/Library/CloudStorage/OneDrive-Illumina,Inc./Documents/Applied_Bioinformatics/modules/machine_learning_for_metabolomics/tue/pracs/"

enose <- read.table(paste0(data.file.path, "EnoseAllSamples.csv"),
    header = TRUE, sep = ",", row.names = 1
)

sensory <- read.table(paste0(data.file.path, "SensoryAllSamples.csv"),
    header = TRUE, sep = ",", row.names = 1
)

# merge the two datasets by row names (sample names)
enose_sensory <- merge(enose, sensory, by = "row.names")
# set row names to the first column (which contains the sample names)
rownames(enose_sensory) <- enose_sensory[, 1]

# remove the first column (which contains the sample names)
enose_sensory <- as.data.frame(enose_sensory[, -1])

####################################################
############### Data Exploration ###################
####################################################

# use mixOmics to perform PCA on the combined dataset
pca_combined <- pca(enose_sensory, ncomp = 3, scale = TRUE)

# plot the PCA results
quartz()
plotIndiv(pca_combined,
    ind.names = enose_sensory$sensory,
    group = enose_sensory$sensory,
    style = "ggplot2", legend = TRUE,
    title = "PCA of combined eNose and Sensory data"
)

####################################################
############### Data Preparation ###################
####################################################

# extract class labels from sensory column
enose_sensory$sensory <- factor(enose_sensory$sensory,
    levels = c("1", "2", "3")
)

# check the proportions of each class to see if the dataset is balanced
round(prop.table(table(enose_sensory$sensory)), digits = 1)

# create training and test set
set.seed(8)
train_index <- createDataPartition(enose_sensory$sensory,
    p = 0.7, list = FALSE, times = 1
)

# create training and test sets
train_set <- enose_sensory[train_index, ]
train_class_labels <- train_set[, ncol(train_set)]

test_set <- enose_sensory[-train_index, ]
test_class_labels <- test_set[, ncol(test_set)]

####################################################
################ Desision Trees ####################
####################################################

# train a decision tree model using rpart
model_tree <- rpart(sensory ~ ., data = train_set)

# plot the decision tree
quartz()
rpart.plot(model_tree, main = "Decision Tree for Sensory Classification")

####################################################
################ Model Evaluation ##################
####################################################

# make predictions on the test set
predicted <- predict(model_tree, newdata = test_set, type = "class")

# create a confusion matrix
conf_matrix <- confusionMatrix(predicted, test_class_labels)

# print the confusion matrix and performance metrics
print(conf_matrix)

####################################################
################ Model Improvment ##################
####################################################

# perform cross-validation to find the optimal complexity parameter (cp) for pruning the tree
pruned_tree <- prune(model_tree, cp = 0.04)

# create a function that can calculate the prediction accuracy for the different values of cp and store the results in a vector
cp_values <- seq(0.01, 0.1, by = 0.01)

prediction_accuracy <- function(cp) {
    pruned_tree <- prune(model_tree, cp = cp)
    predicted <- predict(pruned_tree, newdata = test_set, type = "class")
    conf_matrix <- confusionMatrix(predicted, test_class_labels)
    return(conf_matrix$overall["Accuracy"])
}

accuracy_results <- sapply(cp_values, prediction_accuracy)

# plot the accuracy results against the cp values
quartz()
plot(cp_values, accuracy_results,
    type = "b",
    xlab = "Complexity Parameter (cp)",
    ylab = "Prediction Accuracy",
    main = "Prediction Accuracy vs Complexity Parameter"
)
