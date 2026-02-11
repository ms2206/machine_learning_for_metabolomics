# Practical B: Metabolomics using Machine Learning
# load necessary libraries
library("caret")
library("gmodels")
library("kernlab")
library("LiblineaR")

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
############### Data Preperation ###################
####################################################

# extract class labels from sensory column
enose_sensory$sensory <- factor(enose_sensory$sensory,
    levels = c("1", "2", "3")
)


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
test_class_labels <- factor(test_class_labels, levels = c("1", "2", "3"))

####################################################
################ Model Training ####################
####################################################

# train a support vector machine (SVM) model using the training set
svm_model <- train(
    sensory ~ .,
    data = train_set, method = "svmLinear2"
)

svm_model

####################################################
################# Model evaluation #################
####################################################

# make predictions on the test set using the trained SVM model
svm_predictions <- predict(svm_model, newdata = test_set)
svm_predictions <- factor(svm_predictions, levels = c("1", "2", "3"))


# evaluate the performance of the SVM model using a CrossTable
CrossTable(test_class_labels, svm_predictions,
    prop.chisq = FALSE, prop.t = FALSE, prop.r = FALSE,
    dnn = c("Actual", "Predicted")
)

# confusion matrix
confusionMatrix(svm_predictions, test_class_labels, positive = "3")

# generate 100 random data points for each class
# Generate 100 random data points for each class
set.seed(123)

# Set vector to store the accuracies of the models
accuracies <- c()

for (i in 1:100) {
    # Create new training and test indices for each iteration
    train_index <- createDataPartition(enose_sensory$sensory, p = 0.7, list = FALSE, times = 1)

    # Create training and test sets
    train_set <- enose_sensory[train_index, ]
    test_set <- enose_sensory[-train_index, ]
    test_class_labels <- factor(test_set[, ncol(test_set)], levels = c("1", "2", "3"))

    # Train the SVM model
    model_svm <- train(
        sensory ~ .,
        data = train_set, method = "svmLinear2"
    )

    # Make predictions on the test set
    kernal_predictions <- predict(model_svm, newdata = test_set)
    kernal_predictions <- factor(kernal_predictions, levels = c("1", "2", "3"))

    # Calculate accuracy of the model
    confusion_matrix <- confusionMatrix(kernal_predictions, test_class_labels, positive = "3")
    accuracy <- confusion_matrix$overall["Accuracy"]
    accuracies <- c(accuracies, accuracy)
}

# Print the accuracies
print(accuracies)

# plot the distribution of accuracies
hist(accuracies, main = "Distribution of SVM Accuracies", xlab = "Accuracy", ylab = "Frequency", col = "lightblue")

####################################################
################# Model Improvment #################
####################################################

svm_methods <- c("lssvmLinear", "lssvmPoly", "svmRadial", "svmLinear3")

# create a function to train and evaluate SVM models with different kernels which applies 100 models for each kernel and stores the accuracies in a list
accuracy_models <- function(data, svc_method_name) {
    accuracies <- c()
    for (i in 1:100) {
        classes <- data[, ncol(data)]
        train_index <- createDataPartition(classes, p = 0.7, list = FALSE, times = 1)
        train_index <- as.matrix(train_index)
        train_set <- data[train_index, ]
        test_set <- data[-train_index, ]
        test_class_labels <- factor(test_set[, ncol(test_set)], levels = c("1", "2", "3"))
        model_svm <- train(
            sensory ~ .,
            data = train_set, kernel = svc_method_name
        )
        kernal_predictions <- predict(model_svm, newdata = test_set)
        kernal_confusion_matrix <- confusionMatrix(kernal_predictions, test_class_labels, positive = "3")
        accuracy <- kernal_confusion_matrix$overall["Accuracy"]
        accuracies <- c(accuracies, accuracy)
    }
    return(accuracies)
}

# create a list to store the accuracies of the different models
model_accuracies <- list()
# loop through the different SVM methods and store the accuracies in the list
for (method in svm_methods) {
    model_accuracies[[method]] <- accuracy_models(enose_sensory, method)
}
# print the accuracies of the different models
print(model_accuracies)

# find the mean accuracy of each model
mean_accuracies <- sapply(model_accuracies, mean)
print(mean_accuracies)
