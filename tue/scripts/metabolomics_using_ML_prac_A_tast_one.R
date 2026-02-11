# import libraries
library("class")
library("gmodels")
library("caret")
library("rpart")
library("rpart.plot")
library("mixOmics")
library("ggplot2")

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
############ 1-K Nearest Neightbours ###############
####################################################

# remove class labels from training and test sets for KNN
train_set_knn <- train_set[, -ncol(train_set)]
test_set_knn <- test_set[, -ncol(test_set)]

# perform KNN with k = 3
model_k3 <- knn(
    train = train_set_knn, test = test_set_knn,
    cl = train_class_labels, k = 3
)

# take a look at the summary of the model
summary(model_k3)


####################################################
################# Model evaluation #################
####################################################

# create cross table of predicted vs actual class labels
cross_table <- CrossTable(test_class_labels, model_k3,
    prop.chisq = FALSE, prop.t = FALSE, prop.c = FALSE,
    prop.r = FALSE
)

# create confusion matrix
confusion_matrix <- confusionMatrix(model_k3,
    test_class_labels,
    positive = "3"
)

# print confusion matrix
print(confusion_matrix)

####################################################
################# Model Improvment #################
####################################################

k_results <- function(n, train, test, train_class_labels, test_class_labels) {
    model_k <- knn(
        train = train, test = test,
        cl = train_class_labels, k = n
    )

    confusion_matrix <- confusionMatrix(model_k,
        test_class_labels,
        positive = "3"
    )

    matrix_accuracy <- confusion_matrix$overall["Accuracy"]

    return(matrix_accuracy)
}

# test different values of k
k_values <- seq(1, 20, by = 1)
accuracy_results <- sapply(k_values, k_results,
    train = train_set_knn, test = test_set_knn,
    train_class_labels = train_class_labels,
    test_class_labels = test_class_labels
)

# create a data frame of k values and their corresponding accuracies
accuracy_df <- data.frame(k = k_values, accuracy = accuracy_results)

# plot the accuracy results
ggplot(accuracy_df, aes(x = k, y = accuracy)) +
    geom_line() +
    geom_point() +
    labs(x = "Number of Neighbors (k)", y = "Accuracy", title = "KNN Accuracy for Different k Values") +
    theme_minimal()

####################################################
##################### Sacaling #####################
####################################################

pre_process_values <- preProcess(train_set_knn, method = c("center", "scale"))
train_transformed <- predict(pre_process_values, train_set_knn)
test_transformed <- predict(pre_process_values, test_set_knn)

# perform KNN with k = 3 on scaled data
model_k3_scaled <- knn(
    train = train_transformed, test = test_transformed,
    cl = train_class_labels, k = 3
)

# use k results function to test different k values on scaled data
accuracy_results_scaled <- sapply(k_values, k_results,
    train = train_transformed, test = test_transformed,
    train_class_labels = train_class_labels,
    test_class_labels = test_class_labels
)

# create a data frame of k values and their corresponding accuracies for scaled data
accuracy_df_scaled <- data.frame(k = k_values, accuracy = accuracy_results_scaled)

# plot the accuracy results for scaled data
ggplot(accuracy_df_scaled, aes(x = k, y = accuracy)) +
    geom_line() +
    geom_point() +
    labs(x = "Number of Neighbors (k)", y = "Accuracy", title = "KNN Accuracy for Different k Values (Scaled Data)") +
    theme_minimal()
