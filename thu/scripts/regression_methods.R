# load libraies
library(tidyverse)
library(ggplot2)
library(caret)
library(MASS)


####################################################
##################### Load Data  ###################
####################################################

# load data
data_folder <- "/Users/mspriggs/Library/CloudStorage/OneDrive-Illumina,Inc./Documents/Applied_Bioinformatics/modules/machine_learning_for_metabolomics/thu/prac/Regression_Practical"

# read in data
enose <- read.table(file.path(data_folder, "EnoseAllSamples.csv"), sep = ",", header = TRUE)
bacterial_counts <- read.table(file.path(data_folder, "MicroCounts.csv"), sep = ",", header = TRUE)

# merge data on Samples_0
data <- merge(enose, bacterial_counts, by = "row.names")
rownames(data) <- data[, 1]
data <- as.data.frame(data[, -1])


# conver CFC to log scale
data$log_CFC <- log10(data$CFC)

####################################################
########### Simple Descriptive Plots ###############
####################################################

# histogram
ggplot(data, aes(x = log_CFC)) +
    geom_histogram(binwidth = 0.2, fill = "gray", color = "black") +
    labs(title = "Log Histogram CFC Counts", x = "Bacterial Counts", y = "Frequency") +
    theme_minimal()

####################################################
################# Data Splitting ###################
####################################################

# Data splitting with createDataPartition
set.seed(123)

train_index <- createDataPartition(data$CFC, p = 0.7, list = FALSE)
train_data <- data[train_index, ]
test_data <- data[-train_index, ]

####################################################
########### Multiple Linear Regression #############
####################################################

# fit linear regression model with non-log transformed CFC
lm(CFC ~ ., data = train_data) %>% summary()

# fit linear regression model with log transformed CFC
lm(log_CFC ~ ., data = train_data) %>% summary()

# use varImp to get variable importance
lm_model <- lm(log_CFC ~ DF1 + DF1.1 + DF3 + DF4 + DF5, data = train_data)
lm_model %>% summary()
varImp(lm_model)


####################################################
############### Assessing Fit Quality ##############
####################################################

# You can use the predict() function to produce predicted bacterial count values for your test data.
predicted_counts <- predict(lm_model, newdata = test_data)

# You can then visualise the result of fitting by creating a scatter plot of predicted vs actual values:
ggplot(test_data, aes(x = log_CFC, y = predicted_counts)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = "Predicted vs Actual Log CFC Counts", x = "Actual Log CFC", y = "Predicted Log CFC") +
    theme_minimal()

# calculate RMSE
rmse <- RMSE(predicted_counts, test_data$log_CFC)
print(paste("RMSE:", rmse))

####################################################
#################### Model Tuning ##################
####################################################

# tune with lmAIC
model_fit <- train(CFC ~ . - Samples_0 - log_CFC, data = train_data, method = "lmStepAIC")
summary(model_fit)

# tune with knn
model_fit <- train(CFC ~ . - Samples_0 - log_CFC,
    method = "knn", data = train_data,
    tuneGrid = expand.grid(k = 1:20)
)
print(model_fit)
