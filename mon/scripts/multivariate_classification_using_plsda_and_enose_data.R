# import libraries
library("rgl")
library("matlab")
library("R.matlab")
library("ptw")
library("dplyr")
library("mixOmics")

# read data
data.file.path <- "/Users/mspriggs/Library/CloudStorage/OneDrive-Illumina,Inc./Documents/Applied_Bioinformatics/modules/machine_learning_for_metabolomics/mon/prac/Practicals-Monday/data2.csv" # nolint
data <- read.table(data.file.path,
    header = TRUE, sep = ",", row.names = 1
)

# extract column names
sensors <- colnames(data)

# extract row names
samples <- rownames(data)

# extract class labels from Class column
class_labels <- data$Class

# resp contains data without class labels
resp <- data[, -1]
resp <- as.matrix(resp)

####################################################
####### exploratory data analysis with PCA #########
####################################################

# pca
pca_enose <- pca(resp, ncomp = 6, scale = TRUE)

# define class vector
class_labels <- as.factor(class_labels)


plotIndiv(pca_enose,
    ind.names = samples,
    group = class_labels,
    style = "ggplot2", legend = TRUE, title = "PCA of eNose data"
)


# extract the varience from the pca object
var_1 <- pca_enose$prop_expl_var$X[1]
var_2 <- pca_enose$prop_expl_var$X[2]

# convert to percentage
var_1 <- round(var_1 * 100, 2)
var_2 <- round(var_2 * 100, 2)


biplot(pca_enose,
    xlab = paste("PC1 ", var_1, "%"),
    ylab = paste("PC2 ", var_2, "%"),
    group = class_labels,
    col.per.group = c("red", "blue", "green", "orange", "purple"),
    title = "PCA of eNose data"
)

####################################################
#################### PLS-DA ########################
####################################################

# create training and test sets
set.seed(123)
train_index <- resp[1:20, ] # use the first 20 samples for training
# create class vector for training set
train_class_labels <- class_labels[1:20]
train_class_labels <- as.factor(train_class_labels)

# create test set
test_index <- resp[21:25, ] # use the last 5 samples for testing
# create class vector for test set
test_class_labels <- class_labels[21:25]
test_class_labels <- as.factor(test_class_labels)

# create plsda model
plsda_train_1 <- plsda(train_index, train_class_labels, ncomp = 10)

plotIndiv(plsda_train_1,
    ind.names = TRUE,
    comp = c(1, 2),
    ellipse = TRUE,
    style = "lattice",
    cex = c(rep(1, 5))
)

# variable importance in projection (VIP) scores
train_vip <- vip(plsda_train_1)

# visualise VIP scores
ggplot(
    data = train_vip,
    aes(x = rownames(train_vip), y = comp1, fill = rownames(train_vip))
) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Variable", y = "VIP Score", title = "VIP Scores for PLS-DA Model") +
    theme(legend.position = "none")

# create a correlation plot for the variables
plotVar(plsda_train_1,
    cutoff = 0.7,
    title = "Correlation plot for PLS-DA model"
)

# you can produce a clustered image map (CIM) to visualise the correlation between variables and samples
quartz()
cim(plsda_train_1,
    comp = 1,
    title = "Clustered Image Map for PLS-DA model"
)


####################################################
################ Cross Validation ##################
####################################################

# use k-fold method with k=4
set.seed(2543)
perf_plsda <- perf(plsda_train_1,
    validation = "Mfold",
    folds = 4, progressBar = FALSE, nrepeat = 10
)

# visualise the performance of the model

plot(perf_plsda,
    col = color.mixo(1:3),
    sd = TRUE, legend.position = "horizontal"
)

# set up a grid of keepX values to test
list_keep_x <- seq(4, 24, by = 2)

# tune the model to find the optimal number of variables to keep
tune_plsda_train_1 <- tune.splsda(
    train_index,
    train_class_labels,
    ncomp = 4,
    validation = "Mfold",
    folds = 4,
    dist = "max.dist",
    progressBar = FALSE,
    measure = "BER",
    test.keepX = list_keep_x,
    nrepeat = 10
)

# we can then extract the average error rate
error <- tune_plsda_train_1$error.rate

ncomp <- tune_plsda_train_1$choice.ncomp$ncomp

ncomp

# calculate the optimal number of variable for the model
select_keep_x <- tune_plsda_train_1$choice.keepX[1:ncomp]

select_keep_x
