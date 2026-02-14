# Regression Methods practical
#
# Tomasz Kurowski, January 2024

# Clear workspace
rm(list=ls())

# Close any open graphics devices
graphics.off()

# Load dependencies
library(caret)
library(MASS)

# Input file paths
enose.file <- "./EnoseAllSamples.csv"
micro.file <- "./MicroCounts.csv"

enose <- read.table(enose.file, header=TRUE, sep=",", row.names=1)
microbial.counts <- read.table(micro.file, header=TRUE, sep =",", row.names=1)

# Merge CFC counts with enose data
all.data <- merge(enose, microbial.counts, by="row.names")
rownames(all.data) <- all.data[,1]

# Remove row names column
all.data <- as.data.frame(all.data[,-1])

# Converting bacterial counts to log10 scale
all.data$CFC <- log10(all.data$CFC)

# Plotting bacterial count distribution
hist(all.data$CFC, breaks=30, freq=F, main="Bacterial count distribution",
     xlab="log10 CFU/g", xlim=c(3, 10))

# Split training and test set
set.seed(90)
train.index <- createDataPartition(all.data$CFC, p = .7,
                                   list = FALSE, times = 1)
data.train <- all.data[train.index,]
data.test <- all.data[-train.index,]

# Fit model using all variables
model.fit <- lm(CFC ~ ., data=data.train)
summary(model.fit)

#Inspecting variable importance
varImp(model.fit)

# Define function which plots prediction for a model
# (better than re-typing the same code over and over!)
plot.prediction <- function(model.fit, data.test) {
  test.predictions <- predict(model.fit, data.test)
  model.rmse <- RMSE(data.test$CFC, test.predictions)
  plot(test.predictions, data.test$CFC, xlab='Predicted log10 CFU/g',
       ylab='Actual log10 CFU/g', main=paste('RMSE:', model.rmse))
  return(test.predictions)
}

# Predict values for test model
test.predictions <- plot.prediction(model.fit, data.test)

# lm with stepwise feature selection
model.fit <- train(CFC ~ ., method='lmStepAIC', data=data.train)
# 1. Plot predicions
plot.prediction(model.fit, data.test)
# 2. Add the x=y line (intercept a=0, slope b=1)
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
# 3. Add the x=y line (intercept a=1, slope b=1)
abline(a = 1, b = 1, col = "blue", lwd = 2, lty = 2)
# 4. Add the x=y line (intercept a=-1, slope b=1)
abline(a = -1, b = 1, col = "blue", lwd = 2, lty = 2)

# k-nearest neighbours
model.fit <- train(CFC ~ ., method='knn', data=data.train,
                   tuneGrid=expand.grid(k=1:20))
plot.prediction(model.fit, data.test)
print(model.fit)

# Plot visualizing the minimisation of the k parameter
plot(model.fit$results$k, model.fit$results$RMSE, xaxt="n",
     ylab="RMSE", xlab="k", main="RMSE for k-Nearest Neighbours regression")
axis(1, at=1:20)

# use k=5

model.fit.opt <- train(CFC ~ ., 
                       data = data.train, 
                       method = "knn", 
                       tuneGrid = data.frame(k = 5))                   
plot.prediction(model.fit.opt, data.test)
#print(model.fit.opt)

plot.prediction(model.fit, data.test)
# 2. Add the x=y line (intercept a=0, slope b=1)
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
# 3. Add the x=y line (intercept a=1, slope b=1)
abline(a = 1, b = 1, col = "blue", lwd = 2, lty = 2)
# 4. Add the x=y line (intercept a=-1, slope b=1)
abline(a = -1, b = 1, col = "blue", lwd = 2, lty = 2)

# OPTIONAL: Loading sensory (classification) data
sensory <- read.table('./SensoryAllSamples.csv', sep=',', header=TRUE,
                      row.names=1)
all.data.sensory <- merge(all.data, sensory, by='row.names')

fresh.counts <- all.data.sensory[(all.data.sensory$sensory==1),]$CFC
semi.fresh.counts <- all.data.sensory[(all.data.sensory$sensory==2),]$CFC
spoiled.counts <- all.data.sensory[(all.data.sensory$sensory==3),]$CFC

# Generate histograms
fresh.h <- hist(fresh.counts, breaks=15, plot=F)
semi.fresh.h <- hist(semi.fresh.counts, breaks=15, plot=F)
spoiled.h <- hist(spoiled.counts, breaks=15, plot=F)

# Identify plot limits
xlim <- c(min(all.data$CFC), max(all.data$CFC))
ylim <- c(0, max(c(fresh.h$density, semi.fresh.h$density, spoiled.h$density)))

# Plot distributions together as one plot
plot(fresh.h$mids, fresh.h$density, main="Bacterial count distribution (by class)",
     type='l', col='red', xlim=xlim, ylim=ylim,
     xlab="log10 CFU/g", ylab="Density")
lines(semi.fresh.h$mids, semi.fresh.h$density, col='green')
lines(spoiled.h$mids, spoiled.h$density, col='blue')
legend('top', legend=c('Fresh', 'Semi-fresh', 'Spoiled'),
       col=c('red', 'green', 'blue'), lty=1)
# OPTIONAL END
