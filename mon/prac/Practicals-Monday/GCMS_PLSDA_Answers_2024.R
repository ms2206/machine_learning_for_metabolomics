# Qualitative Metabolomics: GC-MS & PLS-DA
# Practical
#
# Maria Anastasiadi

################## CLASSIFICATION #######################################


# Clear workspace
rm(list=ls())

# Close any open graphics devices
graphics.off()


#install.packages("BiocManager")
#BiocManager::install("mixOmics")

# Load additional packages
require(matlab, quietly=TRUE)
#require(plsgenomics)
require(R.matlab)
require(mixOmics)
require(rgl)
require(ptw)
source("pretreat.r")

### DERMATOPHYTOSIS

#Extract variables

X = read.table("data2.csv", sep=",", header=TRUE, row.names=1)
X = as.matrix(X)
Samples = row.names(X)
Sensors = colnames(X)

# Extract class vector from X

class1 = X[,1]
class1 = as.vector(class1)

# extract the remaining data.
RESP = X[,2:ncol(X)]
RESP = as.matrix(RESP)


# Create VARB character vector (from A)
#Sensors = colnames(RESP)


##PCA
pca.enose <- pca(RESP, ncomp = 6, scale = TRUE)
plot(pca.enose)
print(pca.enose)
Var1<-100*(pca.enose$prop_expl_var$X[1])
Var1<-round(Var1, digits=2)
Var2<-100*(pca.enose$prop_expl_var$X[2])
Var2<-round(Var2, digits=2)

##Create a biplot, a 2D PCA and a 3D PLCA score plot

#define class vector
class1 = as.factor(class1)
x11()
plotIndiv(pca.enose, ind.names = Samples, 
          group = class1, style = "lattice")
plotIndiv(pca.enose, ind.names = Samples, 
          group = class1, style = "3d")
biplot(pca.enose, xlab=paste("PC1 ", Var1, "%"), ylab=paste("PC2 ", Var2, "%"), 
       group = class1, col.per.group= c("red", "blue", "orange", "green", "gray"))



############################################### CLASSIFICATION #############################################

## For PLS-DA, train the model
train1<-RESP[c(-25:-21),]
test1<-RESP[21:25,]

class<-as.factor(class1[c(-25:-21)])
testCl <-as.factor(class1[21:25])

plsda.train1 <- plsda(train1, class, ncomp = 10)

x11()
plotIndiv(plsda.train1, ind.names = TRUE, comp = c(1, 2),
          ellipse = TRUE, style = "lattice", cex = c(rep(1, 5)))

#Generate VIPs
train.vip <- vip(plsda.train1)
train.vip<-as.matrix(train.vip)

#Visualise VIPs
barplot(train.vip[,1], beside=TRUE, 
        col = topo.colors(n=24), ylim=c(0, 1.7), xlim=c(0,38))
#correlation plot for variables
plotVar(plsda.train1, cutoff = 0.7)

#Performance of plsda model
set.seed(2543) # for reproducibility 
perf.plsda <- perf(plsda.train1, validation = "Mfold", 
                  folds=4, progressBar = FALSE, nrepeat=10) 

#perf.plsda$auc
x11()
auc.plsda = auroc(plsda.train1, roc.comp = 10)

x11()
plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")


##Set list of variables to be tested on each component
list.keepX <- seq(4,24, 2)
list.keepX # to output the grid of values tested

tune.plsda.train1<-tune.splsda(train1, class, ncomp = 4,
                               validation = "Mfold",
                               folds = 4, dist = 'max.dist', progressBar = FALSE,
                               measure = "BER", test.keepX = list.keepX,
                               nrepeat = 10)

##Extract the errors
error <- tune.plsda.train1$error.rate
ncomp <- tune.plsda.train1$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.plsda.train1$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.plsda.train1, col = color.jet(4))

#Create new optimised model
splsda.train1.opt <- splsda(train1, class, ncomp = ncomp, keepX = select.keepX)

plotIndiv(splsda.train1.opt, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - final result")

#############################################################################
#plsda.train1.opt <- plsda(train1, class, ncomp = 7)
#x11()
#plotIndiv(splsda.train1.opt, ind.names = TRUE, ellipse = TRUE, legend = TRUE)

##Produce a Clustered Image Map
x11()
cim(splsda.train1.opt)

#Warning: For some the visualisation doesn't work. In this case save it as image
x11()
cim(splsda.train1.opt, save = "jpeg", name = "plot.1")

##Find VIP
train.vip1 <- vip(splsda.train1.opt)
train.vip1<-as.matrix(train.vip1)

#Remove empty rows
train.vip1<-train.vip1[rowSums(train.vip1[])>0,]

x11()
barplot(train.vip1[,1],
        beside = TRUE, col=topo.colors(n=24),
        ylim = c(0, 1.7), xlim = c(0, 38), legend = rownames(train.vip1),cex.names = 1,
        main = "Variable Importance in the Projection for LV1", font.main = 4)

barplot(train.vip1[,2],
        beside = TRUE, col=topo.colors(n=24),
        ylim = c(0, 1.7), xlim = c(0, 38), legend = rownames(train.vip1),cex.names = 1,
        main = "Variable Importance in the Projection for LV2", font.main = 4)

# Finally predict with the independent test set
test.predict1 <- predict(splsda.train1.opt, test1, dist = "max.dist")
# store prediction for the 5th component
prediction <- test.predict1$class$max.dist[,2] 

table(factor(prediction, levels=levels(class1)), testCl)

# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = class1[21:25], predicted = prediction)
confusion.mat

###################### GC-MS #############################
# clear memory
rm(list=ls())

#Close graphics
graphics.off()

#Load the required packages
require(R.matlab)
require(matlab)
require(plsgenomics)
require(mixOmics)


# Load data

FILE="NA_BWGT_FAE_CTRL_CD.mat"
# Use readMat from R.matlab package
DATA = readMat(FILE)

# Extract DATA from structure array
XTIC = as.matrix(DATA$XTIC)
colnames(XTIC)<-sprintf("X%s",seq(1:ncol(XTIC)))
CLASS = as.vector(DATA$CLASS)
SAM = as.character(unlist(DATA$SAM))
RT = DATA$RTscan


##PCA OF RAW DATA

pca.raw <- pca(XTIC, ncomp = 4, scale = TRUE)
#plot(pca.raw)

# samples representation
x11()
plotIndiv(pca.raw, ind.names = SAM, 
          group = as.factor(CLASS), style = "lattice", legend = TRUE)
#Compare raw chromatograms for alignment
x11()
matplot(XTIC[1,],type="l",main="Chromatograms", xlab="RT", ylab="TIC" )
matplot(XTIC[22,],type="l",col=2,add=TRUE)
leg.txt = c(paste("chromatogram", SAM[9]), paste("\n chromatogram", SAM[1]))
legend("topleft", leg.txt, pch="_", col=c(1,2), bty="n")


x11()
matplot(XTIC[20,],type="l",main="Chromatograms", xlab="RT", ylab="TIC" )
matplot(XTIC[1,],type="l",col=2,add=TRUE)
leg.txt = c(paste("chromatogram", SAM[6]), paste("\n chromatogram", SAM[1]))
legend("topleft", leg.txt, pch="_", col=c(1,2), bty="n")

### PTW
ref <- XTIC[1,]
samp <- XTIC[2:24,]
gaschrom.ptw <-ptw(ref, samp, warp.type = "individual", verbose = TRUE, 
                   optim.crit = "WCC",  trwdth = 100, init.coef = c(0, 1, 0))
summary(gaschrom.ptw)

##PCA OF PTW DATA
XTIC.ptw <- as.matrix(gaschrom.ptw[["warped.sample"]])
XTIC.ptw <- rbind(XTIC[1,], XTIC.ptw)

library(dplyr)
XTIC.ptw <- XTIC.ptw %>% replace(is.na(.), 0)
pca.ptw <- pca(XTIC.ptw, ncomp = 4, scale = TRUE)
#plot(pca.raw)

# samples representation
x11()
plotIndiv(pca.ptw, ind.names = SAM, 
          group = as.numeric(as.factor(CLASS)), style = "lattice", legend = TRUE)


x11()
matplot(XTIC.ptw[1,],type="l",main="Chromatograms", xlab="RT", ylab="TIC" )
matplot(XTIC.ptw[22,],type="l",col=2,add=TRUE)
leg.txt = c(paste("chromatogram", SAM[9]), paste("\n chromatogram", SAM[1]))
legend("topleft", leg.txt, pch="_", col=c(1,2), bty="n")

XTIC.ptw.out <- XTIC.ptw[-c(7, 24),]

pca.ptw.out <- pca(XTIC.ptw.out, ncomp = 4, scale = TRUE)
#plot(pca.raw)

# samples representation
CLASS1<-as.factor(CLASS[-c(7,24)])
SAM1 <- SAM[-c(7, 24)]
x11()
plotIndiv(pca.ptw.out, ind.names = SAM1, 
          group = as.factor(CLASS1), style = "lattice", legend = TRUE)


##Repeat plsda with aligned chromatograms

# Split into training and test set
cat("\nExtracting test set...\n")
XTrain.ptw = XTIC.ptw.out[-c(3,9,16,19),]
Ctrain.ptw = CLASS1[-c(3,9,16,19)]
SAMtrain.ptw = SAM1[-c(3,9,16,19)]

XTest.ptw = XTIC.ptw.out[c(3,9,16,19),]
Ctest.ptw = CLASS1[c(3,9,16,19)]
SAMtest.ptw = SAM1[c(3,9,16,19)]

# Now perform plsda with mixOmics

plsda.GC_ptw_train <- plsda(XTrain.ptw, Ctrain.ptw, ncomp = 10)
x11()
plotIndiv(plsda.GC_ptw_train, ind.names = TRUE, ellipse = TRUE, legend = TRUE)

#correlation plot for variables
plotVar(plsda.GC_ptw_train, cutoff = 0.7, var.names=FALSE)

#Performance of plsda model
set.seed(2543) # for reproducibility 
perf.plsda.ptw <- perf(plsda.GC_ptw_train, validation = "Mfold", 
                    folds=4, progressBar = FALSE, nrepeat=10) 

x11()
plot(perf.plsda.ptw, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")


##Set list of variables to be tested on each component
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX # to output the grid of values tested

tune.plsda.GC_ptw_train<-tune.splsda(XTrain.ptw, Ctrain.ptw, ncomp = 4,
                                     validation = "Mfold",
                                     folds = 4, dist = 'max.dist', progressBar = FALSE,
                                     measure = "BER", test.keepX = list.keepX,
                                     nrepeat = 10)

##Extract the errors
error <- tune.plsda.GC_ptw_train$error.rate
ncomp <- tune.plsda.GC_ptw_train$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.plsda.GC_ptw_train$choice.keepX[1:2]  # optimal number of variables to select
select.keepX

plot(tune.plsda.GC_ptw_train, col = color.jet(4))

#Create new optimised model
splsda.train.opt1 <- splsda(XTrain.ptw, Ctrain.ptw, ncomp = 2, keepX = select.keepX)

plotIndiv(splsda.train.opt1, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - final result")



# Finally predict
test.predict.ptw <- predict(plsda.GC_ptw_train, newdata=XTest.ptw, dist="max.dist")

# evaluate the prediction accuracy for the first two components
predict.comp2 <- test.predict.ptw$class$max.dist[,2]
#Extract confusion matrix
table(factor(predict.comp2, levels = levels(CLASS1)), Ctest.ptw)

