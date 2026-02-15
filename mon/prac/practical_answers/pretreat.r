pretreat <- function(trDATA,tstDATA,METHOD) {
#
# Data pretreatment function to perform the following:
# 1. Mean-centring; 2. Auto-scaling; 3. Range-scaling (0 to 1); 4. Range-scaling (-1 to 1);
# 5. Normalisation; 6. Pareto scaling
#
# trDATA is training data to be scaled.
# tstDATA is test data to be scaled from parameters of trDATA
# if tstDATA or trDATA = 0 then only trDATA or tstDATA are scaled respectively
# METHOD = scaling method (see numbers above)
#
# Function requires two sub-functions which perform the scaling:
# pretreat1 for the training set; pretreat2 for the test set
#


if (METHOD < 1 | METHOD > 6) {
  cat("\nNo Scaling performed...\n")
  trDATAscaled <- trDATA
  tstDATAscaled <- tstDATA
  } else {
  if (length(trDATA) == 1) {
    cat("\nNo scaling performed. Aborting.\n")
    trDATAscaled <- trDATA
    tstDATAscaled <- tstDATA
    } else {
    # scale the training set
    SCALED <- pretreat1(trDATA,METHOD)
    trDATAscaled <- SCALED$SCALED
    PAR1 <- SCALED$PAR1
    PAR2 <- SCALED$PAR2

    # scale the test set
    if (length(tstDATA) == 1) {
      cat("\nNot scaling test dataset..\n")
      tstDATAscaled <- 0
      } else {
      tstDATAscaled <- pretreat2(tstDATA,METHOD,PAR1,PAR2)
      }
    }
  }
  
  # Export parameters
  return(list(trDATAscaled=trDATAscaled, tstDATAscaled=tstDATAscaled,METHOD=METHOD))
}

###########################################
### INTENRAL FUNCTIONS ####################
###########################################

pretreat1 <- function(DATA,METHOD) {
#

m <- dim(DATA)[1]
ONE <- matrix(1,m,1)
if (METHOD == 1) {
  PAR1 <- apply(DATA,2,mean)
  PAR2 <- 0
  TREATED <- (DATA - (ONE %*% PAR1))
  TREATED <- RepNAN(TREATED)

  } else if (METHOD == 2) {
  
  PAR1 <- apply(DATA,2,mean)
  PAR2 <- apply(DATA,2,sd)
  
  TREATED <- (DATA - (ONE %*% PAR1))/(ONE %*% PAR2)
  TREATED <- RepNAN(TREATED)
  
  } else if (METHOD == 3) {
  
  PAR1 <- apply(DATA,2,min)
  PAR2 <- apply(DATA,2,max)
	TREATED <- (DATA - (ONE %*% PAR1))/((ONE %*% PAR2) - (ONE %*% PAR1))
	TREATED <- RepNAN(TREATED)
	
	} else if (METHOD == 4) {
	
	PAR1 <- apply(DATA,2,min)
  PAR2 <- apply(DATA,2,max)
  PARN <- CFMINMAX(PAR1,PAR2)
  PAR1 <- PARN$Y1
  PAR2 <- PARN$Y2
  TREATED <- 2 * (DATA - (ONE %*% PAR1))/((ONE %*% PAR2) - (ONE %*% PAR1)) - 1
	TREATED <- RepNAN(TREATED)
	
  } else if (METHOD == 5) {
  
  #PAR1 <- sqrt(sum(DATA^2))
  PAR1 <- sqrt(apply(DATA^2,2,sum))
  PAR2 <- 0
  TREATED <- (DATA/(ONE %*% PAR1))
  TREATED <- RepNAN(TREATED)
  
  } else if (METHOD == 6) {
  
  PAR1 <- apply(DATA,2,mean)
  PAR2 <- sqrt(apply(DATA,2,sd))
  
  TREATED <- (DATA - (ONE %*% PAR1))/(ONE %*% PAR2)
  TREATED <- RepNAN(TREATED)
  
  } else {
  cat("\nWrong method selected or No Scaling performed.\n")
  TREATED <- DATA
  PAR1 <- 0
  PAR2 <- 0
  
  }

  return(list(SCALED=TREATED,PAR1=PAR1,PAR2=PAR2))
}

####################################################################

pretreat2 <- function(DATA,METHOD,PAR1,PAR2) {
#

m <- dim(DATA)[1]
ONE <- matrix(1,m,1)
if (METHOD == 1) {

  TREATED <- (DATA - (ONE %*% PAR1))
  TREATED <- RepNAN(TREATED)

  } else if (METHOD == 2) {
  
  TREATED <- (DATA - (ONE %*% PAR1))/(ONE %*% PAR2)
  TREATED <- RepNAN(TREATED)
  
  } else if (METHOD == 3) {
  
  TREATED <- (DATA - (ONE %*% PAR1))/((ONE %*% PAR2) - (ONE %*% PAR1))
	TREATED <- RepNAN(TREATED)
	
	} else if (METHOD == 4) {
	
  PARN <- CFMINMAX(PAR1,PAR2)
  PAR1 <- PARN$Y1
  PAR2 <- PARN$Y2
  TREATED <- 2 * (DATA - (ONE %*% PAR1))/((ONE %*% PAR2) - (ONE %*% PAR1)) - 1
	TREATED <- RepNAN(TREATED)
	
  } else if (METHOD == 5) {
  
  TREATED <- (DATA/(ONE %*% PAR1))
  TREATED <- RepNAN(TREATED)
  
  } else if (METHOD == 6) {
  
  TREATED <- (DATA - (ONE %*% PAR1))/(ONE %*% PAR2)
  TREATED <- RepNAN(TREATED)
  
  } else {
  cat("\nWrong method selected or No Scaling performed.\n")
  TREATED <- DATA

  
  }

  return(TREATED)
}

###########################################################################

CFMINMAX <- function(X1,X2) {
#
# Compares whether there are similar minima and maxima.
# If so they are not transformed.
#

EQ <- X1==X2
NEQ <- !EQ
if (sum(EQ) != 0) {
  cat("\nWarning: Some maxima and minima are equal. These inputs will not be transformed.\n")
  Y1 <- (X1*NEQ) - (1*EQ) # Minima
  Y2 <- (X2*NEQ) + (1*EQ) # Maxima
  } else {
  Y1 <- X1
  Y2 <- X2
  }
  return(list(Y1=Y1,Y2=Y2))
}

############################################################################

RepNAN <-function(X) {
#
# Replaces any NaNs with zeros.
#
#
 if (length(which(is.na(X))) > 0) {
 EE <- which(is.na(X))
 X[EE] <- 0 # Replace with zero
 Y <- X
 } else {
 Y <- X
 }
 return(Y)
}
 
