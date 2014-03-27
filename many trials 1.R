# Single trial 
rm(list=ls())

# check how much time it takes for the code to run
# Start the clock!
ptm <- proc.time()

# Initial diary (200 subjects)
diaryMatrix <- matrix(data=NA, nrow=200, ncol=200)
respRatePlacebo <- effectSizePlacebo <- rep(NA, 500)
df <- matrix(NA, 500, 20)
set.seed(-1025)
# i = trail #, j = subject ID 
##### Trial for everyone with >= 4 seizures per month
for (i in 1:500){
  for (j in 1:200) {
    diaryMatrix[j, 1] <- j  # ID
    repeat {
      diaryDays <- rnbinom(n = 2*28,  
                           mu = (diaryMatrix[j, 2] <- runif(1, min=0.1, max=0.9)),
                           size = (diaryMatrix[j, 3] <- runif(1, min=1, max=40))) 
      sumBaselineMonth1 <- sum(diaryDays[1:28])
      sumBaselineMonth2 <- sum(diaryDays[29:46])
      # creating an indicator value to use as a condition for the loop:
      # if eather month1 or month2 (or both) have less than 4 seizures per month
      # the loop will keep running 
      sumVal <- ifelse((sumBaselineMonth1 < 4 | sumBaselineMonth2 < 4) 
                       | (sumBaselineMonth1 < 4 & sumBaselineMonth2 < 4), 1, 0)
      if(sumVal == 0) {
        # here the code stacks all iterations of the 'i' loop to a matrix
        diaryMatrix[j, 4:59] <- diaryDays
        break
      }
    }
    # Group assignment (random)
    diaryMatrix[, 60] <- ifelse(diaryMatrix[,1] %in% (sample(nrow(diaryMatrix), size=0.5*dim(diaryMatrix)[1])), 1,0)
    # 20% responders for placebo group (Column 61 = 0)
    r0 <- (diaryMatrix[which(diaryMatrix[,60]==0),]
           [sample(nrow(diaryMatrix[which(diaryMatrix[,60]==0),]),
                   size=(0.2*nrow(diaryMatrix[which(diaryMatrix[,60]==0),]))),1])
    diaryMatrix[,61] <- ifelse(diaryMatrix[,1] %in% r0, 1, 0)
    # 50% responders for treatment group (responder = 2)
    #  r1 <- (diaryMatrix[which(diaryMatrix[,60]==1),]
    #       [sample(nrow(diaryMatrix[which(diaryMatrix[,60]==1),]), 
    #               size=(nrow(diaryMatrix[whichdiaryMatrix[,60]==1),])*0.5)),1])
    #  diaryMatrix[,61] <- ifelse(diaryMatrix[,1] %in% r1, 2, diaryMatrix[,61])
  }
  
  # 50% improvement for 20% of subjects in placebo group for 4x28 days 
  # mean of the rows as a parameter 
  for (j in 1:200) {
    if(diaryMatrix[j,61] %in% c(1,2)) {
      diaryDays <- rnbinom(n = 4*28,  
                           mu = (diaryMatrix[j, 62] <- 0.5*mean(diaryMatrix[j, 4:31])),
                           size = (diaryMatrix[j, 63] <- diaryMatrix[j, 3]))
      diaryMatrix[j, 64:175] <- diaryDays
    } else {
      diaryDays <- rnbinom(n = 4*28,  
                           mu = (diaryMatrix[j, 62] <- mean(diaryMatrix[j, 4:31])),
                           size = (diaryMatrix[j, 63] <- diaryMatrix[j, 3]))
      diaryMatrix[j, 64:175] <- diaryDays
    }
    # those with total of 0 seizures at baseline were replaced with 0.000001 values, 
    # so we can do division by zero at the effect size calculation
    #  diaryMatrix[j,32] <- ifelse(diaryMatrix[j,32] == 0, 0.0001, diaryMatrix[j,32])
    # dataset for only placebo group
    diaryPlacebo <- diaryMatrix[which(diaryMatrix[,60]==0),]
  }
  
  ### Responder rate:
  responder <- ifelse((rowSums(diaryPlacebo[,64:175])/112)*28 
                      < 0.5*(28*rowSums((diaryPlacebo[,4:59]/56))), 1,0)
  respRatePlacebo[i] <- 100*(sum(responder)/nrow(diaryPlacebo))
  
  
  ### % Reduction of seizures
  nSeizBase56Days <- rowSums(diaryPlacebo[,4:59])
  nSeizMaint112Days <- rowSums(diaryPlacebo[,64:175])
  
  seizureReducPlacebo <- (100*(((nSeizMaint112Days/112)*28) - 28*(nSeizBase56Days/56))
                          / 28*(nSeizBase56Days/56))  
  effectSizePlacebo[i] <- mean(seizureReducPlacebo)
  # dataset with the results only for N simulations 
  df[i,1:2] <- cbind(respRatePlacebo[i], effectSizePlacebo[i])
}



############### Low Baseline Seizure Rate (BSR)
set.seed(-1025)
# choose pct of subjects with low BSR and max number of low seizres (0<nSeiz<4)
pct <- 0.1
nSeizures <- 3
for (i in 1:500){
  for (j in 1:200) {
    diaryMatrix[j, 1] <- j  # ID
    repeat {
      diaryDays <- rnbinom(n = 2*28,  
                           mu = (diaryMatrix[j, 2] <- runif(1, min=0.1, max=0.9)),
                           size = (diaryMatrix[j, 3] <- runif(1, min=1, max=40))) 
      sumBaselineMonth1 <- sum(diaryDays[1:28])
      sumBaselineMonth2 <- sum(diaryDays[29:46])
      sumVal <- ifelse((sumBaselineMonth1 < 4 | sumBaselineMonth2 < 4) 
                       | (sumBaselineMonth1 < 4 & sumBaselineMonth2 < 4), 1, 0)
      if(sumVal == 0) {
      diaryMatrix[j, 4:59] <- diaryDays
      break
      }
    }
    # Those who have 1 in lowBaseSeiz will later have < 4 sezires / 28 days at baseline
    # V176 - those who have low baseline seizure rate will have value of 1
    diaryMatrix[j, 176] <- sample(c(1, 0), 1, 
                                           replace=TRUE, prob=c(pct, (1-pct))) 
    if (diaryMatrix[j, 176] == 1) {
        repeat {
          diaryDays1 <- rnbinom(n=28, mu=(diaryMatrix[j,2] <- runif(1, 0.1, 0.99)), 
                                size=(diaryMatrix[j,3] <- runif(1, 1, 99)))
          val1 <- sum(diaryDays1)
          if (val1 <= nSeizures & val1 >0) { 
            diaryMatrix[j, 4:31] <- diaryDays1
            break
          }
        }
        repeat {
          diaryDays2 <- rnbinom(n=28, mu=(diaryMatrix[j,2] <- runif(1, 0.1, 0.99)), 
                                size=(diaryMatrix[j,3] <- runif(1, 1, 99)))
          val2 <- sum(diaryDays2)
          if (val2 <= nSeizures & val2 > 0) { 
            diaryMatrix[j, 32:59] <- diaryDays2
            break
          }
        }
    # Group assignment (random)
    diaryMatrix[, 60] <- ifelse(diaryMatrix[,1] %in% (sample(nrow(diaryMatrix), size=0.5*dim(diaryMatrix)[1])), 1,0)
    # 20% responders for placebo group (Column 61 = 0)
    r0 <- (diaryMatrix[which(diaryMatrix[,60]==0),]
           [sample(nrow(diaryMatrix[which(diaryMatrix[,60]==0),]),
                   size=(0.2*nrow(diaryMatrix[which(diaryMatrix[,60]==0),]))),1])
    diaryMatrix[,61] <- ifelse(diaryMatrix[,1] %in% r0, 1, 0)
  }
}  
# 50% improvement for 20% of subjects in placebo group for 4x28 days 
# mean of the rows as a parameter 
  for (j in 1:200) {
    if(diaryMatrix[j,61] %in% c(1,2)) {
      diaryDays <- rnbinom(n = 4*28,  
                           mu = (diaryMatrix[j, 62] <- 0.5*mean(diaryMatrix[j, 4:31])),
                           size = (diaryMatrix[j, 63] <- diaryMatrix[j, 3]))
      diaryMatrix[j, 64:175] <- diaryDays
    } else {
      diaryDays <- rnbinom(n = 4*28,  
                           mu = (diaryMatrix[j, 62] <- mean(diaryMatrix[j, 4:31])),
                           size = (diaryMatrix[j, 63] <- diaryMatrix[j, 3]))
      diaryMatrix[j, 64:175] <- diaryDays
    }
    # those with total of 0 seizures at baseline were replaced with 0.000001 values, 
    # so we can do division by zero at the effect size calculation
    #  diaryMatrix[j,32] <- ifelse(diaryMatrix[j,32] == 0, 0.0001, diaryMatrix[j,32])
    # dataset for only placebo group
    diaryPlacebo <- diaryMatrix[which(diaryMatrix[,60]==0),]
  }
  
  ### Responder rate:
  responder <- ifelse((rowSums(diaryPlacebo[,64:175])/112)*28 
                      < 0.5*(28*rowSums((diaryPlacebo[,4:59]/56))), 1,0)
  respRatePlacebo[i] <- 100*(sum(responder)/nrow(diaryPlacebo))
  
  
  ### % Reduction of seizures
  nSeizBase56Days <- rowSums(diaryPlacebo[,4:59])
  nSeizMaint112Days <- rowSums(diaryPlacebo[,64:175])
  
  seizureReducPlacebo <- (100*(((nSeizMaint112Days/112)*28) - 28*(nSeizBase56Days/56))
                          / 28*(nSeizBase56Days/56))  
  effectSizePlacebo[i] <- mean(seizureReducPlacebo)
  # dataset with the results only for N simulations 
  df[i,3:4] <- cbind(respRatePlacebo[i], effectSizePlacebo[i])
}


# selecting 10% of the sample to have only <= 1 seizure per month
diaryMatrix[j,35] <- sample(c(1, 0), 1, replace=TRUE, prob=c(0.1, 0.9))  
#  sampleBad <- diaryMatrix[sample(nrow(diaryMatrix), size=(nrow(diaryMatrix)*0.1)), 1]
#  c <- ifelse(diaryMatrix[,1] %in% sampleBad, 1, 0)
# only one or less seizures per month for those 10%  
if (diaryMatrix[j,35] == 1){
  diary28Days <- rnbinom(n = 28, mu = (diaryMatrix[j, 2] <- 0.01), 
                         size = (diaryMatrix[j, 3] <- runif(1, min=1, max=3))) 
  diaryMatrix[j, 4:31] <- diary28Days
}
diaryMatrix[, 32] <- apply(diaryMatrix[, 4:31],1,sum)  # N of seizures per 28 days

# Stop the clock
t <- proc.time() - ptm
t/60

# Cleaning the workspace
rm(r0, diary28Days, diaryMain, sumBase, ptm, exT, seizureReducPlacebo)

# checking means vs sums for baseline vs last month of the trial
diaryMatrix[21:40, c(33,34, 32,150)]
diaryMatrix[21:40, c(33,34, 2,36)]


# Graphs 
# dataset for plots
df <- as.data.frame(cbind(effectSizePlacebo, respRatePlacebo))
df$rounded <- round(df[,1], 5)

library(ggplot2)
# Effect size over 500 simulations
p <- ggplot(df, aes(x=effSize)) 
(p + geom_histogram(binwidth=0.2, colour="darkgreen", fill="white") + theme_bw() 
 + geom_vline(aes(xintercept=mean(effSize, na.rm=T)), color="magenta", linetype="dashed", size=1)
 + xlab("Effect Size %") + ggtitle("Effect Size for Placebo Group, NSim=500"))

# Responder rate over 500 simulations
q <- ggplot(df, aes(x=respRate)) 
(q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
 + geom_vline(aes(xintercept=mean(respRate, na.rm=T)), color="magenta", linetype="dashed", size=1)
 + xlab("Responder Rate %") + ggtitle("Responder Rate Placebo, NSim=500"))
q <- ggplot(df, aes(x=respRate_le3_30)) 
(q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
 + geom_vline(aes(xintercept=mean(respRate, na.rm=T)), color="magenta", linetype="dashed", size=1)
 + xlab("Responder Rate %") + ggtitle("Responder Rate: <=3 seizures, 30%, NSim=500"))
q <- ggplot(df, aes(x=respRate_le3_20)) 
(q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
 + geom_vline(aes(xintercept=mean(respRate_le3_20, na.rm=T)), color="magenta", linetype="dashed", size=1)
 + xlab("Responder Rate %") + ggtitle("Responder Rate: <=3 seizures, 20%, NSim=500"))
q <- ggplot(df, aes(x=respRate_le3_10)) 
(q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
 + geom_vline(aes(xintercept=mean(respRate_le3_10, na.rm=T)), color="magenta", linetype="dashed", size=1)
 + xlab("Responder Rate %") + ggtitle("Responder Rate: <=3 seizures, 10%, NSim=500"))


