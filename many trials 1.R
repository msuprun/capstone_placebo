rm(list=ls())

# check how much time it takes for the code to run
# Start the clock!
ptm <- proc.time()

# Initial diary (200 subjects)
diaryMatrix <- matrix(data=NA, nrow=200, ncol=200)
respRatePlacebo <- effectSizePlacebo <- rep(NA, 3000)
df <- matrix(NA, 3000, 20)
colnames(df) <- c("respRate", "effSize", "respRate_le3_10", "effSize_le3_10",
                  "respRate_le3_20", "effSize_le3_20", "respRate_le3_30", "effSize_le3_30", 
                  "respRate_le2_10", "effSize_le2_10",  "respRate_le2_20", "effSize_le2_20", 
                  "respRate_le2_30", "effSize_le2_30", "respRate_eq1_10", "effSize_eq1_10", 
                  "respRate_eq1_20", "effSize_eq1_20", "respRate_eq1_30",  "effSize_eq1_30")
set.seed(-1025)
# i = trail #, j = subject ID 
##### Trial for everyone with >= 4 seizures per month
for (i in 1:3000){
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
    if(diaryMatrix[j,61] == 1) {
      repeat {
      diaryDays <- rnbinom(n = 4*28,  
                           mu = (diaryMatrix[j, 62] <- 0.5*mean(diaryMatrix[j, 4:59])),
                           size = (diaryMatrix[j, 63] <- diaryMatrix[j, 3]))
      resp <- ifelse(mean(diaryDays) <= 0.5*mean(diaryMatrix[j, 4:59]), 1, 0)
      if (resp == 1) {
        diaryMatrix[j, 64:175] <- diaryDays
        break
      }
     }
    } else {
      diaryDays <- rnbinom(n = 4*28,  
                           mu = (diaryMatrix[j, 62] <- mean(diaryMatrix[j, 4:59])),
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
                      <= 0.5*(28*rowSums((diaryPlacebo[,4:59]/56))), 1,0)
  respRatePlacebo[i] <- 100*(sum(responder)/nrow(diaryPlacebo))
  
  
  ### % Reduction of seizures
  nSeizBase56Days <- rowSums(diaryPlacebo[,4:59])
  nSeizMaint112Days <- rowSums(diaryPlacebo[,64:175])
  
  seizureReducPlacebo <- (100*(((nSeizMaint112Days/112)*28) - 28*(nSeizBase56Days/56))
                          / 28*(nSeizBase56Days/56))  
  effectSizePlacebo[i] <- median(seizureReducPlacebo)
  # dataset with the results only for N simulations 
  df[i,1:2] <- cbind(respRatePlacebo[i], effectSizePlacebo[i])
}



############### Low Baseline Seizure Rate (BSR)
set.seed(-1025)
# choose pct of subjects with low BSR and max number of low seizres (0<nSeiz<4)
pct <- 0.1
nSeizures <- 3
for (i in 1:3000){
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
    r1 <- diaryMatrix[sample(nrow(diaryMatrix), size=(pct*nrow(diaryMatrix)))]
    diaryMatrix[,176] <- ifelse(diaryMatrix[,1] %in% r1, 1, 0)
    if (diaryMatrix[j, 176] == 1) {
        repeat {
          diaryDays1 <- rnbinom(n=28, mu=(diaryMatrix[j,2] <- runif(1, 0.1, 0.99)), 
                                size=(diaryMatrix[j,3] <- runif(1, 1, 99)))
          sumBaselineMonth1 <- sum(diaryDays1)
          if (sumBaselineMonth1 <= nSeizures & sumBaselineMonth1 >0) { 
            diaryMatrix[j, 4:31] <- diaryDays1
            break
          }
        }
        repeat {
          diaryDays2 <- rnbinom(n=28, mu=(diaryMatrix[j,2] <- runif(1, 0.1, 0.99)), 
                                size=(diaryMatrix[j,3] <- runif(1, 1, 99)))
          sumBaselineMonth2 <- sum(diaryDays2)
          if (sumBaselineMonth2 <= nSeizures & sumBaselineMonth2 > 0) { 
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
    if(diaryMatrix[j,61] == 1) {
      repeat {
        diaryDays <- rnbinom(n = 4*28,  
                             mu = (diaryMatrix[j, 62] <- 0.5*mean(diaryMatrix[j, 4:59])),
                             size = (diaryMatrix[j, 63] <- diaryMatrix[j, 3]))
        resp <- ifelse(mean(diaryDays) <= 0.5*mean(diaryMatrix[j, 4:59]), 1, 0)
        if (resp == 1) {
          diaryMatrix[j, 64:175] <- diaryDays
          break
        }
      }
    } else {
      diaryDays <- rnbinom(n = 4*28,  
                           mu = (diaryMatrix[j, 62] <- mean(diaryMatrix[j, 4:59])),
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
                      <= 0.5*(28*rowSums((diaryPlacebo[,4:59]/56))), 1,0)
  respRatePlacebo[i] <- 100*(sum(responder)/nrow(diaryPlacebo))
  
  
  ### % Reduction of seizures
  nSeizBase56Days <- rowSums(diaryPlacebo[,4:59])
  nSeizMaint112Days <- rowSums(diaryPlacebo[,64:175])
  
  seizureReducPlacebo <- (100*(((nSeizMaint112Days/112)*28) - 28*(nSeizBase56Days/56))
                          / 28*(nSeizBase56Days/56))  
  effectSizePlacebo[i] <- median(seizureReducPlacebo)
  # dataset with the results only for N simulations 
  df[i,3:4] <- cbind(respRatePlacebo[i], effectSizePlacebo[i])
}

# Stop the clock
t <- proc.time() - ptm
t/60

# Number of subjects with low BSR
sum(diaryMatrix[,176])
# Checking if low BSR executed properly (each value should be <= nSeizures)
rowSums(diaryMatrix[which(diaryMatrix[,176]==1),4:31])
rowSums(diaryMatrix[which(diaryMatrix[,176]==1),32:59])

# Cleaning the workspace
rm(pct, r0, sumBaselineMonth1, sumBaselineMonth2, sumVal, 
   t, ptm, nSeizures, i, j, diaryDays, diaryDays1, diaryDays2)

# testing ES and RRate 
colMeans(df[,c(1,3,5,7,9,11,13,15,17,19)])
colMeans(df[,c(2,4,6,8,10,12,14,16,18,20)])
t.test(x=df$effSize, y=df$effSize_le3_30, alternative="two.sided")
t.test(x=df$respRate, y=df$respRate_le3_30, alternative="two.sided")


# Graphs 
# dataset for plots
df <- as.data.frame(df)

library(ggplot2)
# Effect size over 500 simulations
p <- ggplot(df, aes(x=effSize)) 
pp <- ((p + geom_histogram(binwidth=0.2, colour="darkgreen", fill="white") + theme_bw() 
 + geom_vline(aes(xintercept=mean(effSize, na.rm=T)), color="magenta", linetype="dashed", size=1)
 + xlab("Effect Size, %") + ggtitle(">=4 seizures")))
p <- ggplot(df, aes(x=effSize_le3_10)) 
p1 <- ((p + geom_histogram(binwidth=0.2, colour="darkgreen", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(effSize_le3_10, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("Effect Size %") + ggtitle("<4 seizures, 10%")))
p <- ggplot(df, aes(x=effSize_le3_20)) 
p2 <- ((p + geom_histogram(binwidth=0.2, colour="darkgreen", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(effSize_le3_20, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("Effect Size %") + ggtitle("<4 seizures, 20%")))
p <- ggplot(df, aes(x=effSize_le3_30)) 
p3 <- ((p + geom_histogram(binwidth=0.2, colour="darkgreen", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(effSize_le3_30, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("Effect Size, %") + ggtitle("<4 seizures, 30%")))

# Responder rate over 500 simulations
q <- ggplot(df, aes(x=respRate)) 
qq <- ((q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
 + geom_vline(aes(xintercept=mean(respRate, na.rm=T)), color="magenta", linetype="dashed", size=1)
 + xlab("Responder Rate %") + ggtitle(">=4 seizures")))
q <- ggplot(df, aes(x=respRate_le3_10)) 
q1 <- ((q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(respRate_le3_10, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("Responder Rate %") + ggtitle("<4 seizures, 10%")))
q <- ggplot(df, aes(x=respRate_le3_20)) 
q2 <- ((q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(respRate_le3_20, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("Responder Rate %") + ggtitle("<4 seizures, 20%")))
q <- ggplot(df, aes(x=respRate_le3_30)) 
q3 <- ((q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
 + geom_vline(aes(xintercept=mean(respRate_le3_30, na.rm=T)), color="magenta", linetype="dashed", size=1)
 + xlab("Responder Rate, %") + ggtitle("<4 seizures, 30%"))
 + geom_text(aes(x=15, y=650, label=paste("Mean =",format(mean(respRate_le3_30, na.rm=T),digits=5)), 
                 family="Courier", fontface="plain", lineheight=10)))
q3
# Several plots in one window 
library(gridExtra)
grid.arrange(qq,q1,q2,q3, main="Responder Rate Placebo, N = 3000")
grid.arrange(pp,p1,p2,p3, main="Effect Size Placebo, N = 3000")

test <- diaryMatrix[which(diaryMatrix[,61] == 1),]
r.test <- ifelse((rowSums(test[,64:175])/112)*28 
                    <= 0.5*(28*rowSums((test[,4:59]/56))), 1,0)
cbind(rowMeans(test[, 4:59]), rowMeans(test[,64:175])*2, rowMeans(test[,64:175]))