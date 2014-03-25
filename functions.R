### All functions for the placebo response simulations based on the baseline 
# seizure rate.

# Trial: 2 months (28 days each) of baseline diary, 
# then 4 months of maintanence

## Functions are here for convenience, they will be repeated in each code, where
# they are used. 
 
rm(list=ls())

########################### TRIAL DATASET #######################################
# creating dataset for the trial
baselineDiary <- data.frame(matrix(data=NA, nrow=200, ncol=61))
colnames(baselineDiary)[1:3] <- c("ptID", "baselineMu", "baselineSize")
colnames(baselineDiary)[4:59] <- paste('baselineDay', 1:56, sep="") 
colnames(baselineDiary)[60:61] <- c("sumBaselineMonth1", "sumBaselineMonth2")

#################################################################################
########################### BASELINE DIARY ######################################
#################################################################################
### This function creates a matrix for the entire trial and simulates baseline 
# diary values. 
# nSubjects <- total number of subjects per trial
# nDays <- number of days will be considered as "month" (nDays=28)
# nMonths <- number of months required for the diary
# if we need two months of observations = nMonths*nDays (=2*28)
set.seed(-1025)
baseline <- function(nSubjects, nDays, nMonths){ 
# simulating baseline diary  
for (j in 1:nSubjects) {
  baselineDiary$ptID[j] <- j  
# repeating the loop until all subjects have at least 4 seizures per 28 days 
  repeat {
    diaryDays <- rnbinom(n = nMonths*nDays,  
                           mu = (baselineDiary$baselineMu[j] 
                                 <- runif(1, min=0.1, max=0.9)),
                           size = (baselineDiary$baselineSize[j] 
                                   <- runif(1, min=1, max=40))) 
    sumBaselineMonth1 <- sum(diaryDays[1:28])
    sumBaselineMonth2 <- sum(diaryDays[29:46])
# creating an indicator value to use as a condition for the loop:
# if eather month1 or month2 (or both) have less than 4 seizures per month
# the loop will keep running 
    sumVal <- ifelse((sumBaselineMonth1 < 4 | sumBaselineMonth2 < 4) 
                     | (sumBaselineMonth1 < 4 & sumBaselineMonth2 < 4), 1, 0)
    # if total N of sezires for 28 dyas < 4 then re-run 
    if(sumVal == 0) {
      # here the code stacks all iterations of the 'i' loop to a matrix
      baselineDiary[j, 4:59] <- diaryDays
      baselineDiary$sumBaselineMonth1[j] <- sumBaselineMonth1
      baselineDiary$sumBaselineMonth2[j] <- sumBaselineMonth2
      break
      }
    }
  }
  return(baselineDiary)
}

baselineDiary <- baseline(nSubjects=200, nDays=28, nMonths=2)


  # selecting 10% of the sample to have only <= 1 seizure per month
  diaryMatrix[j,60] <- sample(c(1, 0), 1, replace=TRUE, prob=c(0.1, 0.9))  
  #  sampleBad <- diaryMatrix[sample(nrow(diaryMatrix), size=(nrow(diaryMatrix)*0.1)), 1]
  #  c <- ifelse(diaryMatrix[,1] %in% sampleBad, 1, 0)
  # only one or less seizures per month for those 10%  
  if (diaryMatrix[j,60] == 1){
    diary28Days <- rnbinom(n = 2*28, mu = (diaryMatrix[j, 2] <- 0.01), 
                           size = (diaryMatrix[j, 3] <- runif(1, min=1, max=3))) 
    diaryMatrix[j, 4:59] <- diary28Days
  }
  #  diaryMatrix[, 61] <- apply(diaryMatrix[, 4:31],1,sum)  # N of seizures per 28 days first month
  # Group assignment (random)
  ##  diaryMatrix[, 61] <- ifelse(diaryMatrix[,1] %in% (sample(nrow(diaryMatrix), size=0.5*dim(diaryMatrix)[1])), 1,0)
  # 20% responders for placebo group (Column 62 = 0)
  ##  r0 <- diaryMatrix[which(diaryMatrix[,62]==0),][sample(nrow(diaryMatrix[which(diaryMatrix[,33]==0),]), size=(nrow(diaryMatrix[which(diaryMatrix[,33]==0),])*0.2)), 1]
  # Column #62 - whether subject responded to placebo
  ##  diaryMatrix[,62] <- ifelse(diaryMatrix[,1] %in% r0, 1, 0)
  # 50% responders for drug group (Column x = 1)
  #  r1 <- diaryMatrix[which(diaryMatrix[,x]==1),][sample(nrow(diaryMatrix[which(diaryMatrix[,33]==1),]), size=0.5*dim(diaryMatrix[which(diaryMatrix[,33]==1),])[1]), 1]
  #  diaryMatrix[,x] <- ifelse(diaryMatrix[,1] %in% r1[,1], 2, diaryMatrix[,34])
}

# 50% improvement for 20% of subjects in placebo group for 4x28 days 
# mean of the rows as a parameter 
for (j in 1:200) {
  if(diaryMatrix[j,61] == 1) {
    diaryMain <- rnbinom(n = 4*28,  
                         mu = (diaryMatrix[j, 36] <- 0.5*mean(diaryMatrix[j, 4:31])),
                         size = (diaryMatrix[j, 37] <- diaryMatrix[j, 3]))
    diaryMatrix[j, 38:149] <- diaryMain
  } else {
    diaryMain <- rnbinom(n = 4*28,  
                         mu = (diaryMatrix[j, 36] <- mean(diaryMatrix[j, 4:31])),
                         size = (diaryMatrix[j, 37] <- diaryMatrix[j, 3]))
    diaryMatrix[j, 38:149] <- diaryMain
  }
  # those with total of 0 seizures at baseline were replaced with 0.000001 values, 
  # so we can do division by zero at the effect size calculation
  diaryMatrix[j,32] <- ifelse(diaryMatrix[j,32] == 0, 0.0001, diaryMatrix[j,32])
  # calculating # of seizures for the last 28 days of the trial
  diaryMatrix[j,150] <- sum(diaryMatrix[j,122:149])
  # If # seizures/28d at the end of the study is <= 50% of # of seizures/28d
  # at baseline, => responder
  diaryMatrix[,151] <- ifelse(diaryMatrix[,150]/2 >= diaryMatrix[,32], 1,0)
}


### Responder rate:
# Responder - a patient with >= 50% reduction in monthly # of seizures.
# If # seizures/28d at the end of the study is <= 50% of # of seizures/28d
# at baseline, => responder. 
# Outcome is # of responders in the study out of the total sample size 
# Responder Rate: # of responders to placebo / total N in placebo 
respRatePlacebo <- (sum(diaryMatrix[which(diaryMatrix[,33]==0),151])
                    /dim(diaryMatrix[which(diaryMatrix[,33]==0),])[1])


df <- matrix(NA, 1, 2)
x1 <- function(x){
  vals <- rnorm(x)
  mu <- mean(vals)
  std <- sd(vals)
  df[1,1] <-mu
  df[1,2] <- sd
  return(list(sd,mu,vals,df))
}

x1(5)
