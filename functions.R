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
