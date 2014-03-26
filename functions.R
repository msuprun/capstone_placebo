### All functions for the placebo response simulations based on the baseline 
# seizure rate.

# Trial: 2 months (28 days each) of baseline diary, 
# then 4 months of maintanence

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
      baselineDiary[j, 4:(4+nMonths*nDays-1)] <- diaryDays
      baselineDiary$sumBaselineMonth1[j] <- sumBaselineMonth1
      baselineDiary$sumBaselineMonth2[j] <- sumBaselineMonth2
      break
      }
  }
# Group assignment: if a value from Std Norm Dist > 0 then group1, group0 otherwise
  baselineDiary$group[j] <- ifelse(rnorm(1) > 0, 1, 0)
# 20% responders for placebo group (responder = 1)
  r0 <- (baselineDiary[which(baselineDiary$group==0),]
         [sample(nrow(baselineDiary[which(baselineDiary$group==0),]), 
                 size=(nrow(baselineDiary[which(baselineDiary$group==0),])*0.2)),1])
  baselineDiary$responder <- ifelse(baselineDiary$ptID %in% r0, 1, 0)
# 50% responders for treatment group (responder = 2)
  r1 <- (baselineDiary[which(baselineDiary$group==1),]
       [sample(nrow(baselineDiary[which(baselineDiary$group==1),]), 
               size=(nrow(baselineDiary[which(baselineDiary$group==1),])*0.5)),1])
  baselineDiary$responder <- ifelse(baselineDiary$ptID %in% r1, 2, baselineDiary$responder)
  }
  return(baselineDiary)
}

#################################################################################
########################### MAINTENANCE DIARY ###################################
#################################################################################
# 50% improvement for 20% of subjects in placebo group 
# and for 50% of subjects in treatment group for 4x28 days 
# mean and size of the rows as parameters 
set.seed(-1025)
maintenance <- function(nSubjects, nDays, nMonths){
for (j in 1:nSubjects) {
  if(baselineDiary$responder[j] %in% c(1,2)) {
    diaryDays <- rnbinom(n = nMonths*nDays,  
                         mu = (maintDiary$maintMu[j] <- 0.5*rowMeans(baselineDiary[j, 4:31])),
                         size = (maintDiary$maintSize[j] <- baselineDiary$baselineSize[j]))
    maintDiary[j, 3:(3+nMonths*nDays-1)] <- diaryDays
  } else {
    diaryDays <- rnbinom(n = nMonths*nDays,  
                         mu = (maintDiary$maintMu[j] <- rowMeans(baselineDiary[j, 4:31])),
                         size = (maintDiary$maintSize[j] <- baselineDiary$baselineSize[j]))
    maintDiary[j, 3:(3+nMonths*nDays-1)] <- diaryDays
  }
# Number of seizures during the last month of maintenance
  maintDiary$sumSeizLastMonth[j] <- rowSums(maintDiary[j,87:114])
 }
  return(maintDiary)
}


#################################################################################
########################### SMALL N OF BASELINE SEIZURES ########################
#################################################################################

## Selecting % of the sample to have only small N (0<n<4) of seizures per month.
# selection is made by drawing at random from (1,0) with specified propbabilities
# (pct). If 1, then subject will have low baseline seizure rate, 0 - otherwise
## Each subject should have at least 1 seizure
set.seed(-1025)
lowSeizure <- function(pct, nSeizures, nSubjects, nMonths, nDays){
  for (j in 1:nSubjects){
# Those who have 1 in lowBaseSeiz will later have < 4 sezires / 28 days at baseline
  baselineDiary$lowBaseSeiz[j] <- sample(c(1, 0), 1, 
                                         replace=TRUE, prob=c(pct, (1-pct))) 
  if (baselineDiary$lowBaseSeiz[j] == 1) {
  if (nSeizures %in% c(1,2,3)) {
    repeat {
      diaryDays1 <- rnbinom(n=1*nDays, mu=(baselineDiary$baselineMu[j] <- runif(1, 0.1, 0.99)), 
                           size=(baselineDiary$baselineSize[j] <- runif(1, 1, 99)))
      val1 <- sum(diaryDays1)
      if (val1 <= nSeizures & val1 >0) { 
        baselineDiary[j, 4:31] <- diaryDays1
        break
        }
      }
    repeat {
      diaryDays2 <- rnbinom(n=1*nDays, mu=(baselineDiary$baselineMu[j] <- runif(1, 0.1, 0.99)), 
                            size=(baselineDiary$baselineSize[j] <- runif(1, 1, 99)))
      val2 <- sum(diaryDays2)
      if (val2 <= nSeizures & val2 > 0) { 
        baselineDiary[j, 32:59] <- diaryDays2
        break
      }
    }
    } else {
      stop("nSeizures should be 1, 2, or 3")
    }
  } 
}
  return(baselineDiary)
}
  
#  sampleBad <- diaryMatrix[sample(nrow(diaryMatrix), size=(nrow(diaryMatrix)*0.1)), 1]
#  c <- ifelse(diaryMatrix[,1] %in% sampleBad, 1, 0) 


