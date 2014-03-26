### All functions for the placebo response simulations based on the baseline 
# seizure rate.

# Trial: 2 months (28 days each) of baseline diary, 
# then 4 months of maintanence

## Functions are here for convenience, they will be repeated in each code, where
# they are used. 
 
rm(list=ls())

########################### TRIAL DATASETS #######################################
# creating dataset for baseline diary
baselineDiary <- data.frame(matrix(data=NA, nrow=200, ncol=6))
colnames(baselineDiary)[1:3] <- c("ptID", "baselineMu", "baselineSize")
colnames(baselineDiary)[4:59] <- paste('baselineDay', 1:56, sep="") 
colnames(baselineDiary)[60:61] <- c("sumBaselineMonth1", "sumBaselineMonth2")
colnames(baselineDiary)[62:63] <- c("group", "responder")

# creating dataset for maintanence diary
maintDiary <- data.frame(matrix(data=NA, nrow=200, ncol=115))
colnames(maintDiary)[1:2] <- c("maintMu", "maintSize")
colnames(maintDiary)[3:114] <- paste('maintDay', 1:112, sep="") 
colnames(maintDiary)[115] <- "sumSeizLastMonth"


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

baselineDiary <- baseline(nSubjects=200, nDays=28, nMonths=2)

#################################################################################
########################### MAINTENANCE DIARY ###################################
#################################################################################
# 50% improvement for 20% of subjects in placebo group 
# and for 50% of subjects in treatment group for 4x28 days 
# mean and size of the rows as parameters 
maint <- function(nSubjects, nDays, nMonths){
for (j in 1:nSubjects) {
  if(baselineDiary$responder[j] == 1) {
    diaryDays <- rnbinom(n = nMonths*nDays,  
                         mu = (maintDiary$maintMu[j] <- 0.5*rowMeans(baselineDiary[j, 4:59])),
                         size = (maintDiary$maintSize[j] <- baselineDiary$baselineSize[j]))
    maintDiary[j, 3:114] <- diaryDays
  } else {
    diaryDays <- rnbinom(n = nMonths*nDays,  
                         mu = (maintDiary$maintMu[j] <- rowMeans(baselineDiary[j, 4:59])),
                         size = (maintDiary$maintSize[j] <- baselineDiary$baselineSize[j]))
    maintDiary[j, 3:114] <- diaryDays
    }
  if(baselineDiary$responder[j] == 2) {
    diaryDays <- rnbinom(n = nMonths*nDays,  
                         mu = (maintDiary$maintMu[j] <- 0.5*rowMeans(baselineDiary[j, 4:59])),
                         size = (maintDiary$maintSize[j] <- baselineDiary$baselineSize[j]))
    maintDiary[j, 3:(3+nMonths*nDays-1)] <- diaryDays
  } else {
    diaryDays <- rnbinom(n = nMonths*nDays,  
                         mu = (maintDiary$maintMu[j] <- rowMeans(baselineDiary[j, 4:59])),
                         size = (maintDiary$maintSize[j] <- baselineDiary$baselineSize[j]))
    maintDiary[j, 3:(3+nMonths*nDays-1)] <- diaryDays
  }
# Number of seizures during the last month of maintenance
  maintDiary$sumSeizLastMonth[j] <- rowSums(maintDiary[j,87:114])
  }
  return(maintDiary)
}

maintDiary <- maint(nSubjects=200, nDays=28, nMonths=4)



# calculating # of seizures for the last 28 days of the trial
diaryMatrix[j,150] <- sum(diaryMatrix[j,122:149])
# If # seizures/28d at the end of the study is <= 50% of # of seizures/28d
# at baseline, => responder
diaryMatrix[,151] <- ifelse(diaryMatrix[,150]/2 >= diaryMatrix[,32], 1,0)
### Responder rate:
# Responder - a patient with >= 50% reduction in monthly # of seizures.
# If # seizures/28d at the end of the study is <= 50% of # of seizures/28d
# at baseline, => responder. 
# Outcome is # of responders in the study out of the total sample size 
# Responder Rate: # of responders to placebo / total N in placebo 
respRatePlacebo <- (sum(diaryMatrix[which(diaryMatrix[,33]==0),151])
                    /dim(diaryMatrix[which(diaryMatrix[,33]==0),])[1])

