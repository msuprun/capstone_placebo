# Single trial 
rm(list=ls())
# Initial diary (200 subjects)
diaryMatrix <- matrix(data=NA, nrow=200, ncol=200)

set.seed(-1025)
# i = trail #, j = ID #
for (j in 1:200) {
  diaryMatrix[j, 1] <- j  # ID
  repeat {
    diary28Days <- rnbinom(n = 2*28,  
                           mu = (diaryMatrix[j, 2] <- runif(1, min=0.1, max=0.9)),
                           size = (diaryMatrix[j, 3] <- runif(1, min=1, max=40))) 
    sumBase1 <- sum(diary28Days[1:28])
    sumBase2 <- sum(diary28Days[29:46])
    sumVal <- ifelse(sumBase1 < 4 | sumBase2 < 4, 1, 0)
# if total N of sezires for 28 dyas < 4 then re-run 
    if(sumVal == 0) {
# here the code stacks all iterations of the 'i' loop to a matrix
      diaryMatrix[j, 4:59] <- diary28Days
      break
    }
  }
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

rm(r0, diary28Days, diaryMain, sumBase)

# see how many had 1 seizre 
rowSums(diaryMatrix[which(diaryMatrix[,60] == 1),4:31])
rowSums(diaryMatrix[which(diaryMatrix[,60] == 1),4:31])

### Responder rate:
# Responder - a patient with >= 50% reduction in monthly # of seizures.
# If # seizures/28d at the end of the study is <= 50% of # of seizures/28d
# at baseline, => responder. 
# Outcome is # of responders in the study out of the total sample size 
# Responder Rate: # of responders to placebo / total N in placebo 
respRatePlacebo <- (sum(diaryMatrix[which(diaryMatrix[,33]==0),151])
                    /dim(diaryMatrix[which(diaryMatrix[,33]==0),])[1])


### % Reduction of seizures
# 100*(#of seiz/28d at the end - #of seiz/28d at baseline)
# /(#of seiz/28d at baseline).
# average across subjects = effect size 
# only for placebo group
seizureReducPlacebo <- round(100*((diaryMatrix[which(diaryMatrix[,33]==0),150] 
                             - diaryMatrix[which(diaryMatrix[,33]==0),32])
                            /(diaryMatrix[which(diaryMatrix[,33]==0),32])),5)  
effectSizePlacebo <- mean(seizureReducPlacebo)





# Mean to generate data vs actual mean for 28 days 
# 50 % impovement for the observed mean 
cbind(diaryMatrix[which(diaryMatrix[,34]==1), 2], 
      rowMeans(diaryMatrix[which(diaryMatrix[,34]==1), 4:31]))


# Make pretty dataset
diaryPretty <- data.frame(diaryMatrix)
colnames(diaryPretty)[1:3] <- c("ID", "Mu_Diary", "Size_Diary")
names(diaryPretty)[4:31] <- paste0("Day", 1:28) 
colnames(diaryPretty)[32:36] <- c("SumDiary", "Group", "PlaceboResp", "Mu_Maintenance", "Size_Maintenance")
names(diaryPretty)[37:148] <- paste0("MDay", 1:112) 

# 10% of pt's with 1 sezire per month 
# 20% responded
# 10% same distribution (not in the responder group) 5 to 35 by 5% 
# everyone else 
# for 1000 times 
# distrinbution (almost normal) of 1000 saples for diaries -> increase samples until get normal
# calculate response rate for each patient (mean before and after)
#responder rate = number of pts who respondand percent seizure reduction
# dif b/w number per month at base line and at the end / number suzires 

    