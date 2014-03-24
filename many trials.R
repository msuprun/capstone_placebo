# Single trial 
rm(list=ls())
# Initial diary (200 subjects)
diaryMatrix <- matrix(data=NA, nrow=2000, ncol=40)
row.num <- 1
set.seed(-1025)
# i = trail #, j = ID #
for (i in 1:4){
  for (j in 1:500) {
    diaryMatrix[row.num, 1] <- j  # ID
    repeat {
      diary28Days <- rnbinom(n = 28,  
                             mu = (diaryMatrix[row.num, 2] <- runif(1, min=0.1, max=0.9)),
                             size = (diaryMatrix[row.num, 3] <- runif(1, min=1, max=40))) 
      sumVal <- sum(diary28Days)
      # if total N of sezires for 28 dyas < 4 then re-run 
      if(sumVal >= 4) {
        # here the code stacks all iterations of the 'i' loop to a matrix
        diaryMatrix[row.num, 4:31] <- diary28Days
        break
      }
    }
    diaryMatrix[, 32] <- apply(diaryMatrix[, 4:31],1,sum)  # N of seizures per 28 days
    # Group assignment (random)
    diaryMatrix[, 33] <- ifelse(diaryMatrix[,1] %in% (sample(nrow(diaryMatrix), size=0.5*dim(diaryMatrix)[1])), 1,0)
    # 20% responders for placebo group (Column 33 = 0)
    r0 <- diaryMatrix[which(diaryMatrix[,33]==0),][sample(nrow(diaryMatrix[which(diaryMatrix[,33]==0),]), size=(nrow(diaryMatrix[which(diaryMatrix[,33]==0),])*0.2)), 1]
    diaryMatrix[,34] <- ifelse(diaryMatrix[,1] %in% r0, 1, 0)
    # 50% responders for drug group (Column 33 = 1)
    #  r1 <- diaryMatrix[which(diaryMatrix[,33]==1),][sample(nrow(diaryMatrix[which(diaryMatrix[,33]==1),]), size=0.5*dim(diaryMatrix[which(diaryMatrix[,33]==1),])[1]), 1]
    #  diaryMatrix[,34] <- ifelse(diaryMatrix[,1] %in% r1[,1], 2, diaryMatrix[,34])
    # Column #34 - whether subject responded to ttt
    # 1 - response in placebo group
    # 2 - response in drug group
    # 0 - no response in either group 
    
    # outputting trial number 
    diaryMatrix[row.num, 35] <- i
    row.num <- row.num + 1
  }
}

library(lattice)
df <- as.data.frame(diaryMatrix[, c(4:31,32,35)])
df$means <- rowMeans(diaryMatrix[,4:31])
histogram(~ means | V30, data=df)

# 50% improvement for 20% of subjects in placebo group for 4x28 days 
# mean of the rows as a parameter 
for (i in 1:200) {
  if(diaryMatrix[i,34] == 1) {
    diaryMain <- rnbinom(n = 4*28,  
                         mu = (diaryMatrix[i, 35] <- 0.5*mean(diaryMatrix[i, 4:31])),
                         size = (diaryMatrix[i, 36] <- diaryMatrix[i, 3]))
    diaryMatrix[i, 37:148] <- diaryMain
  } else {
    diaryMain <- rnbinom(n = 4*28,  
                         mu = (diaryMatrix[i, 35] <- mean(diaryMatrix[i, 4:31])),
                         size = (diaryMatrix[i, 36] <- diaryMatrix[i, 3]))
    diaryMatrix[i, 37:148] <- diaryMain
  }
}

rm(r0, r1, diary28Days, sumVal, diaryMain)

mean(diaryMatrix[which(diaryMatrix[,34]==1), 4:31])
mean(diaryMatrix[which(diaryMatrix[,34]==1), 37:148])
cbind(diaryMatrix[which(diaryMatrix[,34]==1), 2], 
      diaryMatrix[which(diaryMatrix[,34]==1), 35])
cbind(rowMeans(diaryMatrix[which(diaryMatrix[,34]==1), 4:31]), 
      rowMeans(diaryMatrix[which(diaryMatrix[,34]==1), 37:148]))
cbind(rowMeans(diaryMatrix[which(diaryMatrix[,34]==1), 4:31]), 
      rowMeans(diaryMatrix[which(diaryMatrix[,34]==1), 37:148]))

# Mean to generate data vs actual mean for 28 days 
# 50 % impovement for the observed mean 
cbind(diaryMatrix[which(diaryMatrix[,34]==1), 2], 
      rowMeans(diaryMatrix[which(diaryMatrix[,34]==1), 4:31]))


# Make pretty dataset
diaryPretty <- data.frame(diaryMatrix)
colnames(diaryPretty)[1:3] <- c("ID", "Mu", "Size")
names(diaryPretty)[4:31] <- paste0("Day", 1:28) 
colnames(diaryPretty)[32:36] <- c("SumDiary", "Group", "PlaceboResp", "MuResp", "SizeResp")
names(diaryPretty)[37:148] <- paste0("MDay", 1:112) 

# 10% of pt's with 1 sezire per month 
# 20% responded
# 10% same distribution (not in the responder group) 5 to 35 by 5% 
# everyone else 
# for 1000 times 
# distrinbution (almost normal) of 1000 saples for diaries -> increase samples until get normal


