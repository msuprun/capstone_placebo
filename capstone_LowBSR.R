# check how much time it takes for the code to run
# Start the clock!
ptm <- proc.time()

nSimulations <- 5000
nSubjects <- 200

####################################################################################
####################### All subjects have > 3 seizures #############################
####################################################################################

# Initial diary (200 subjects)
diaryMatrix <- matrix(data=NA, nrow=nSubjects, ncol=200)
respRatePlacebo <- effectSizePlacebo <- rep(NA, nSubjects*0.5)
df <- matrix(NA, nSimulations, 2)

set.seed(-1025)
# i = trail #, j = subject ID 
##### Trial for everyone with >= 4 seizures per month
for (i in 1:nSimulations){
  for (j in 1:nSubjects) {
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
  for (j in 1:nSubjects) {
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
  df[i,] <- cbind(respRatePlacebo[i], effectSizePlacebo[i])
}


####################################################################################
####################### All subjects have < 4 seizures #############################
####################################################################################

# nSimulations = number of simulations
# nSubjects = number of subjects per simulations (trial)
# nSeizures = number of sezires 0 < n < 4 
# pct = percent of subjects with low BSR 

LowBSR <- function(nSimulations, nSubjects, nSeizures, pct) {
  
  a <- length(nSeizures)
  b <- length(pct)
  
  # dataset for a single trial  
  diaryMatrix <- matrix(NA, nrow=nSubjects, ncol=200)
  
  # dataset for reponder rate (RR) and effect size (ES)   
  respRatePlacebo <- effectSizePlacebo <- numeric(nSimulations)
  
  # dataset for the results of all simulations  
  stats <- rep(0, 4)
  
  # z = nSeizures; y = pct 
  for (z in 1:a) {
    for (y in 1:b) {
      
      for (i in 1:nSimulations){
        for (j in 1:nSubjects) {
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
          r1 <- diaryMatrix[sample(nrow(diaryMatrix), size=(pct[y]*nrow(diaryMatrix)))]
          diaryMatrix[,176] <- ifelse(diaryMatrix[,1] %in% r1, 1, 0)
          if (diaryMatrix[j, 176] == 1) {
            repeat {
              diaryDays1 <- rnbinom(n=28, mu=(diaryMatrix[j,2] <- runif(1, 0.1, 0.99)), 
                                    size=(diaryMatrix[j,3] <- runif(1, 1, 99)))
              sumBaselineMonth1 <- sum(diaryDays1)
              if (sumBaselineMonth1 <= nSeizures[z] & sumBaselineMonth1 >0) { 
                diaryMatrix[j, 4:31] <- diaryDays1
                break
              }
            }
            repeat {
              diaryDays2 <- rnbinom(n=28, mu=(diaryMatrix[j,2] <- runif(1, 0.1, 0.99)), 
                                    size=(diaryMatrix[j,3] <- runif(1, 1, 99)))
              sumBaselineMonth2 <- sum(diaryDays2)
              if (sumBaselineMonth2 <= nSeizures[z] & sumBaselineMonth2 > 0) { 
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
        for (j in 1:nSubjects) {
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
        
        newrow <- cbind(y, z, respRatePlacebo[i], effectSizePlacebo[i])
        stats <- rbind(stats, newrow)
      }
    }
  }
  result <- as.data.frame(stats[-1,], row.names=FALSE)
  return(result)
}

pct <- c(0.3, 0.2, 0.1)
nSeizures <- c(3, 2, 1)
set.seed(-1025)
stats <- LowBSR(nSimulations=nSimulations, nSubjects=nSubjects, 
                nSeizures=nSeizures, pct=pct)

# Stop the clock
t <- proc.time() - ptm

#stats <- as.data.frame(bob, row.names=FALSE)
stats$trialNum <- 1:nSimulations
stats$name <- paste0(stats$y, sep="_", stats$z)
stats$y <- stats$z <- NULL
# Long to wide format 
library(reshape)
wide <- reshape(stats, timevar = "name",
                idvar = "trialNum", direction = "wide")

# Final dataset 
statsFinal <- cbind(wide, df)
# Chnaging variable names 
colnames(statsFinal) <- c("trialNum","respRate_le3_30", "effSize_le3_30", "respRate_le3_20", "effSize_le3_20", 
                  "respRate_le3_10", "effSize_le3_10", "respRate_le2_30", "effSize_le2_30",  
                  "respRate_le2_20", "effSize_le2_20", "respRate_le2_10", "effSize_le2_10", 
                  "respRate_eq1_30", "effSize_eq1_30", "respRate_eq1_20", "effSize_eq1_20", 
                  "respRate_eq1_10",  "effSize_eq1_10", "respRate_gt3_all", "effSize_gt3_all")
# Re-arranging column order
statsFinal <- statsFinal[,c(1,20,6,4,2,12,10,8,18,16,14,
                            21,7,5,3,13,11,9,19,17,15)]

# saving stats dataset as csv
write.csv(statsFinal, "stats_final_4-24-2014.csv", row.names=FALSE)
# cleaning workspace (leaving only summary datasets)
#rm(list=setdiff(ls(), c("df", "stats", "nSeizures", 
#                       "nSimulations", "diaryPlacebo", "t")))

# Summary stats
summary(statsFinal[,2:11])  # Responder rate vars
summary(statsFinal[,12:21])  # Effect size vars

# creating pretty descriptive tables 
library(sjPlot)
sjt.df(statsFinal[,2:11], alternateRowColors=TRUE,
       title="Table 1. Responder Rate, N Simulations = 5000", 
       showCommentRow=TRUE, 
       commentString="Responder Rate for different percentages (10,20,30%) of low baseline seizure rates: 
       1. more than 4 seizures per month, 2. less than 4 seizures per month,
       3. less than 3 seizures per month, 4. less than 2 seizures per month")
sjt.df(statsFinal[,12:21], alternateRowColors=TRUE,
       title="Table 2. Effect Size, N Simulations = 5000", 
       showCommentRow=TRUE, 
       commentString="Effect Size for different percentages (10,20,30%) of low baseline seizure rates: 
       1. more than 4 seizures per month, 2. less than 4 seizures per month,
       3. less than 3 seizures per month, 4. less than 2 seizures per month")

# Plots
library(ggplot2)
dev.off()  # to avoid graphics device issues 
# Responder rate over 5000 simulations
q <- ggplot(statsFinal, aes(x=respRate_gt3_all)) 
qq <- ((q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(respRate_gt3_all, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("") + ggtitle(">=4 seizures, 100%"))
        + geom_text(aes(x=15.5, y=1900, label=paste("Mean =",format(mean(respRate_gt3_all, na.rm=T),digits=3)))))
q <- ggplot(statsFinal, aes(x=respRate_le3_10)) 
q1 <- ((q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(respRate_le3_10, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("Responder Rate %") + ggtitle("<=3 seizures, 10%"))
        + geom_text(aes(x=15.5, y=1450, label=paste("Mean =",format(mean(respRate_le3_10, na.rm=T),digits=3)))))
q <- ggplot(statsFinal, aes(x=respRate_le3_20)) 
q2 <- ((q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(respRate_le3_20, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("") + ggtitle("<=3 seizures, 20%"))
        + geom_text(aes(x=16, y=1250, label=paste("Mean =",format(mean(respRate_le3_20, na.rm=T),digits=3)))))
q <- ggplot(statsFinal, aes(x=respRate_le2_30)) 
q3 <- ((q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(respRate_le2_30, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("Responder Rate %") + ggtitle("<=2 seizures, 30%"))
        + geom_text(aes(x=18.5, y=950, label=paste("Mean =",format(mean(respRate_le2_30, na.rm=T),digits=3)))))
q <- ggplot(statsFinal, aes(x=respRate_eq1_20)) 
q4 <- ((q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(respRate_eq1_20, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("") + ggtitle("1 seizure, 20%"))
        + geom_text(aes(x=18.5, y=800, label=paste("Mean =",format(mean(respRate_eq1_20, na.rm=T),digits=3)))))
q <- ggplot(statsFinal, aes(x=respRate_eq1_30)) 
q5 <- ((q + geom_histogram(binwidth=1, colour="midnightblue", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(respRate_eq1_30, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("Responder Rate, %") + ggtitle("1 seizure, 30%"))
        + geom_text(aes(x=19, y=650, label=paste("Mean =",format(mean(respRate_eq1_30, na.rm=T),digits=3)))))

# Effect size over 5000 simulations
p <- ggplot(statsFinal, aes(x=effSize_gt3_all)) 
pp <- ((p + geom_histogram(binwidth=0.2, colour="darkgreen", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(effSize_gt3_all, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("") + ggtitle(">=4 seizures"))
        + geom_text(aes(x=-3, y=750, label=paste("Mean =",format(mean(effSize_gt3_all, na.rm=T),digits=3)))))
p <- ggplot(statsFinal, aes(x=effSize_le3_10)) 
p1 <- ((p + geom_histogram(binwidth=0.2, colour="darkgreen", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(effSize_le3_10, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("Effect Size %") + ggtitle("<=3 seizures, 10%"))
        + geom_text(aes(x=-2.5, y=1100, label=paste("Mean =",format(mean(effSize_le3_10, na.rm=T),digits=3)))))
p <- ggplot(statsFinal, aes(x=effSize_le3_20)) 
p2 <- ((p + geom_histogram(binwidth=0.2, colour="darkgreen", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(effSize_le3_20, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("") + ggtitle("<=3 seizures, 20%"))
        + geom_text(aes(x=-1.8, y=1750, label=paste("Mean =",format(mean(effSize_le3_20, na.rm=T),digits=3)))))
p <- ggplot(statsFinal, aes(x=effSize_le2_30)) 
p3 <- ((p + geom_histogram(binwidth=0.2, colour="darkgreen", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(effSize_le2_30, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("Effect Size, %") + ggtitle("<=2 seizures, 30%"))
        + geom_text(aes(x=-1.6, y=3350, label=paste("Mean =",format(mean(effSize_le2_30, na.rm=T),digits=3)))))
p <- ggplot(statsFinal, aes(x=effSize_eq1_20)) 
p4 <- ((p + geom_histogram(binwidth=0.2, colour="darkgreen", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(effSize_eq1_20, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("") + ggtitle("1 seizure, 20%"))
        + geom_text(aes(x=-1.7, y=3200, label=paste("Mean =",format(mean(effSize_eq1_20, na.rm=T),digits=3)))))
p <- ggplot(statsFinal, aes(x=effSize_eq1_30)) 
p5 <- ((p + geom_histogram(binwidth=0.2, colour="darkgreen", fill="white") + theme_bw() 
        + geom_vline(aes(xintercept=mean(effSize_eq1_30, na.rm=T)), color="magenta", linetype="dashed", size=1)
        + xlab("Effect Size, %") + ggtitle("1 seizure, 30%"))
        + geom_text(aes(x=-1.1, y=4150, label=paste("Mean =",format(mean(effSize_eq1_30, na.rm=T),digits=3)))))

# Several plots in one window 
library(gridExtra)
grid.arrange(qq,q1,q2,q3,q4,q5, main="Responder Rate, N = 5000", nrow=2)
grid.arrange(pp,p1,p2,p3,p4,p5, main="Effect Size, N = 5000", nrow=2)


