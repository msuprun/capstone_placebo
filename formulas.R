rm(list=ls())

baselineDiary <- baseline(nSubjects=200, nDays=28, nMonths=2)
maintDiary <- maint(nSubjects=200, nDays=28, nMonths=4)

#################################################################################
########################### EFFECT SIZE #########################################
#################################################################################
### % Reduction of seizures
# 100*(#of seiz/28d at the end - #of seiz/28d at baseline)
# /(#of seiz/28d at baseline).
# average across subjects = effect size 
# only for placebo group
nSeizBase46Days <- rowSums(baselineDiary[which(baselineDiary$group==0),4:59])
nSeizMaint112Days <- rowSums(maintDiary[which(baselineDiary$group==0),3:114])

seizureReducGr0 <- (100*(((nSeizMaint112Days/112)*28) - 28*(nSeizBase46Days/56))
                    / 28*(nSeizBase46Days/56))  
effectSizePlacebo <- mean(seizureReducGr0)

#################################################################################
########################### RESPONDER RATE ######################################
#################################################################################

responder <- ifelse(seizureReducGr0 >= 0.5*(28*(nSeizBase46Days/56)), 1,0)


sum(baselineDiary[2, 4:59])
sum(maintDiary[2, 3:114])
0.5*((40/56)*28) > (37/112)*28


