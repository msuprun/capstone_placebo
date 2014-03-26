# Outcome measure calculations for epilpsy trials 

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
# Responder - a patient with >= 50% reduction in monthly # of seizures.
# Outcome is # of responders in the study out of the total sample size 
# Responder Rate: # of responders to placebo (V151) / total N in placebo = proportion

responder <- ifelse((nSeizMaint112Days/112)*28 < 0.5*(28*(nSeizBase46Days/56)), 1,0)
responderRate <- 100*(sum(responder)/
                        nrow(baselineDiary[which(baselineDiary$group==0),]))
  

