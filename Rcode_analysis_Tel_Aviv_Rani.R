####################################################################################################
#                                                                                                  #
# Adapted from:                                                                                    #                      
#         DETERMINING THE EFFECTO FO TRAINING DURATION ON THE BEHAVIORAL EXPRESSION OF 
#               HABITIAL CONTROL IN HUMANS: A MULTI-LABORATORY INVESTIVATION
#                                                                                                  #
#                    International consortium for the study of habits                              #
#                                                                                                  #
# * The code was crated by Eva Pool code, validated by Rani Gera and then adapted by Rani Gera for
# his PhD thesis
####################################################################################################


library(car)
library(afex)
library(doBy)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggExtra)
library(BayesFactor)
library(sjstats)
library(jtools)
library(plyr)
library(dplyr)
library(tidyr)
library(metafor)
library(rmcorr)
library(flexmix)
library(psych)
library(emmeans)
library(devtools)
library(effectsize)
library(GPArotation)

#---------------------------------------------------------------------------
#                    PRELIMINARY STUFF 
#---------------------------------------------------------------------------

# Set path
full_path       <- '/Volumes/GoogleDrive/My\ Drive/Experiments/HRS/MULTILAB_HABIT-final'
home_path       <- full_path
figures_path    <- file.path(home_path,'ANALYSIS', 'figures_Rani')
utilities_path  <- file.path(home_path,'ANALYSIS','R')
setwd (home_path)

# source my utilites
source (file.path(utilities_path, 'getChangeIndex.R'))
source (file.path(utilities_path, 'getClassicIndex.R'))
source (file.path(utilities_path, 'makeIndividualDiffPlot.R'))
source (file.path(utilities_path, 'makeSplitGroupPlot.R'))
source (file.path(utilities_path, 'countTrialxCondition.R'))

# get tool
devtools::source_gist("2a1bb0133ff568cbe28d", 
                      filename = "geom_flat_violin.R")

# get database
FULL <- read.delim(file.path(home_path,'DATA/FULL_DATABASE.txt'), header = T, sep ='') # read in dataset
FULL = FULL[FULL$site=='Tel_Aviv',]

# define factors
FULL$site      <- factor(FULL$site)
FULL$ID        <- factor(FULL$ID)
FULL$session   <- factor(FULL$session)
FULL$run       <- factor(FULL$run)
FULL$trial     <- factor(FULL$trial)
FULL$cue       <- factor(FULL$cue)
FULL$prepost   <- factor(FULL$prepost)
FULL$group     <- factor(FULL$group)

# remove the baseline condition from the data
FULL <- subset(FULL, cue == 'Valued' | cue == 'Devalued')

# get the last run of the last training session and all the runs after satiation
DAY1   <- subset(FULL, group == '1-day')
DAY3   <- subset(FULL, group == '3-day')

DAY1 <- ddply(DAY1, .(ID), transform, averagePress  = mean(pressFreq[prepost=="pre"]))
DAY3 <- ddply(DAY3, .(ID), transform, averagePress  = mean(pressFreq[prepost=="pre"]))

C.DAY1 <- subset(DAY1, run == '2' | run == '3')
DAY3   <- subset(DAY3, session == '3') # we want the last day only
C.DAY3 <- subset(DAY3, run == '4' | run == '5')

CHANGE <- rbind(C.DAY1,C.DAY3)


# get variable of interest
CHANGE <- ddply(CHANGE, .(ID), transform, normChangeBehav  = (mean(normPressFreq[prepost=="post" & cue=='Valued']) - mean(normPressFreq[prepost=="pre" & cue=='Valued'])) - (mean(normPressFreq[prepost=="post" & cue=='Devalued']) - mean(normPressFreq[prepost=="pre" & cue=='Devalued'])))
CHANGE <- ddply(CHANGE, .(ID), transform, normChangeLiking = (mean(normLiking[prepost=="post" & cue=='Valued']) - mean(normLiking[prepost=="pre" & cue=='Valued'])) - (mean(normLiking[prepost=="post" & cue=='Devalued']) - mean(normLiking[prepost=="pre" & cue=='Devalued'])))

# code itemxcondition
CHANGE <- ddply(CHANGE, .(ID,prepost), countTrialxCondition)

# get total number of participants included
plyr::count(CHANGE$ID) # note that Caltech2 used a slightly different protocol so there are less repeat per condition

CHANGE$countTrialxCondition <- factor(CHANGE$countTrialxCondition)

# subset by site
C.CALTECH = subset(CHANGE, site == 'Caltech1')
C.CALTECH2= subset(CHANGE, site == 'Caltech2')
C.HAMBURG = subset(CHANGE, site == 'Hamburg')
C.SYDNEY  = subset(CHANGE, site == 'Sydney')
C.TELAVIV = subset(CHANGE, site == 'Tel_Aviv')



#---------------------------------------------------------------------------
#                   MANIPULACTION CHECKS 
#---------------------------------------------------------------------------

#----------------------------- SNACK PLEASANTNESS-------------------------------------

SNACK.means <- aggregate(CHANGE$outcomeliking, by = list(CHANGE$ID, CHANGE$prepost, CHANGE$cue, CHANGE$group, CHANGE$site), FUN='mean') # extract means
colnames(SNACK.means) <- c('ID','prepost','cue','group', 'site','outcomeliking')
SNACK.index <- ddply(SNACK.means, .(ID, prepost), transform, outcomeliking  = outcomeliking-outcomeliking[cue=="Devalued"])
SNACK.index <- subset(SNACK.index, cue!='Devalued')

#--------- tel-aviv
#main
snack.t.stat <- aov_car(outcomeliking ~ group*prepost + Error (ID/prepost), data = subset(SNACK.index, site == 'Tel_Aviv'), anova_table = list(correction = "GG", es = "pes"))
omega_squared(snack.t.stat)

#----------------------------- HUNGER-------------------------------------
HUNGER.means <- aggregate(CHANGE$hunger, by = list(CHANGE$ID, CHANGE$group, CHANGE$prepost,CHANGE$site), FUN='mean') # extract means
colnames(HUNGER.means) <- c('ID','group','prepost','site','hunger')
ggline(HUNGER.means, x = "prepost", y = "hunger", color = "group", add = c("mean_se","jitter"), order = c("pre", "post"), ylab = "Hunger", xlab = "Time relative to devaluation") # plot

#--------- tel-aviv
#main
hunger.t.stat <- aov_car(hunger ~ group*prepost + Error (ID/prepost), data = subset(HUNGER.means, site == 'Tel_Aviv'), anova_table = list(correction = "GG", es = "pes"))
omega_squared(hunger.t.stat)

#---------------------------------------------------------------------------
#                  OUTCOME DEVALUATION CHANGES BY SITE
#---------------------------------------------------------------------------

# ---------------------------- TEL-AVIV  ----------------------------------

# get database
C.TELAVIV = subset(CHANGE, site == 'Tel_Aviv')
TELAVIV.index <- getChangeIndex(C.TELAVIV)# aggregate based on pre-post

# stat
int.t.stat <- aov_car(pressFreq ~ group*cue + Error (ID/cue), data = TELAVIV.index)
omega_squared(int.t.stat)


#---------------------------------------------------------------------------
#                  OUTCOME DEVALUATION CHANGES DISTRUBUTIONS ALL
#---------------------------------------------------------------------------

#---------------------------- FORMAT DATABASE ------------------------------

# first we need the aggregated data
CHANGE.means <- aggregate(CHANGE$normChangeBehav, by = list(CHANGE$ID, CHANGE$group, CHANGE$site), FUN='mean') # extract means
colnames(CHANGE.means) <- c('ID','group','site', 'normChangeBehav')


#---------------------------- FLEXMIX TO IDENTIFY CLUSTERS -----------------

#  what is the number of clusters that better explains the data
n_clusters <- stepFlexmix(normChangeBehav ~ group, data = CHANGE.means, control = list(verbose = 0), k = 1:5, nrep = 200)
getModel(n_clusters, "BIC")

# get cluster size
getModel(n_clusters, which = 1)
getModel(n_clusters, which = 2)
getModel(n_clusters, which = 3)
getModel(n_clusters, which = 4)
getModel(n_clusters, which = 5)

# the we do the analysis specifying the number of cluster we found with step flex
mixlm <- flexmix(normChangeBehav ~ group, data = CHANGE.means, k = 2)

print(table(clusters(mixlm), CHANGE.means$group))
CHANGE.means$Cluster = factor(clusters(mixlm)) # create a variable based on the clustering


#---------------------------- FIGURE 4 -----------------

# rename variables for plot
CHANGE.means$group     <- dplyr::recode(CHANGE.means$group, "1-day" = "Moderate training", "3-day" = "Extensive training" )
CHANGE.means$Cluster   <- dplyr::recode(CHANGE.means$Cluster, "2" = "Outcome-insensitive", "1" = "Outcome-sensitive" )


pp <-  ggplot(CHANGE.means, aes(normChangeBehav, fill = Cluster)) +
  geom_histogram(aes(y=..density..),alpha=0.2,binwidth=0.2)+
  geom_density(alpha = 0.5)+
  xlab('Behavioral adaptation index')+
  ylab('Density')+
  facet_grid(~group)+
  scale_fill_manual(values=c("#F5793A", "#C9C9DB")) +
  theme_bw()

ppp <-  pp + theme_bw(base_size = 17, base_family = "Helvetica")+
  theme(strip.text.x = element_text(size = 14, face = "bold"),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))

ggsave(file.path(figures_path,'histograms_clusters_Tel_Aviv.png'), ppp, dpi = 100)

pdf(file.path(figures_path,'Figure_4_histograms_clusters_Tel_Aviv.pdf'))
print(ppp)
dev.off()



#---------------------------------------------------------------------------
#                       INDIVIDUAL DIFFERENCES 
#---------------------------------------------------------------------------

#----------------------------  DATA REDUCTION  ----------------------------------

# prepare database for the FA
Q_EFA.means.ID <- aggregate(ANXIETY ~ ID * TICS_SOOV * TICS_PREPE * TICS_WODI * TICS_EXWO * TICS_LACK * TICS_SOTE * TICS_SOIS * TICS_WORY * TICS_WOOV * BIS_motor * BIS_attentional * BIS_nonplanning,
                            data = CHANGE, FUN = mean, na.action = na.pass) # we do not include the total scales
Q_EFA.means <- Q_EFA.means.ID
Q_EFA.means$ID <- NULL

# quick look at the covarivance structure
r.subscale = cor(Q_EFA.means, use = "pairwise.complete.obs")
cor.plot(Q_EFA.means,numbers=TRUE,main="correlation matrix")
names(Q_EFA.means)[names(Q_EFA.means) == 'V1'] <- 'STAI'

# check distributions before proceeding with FA
describe (Q_EFA.means)
pairs.panels(na.omit(Q_EFA.means))

# determine the number of factors
nFactor  <- fa.parallel(Q_EFA.means, fm = "ml")


# apply EFA with oblimin
quest.1.efa <- fa(r = Q_EFA.means, nfactors = 4, rotate = "oblimin", fm = "ml")

print(quest.1.efa$loadings,cutoff = 0.0)

# create figure with EFA solution
fa.diagram(quest.1.efa)

# save the plot in the figures folder
dev.print(pdf, file.path(figures_path,'Figure_EFA_oblimin.pdf'))
dev.off()

# calculate the factors loadings
s = factor.scores (Q_EFA.means, quest.1.efa) # 
s

#---------------------------- USE FACTOR AS AS MODERATOR IN THE MAIN ANALYSIS ----------

# merge with the FULL database
axes <- s$scores

# combine it with the participants ID
dat <- cbind(Q_EFA.means.ID, axes)
EFA_CHANGE <- join (CHANGE,dat, type = "full")

# run full model for each factor individually

# stress work
inter.work = lmer(normPressFreq~ group*cue*prepost*ML1 + itemxcondition + (1+cue*prepost+itemxcondition|ID), data = EFA_CHANGE, REML=FALSE, control = lmerControl(optimizer ="bobyqa"))
summary(inter.work)
Confint(inter.work, level = 0.95) 

# ----- assumptions check
plot(fitted(inter.work),residuals(inter.work)) 
qqnorm(residuals(inter.work))
hist(residuals(inter.work))

# stress social
inter.social = lmer(normPressFreq~ group*cue*prepost*ML3 + itemxcondition + (1+cue*prepost+itemxcondition|ID), data = EFA_CHANGE, REML=FALSE, control = lmerControl(optimizer ="bobyqa"))
summary(inter.social)
Confint(inter.social, level = 0.95) 

# ----- assumptions check
plot(fitted(inter.social),residuals(inter.social)) 
qqnorm(residuals(inter.social))
hist(residuals(inter.social))

# stress affect
inter.affect = lmer(normPressFreq~ group*cue*prepost*ML4 + itemxcondition + (1+cue*prepost+itemxcondition|ID), data = EFA_CHANGE, REML=FALSE, control = lmerControl(optimizer ="bobyqa"))
summary(inter.affect)
Confint(inter.affect, level = 0.95) 


# ----- assumptions check
plot(fitted(inter.affect),residuals(inter.affect)) 
qqnorm(residuals(inter.affect))
hist(residuals(inter.affect))

# implusivity
inter.implusivity = lmer(normPressFreq~ group*cue*prepost*ML2 + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = EFA_CHANGE, REML=FALSE,control = lmerControl(optimizer ="bobyqa"))
summary(inter.implusivity)
Confint(inter.implusivity, level = 0.95) 

# ----- assumptions check
plot(fitted(inter.implusivity),residuals(inter.implusivity)) 
qqnorm(residuals(inter.implusivity))
hist(residuals(inter.implusivity))

# test and different points of the model to understand interaction

# Stress affective -1 SD people low on affectiv stress have effect of overtraining
EFA_CHANGE$AFF_pSD <- scale(EFA_CHANGE$ML4, scale = T) + 1 # here I'm going to test at - 1SD (so people that are low in anxiety)
sslop.pSD = lmer(normPressFreq~ group*cue*prepost*AFF_pSD + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = EFA_CHANGE, REML=FALSE,control = lmerControl(optimizer ="bobyqa"))
summary(sslop.pSD)
Confint(sslop.pSD, level = 0.95) 

# Stress Affective +1 SD people high on affective stress have effect of overtraining
EFA_CHANGE$AFF_mSD <- scale(EFA_CHANGE$ML4, scale = T) - 1 # here I'm going to test at + 1SD (so people that are high in anxiety)
sslop.mSD = lmer(normPressFreq ~ group*cue*prepost*AFF_mSD + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = EFA_CHANGE, REML=FALSE, control = lmerControl(optimizer ="bobyqa"))
summary(sslop.mSD)
Confint(sslop.mSD, level = 0.95) 

#---------------------------- FIGURE 5 ---------------------------

# this tests the model predictions as we do in lmer but does not allow to display distributions
AFF.means <- aggregate(EFA_CHANGE$normChangeBehav, by = list(EFA_CHANGE$ID, EFA_CHANGE$group, EFA_CHANGE$site, EFA_CHANGE$AFF_pSD, EFA_CHANGE$AFF_mSD, EFA_CHANGE$ML4), FUN='mean', na.rm = T) # extract means
colnames(AFF.means) <- c('ID','group','site', 'AFF_pSD', 'AFF_mSD','AFF', 'normChangeBehav')

# to assess to imapct of the extream let's run a robust analysis
# AFF.means$normChangeBehav <- rank(AFF.means$normChangeBehav)


# ADJUSTED MEANS in case we want see the estimations from the model
acqC1.aov      <- aov_car(normChangeBehav  ~ group*AFF +Error(ID), data = AFF.means, observed = c("AFF"), factorize = F, anova_table = list(es = "pes"))
acqC1.adjmeans <- emmeans(acqC1.aov, specs = c("group"), by = "AFF", at = list(AFF= c(-1, 1)))
acqC1.adjmeans


acqC1.low.aov      <- aov_car(normChangeBehav  ~ group*AFF_pSD +Error(ID), data = AFF.means, observed = c("AFF"), factorize = F, anova_table = list(es = "pes"))
acqC1.high.aov     <- aov_car(normChangeBehav  ~ group*AFF_mSD +Error(ID), data = AFF.means, observed = c("AFF"), factorize = F, anova_table = list(es = "pes"))


AFF.means$group           <- dplyr::recode(AFF.means$group, "1-day" = "Moderate", "3-day" = "Extensive" )


pp <- ggplot(AFF.means, aes(x = AFF, y = normChangeBehav, fill = group, color = group)) +
  geom_point(alpha = .2, position = position_jitterdodge(jitter.width = .0, jitter.height = 0)) +
  geom_smooth(method = lm, level = .95, alpha = .1, fullrange=TRUE) +
  ylab('Behavioral adaptation index')+
  xlab('Stress Affect')+
  annotate("rect", xmin=0.95, xmax=1.05, ymin=min(AFF.means$normChangeBehav), ymax=max(AFF.means$normChangeBehav), alpha=0.3, fill="gray") +
  annotate("rect", xmin=-0.95, xmax=-1.05, ymin=min(AFF.means$normChangeBehav), ymax=max(AFF.means$normChangeBehav), alpha=0.3, fill="gray") +
  scale_fill_manual(values=c("#56B4E9", "#0F2080")) +
  scale_color_manual(values=c("#56B4E9", "#092C48")) +
  scale_x_continuous(breaks=seq(-2.5,2.5,0.5)) +
  theme_bw()


theme_continous_plot <- theme_bw(base_size = 18, base_family = "Helvetica")+
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 22),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 22))

ppp <- pp + theme_continous_plot


pppp <- ggMarginal(ppp + theme(legend.position = c(1, 1), legend.justification = c(1, 1), legend.background = element_rect(color = "white")), 
                  type = "density", groupFill = T, color = NA, alpha = .2)

pdf(file.path(figures_path,'Figure_5_IndividualDifferences_pannelA.pdf'))
print(pppp)
dev.off()


adj_meanCR  <- c(0.511, 0.127,  0.261, 0.403)
adj_lowerCL <- c(0.3157, -0.0689, 0.0629, 0.2094)
adj_upperCL <- c(0.707, 0.322, 0.459, 0.596)
adj_group   <- c("Moderate", "Extensive", "Moderate", "Extensive")
adj_SD      <- c("Lower Stress Affect (-1 SD)", "Lower Stress Affect (-1 SD)", "Higher Stress Affect (+1 SD)", "Higher Stress Affect (+1 SD)")
adj_means   <- data.frame(adj_meanCR, adj_lowerCL, adj_upperCL, adj_group, adj_SD)

adjmeans_plot <- ggplot(data = adj_means, aes(x = factor(adj_group, levels = c("Moderate","Extensive")), y = adj_meanCR, 
                                              color = adj_group,
                                              fill = adj_group)) + 
  geom_crossbar(aes(y = adj_meanCR, ymin =adj_lowerCL, ymax = adj_upperCL), width = 0.85 , alpha = 0.1) +
  facet_grid(~ factor(adj_SD, levels = c("Lower Stress Affect (-1 SD)","Higher Stress Affect (+1 SD)"))) +
  ylab('Behavioral adaptation index')+
  xlab('Amount of Training')+
  ylim(min= -1.5, max = 3)+
  scale_fill_manual(values=c("#0F2080","#56B4E9" )) +
  scale_color_manual(values=c( "#092C48","#56B4E9")) +
  theme_bw()


theme_means_plots <- theme_bw(base_size = 18, base_family = "Helvetica")+
  theme(strip.text.x = element_text(size = 16, face = "bold"),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))

ppp <- adjmeans_plot+ theme_means_plots 


pdf(file.path(figures_path,'Figure_5_IndividualDifferences_pannelB.pdf'))
print(ppp)
dev.off()
