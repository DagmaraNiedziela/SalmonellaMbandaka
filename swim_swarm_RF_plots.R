library(readxl) 
library(ggplot2) 
library(lsmeans)
library(multcomp)
library(multcompView)
install.packages(c("lsmeans","multcomp","multcompView"))
sessionInfo()

# Plots ####
swim_swarm <- read_excel("swim_swarm_RF_data.xlsx", sheet = 3)
swim_swarm
cfs_numbers <- read.csv("CFS_numbers.csv") 
colnames(cfs_numbers)
cfs_numbers$Isolate_no <- as.character(cfs_numbers$Isolate_no)
# "Original_number" "CFS_number_old"  "CFS_number_new"
colnames(cfs_numbers) <- c("Original_number","Isolate_no","CFS_number_new")
styphi <- c("S. Typhimurium ST4/74","S. Typhimurium ST4/74","S. Typhimurium ST4/74")
smanda <- c("S. Mbandaka NCTC7892","S. Mbandaka NCTC7892","S. Mbandaka NCTC7892")

cfs_numbers2 <- rbind(cfs_numbers,styphi,smanda)

library(dplyr)
swim_swarm2 <- left_join(swim_swarm,cfs_numbers2,by="Isolate_no")
swim_swarm2

ggplot(data = subset(swim_swarm2, Motility_type == "swim")) + 
  geom_boxplot(mapping = aes(x = CFS_number_new, y = Growth, fill = CFS_number_new, group = CFS_number_new)) + 
  facet_wrap(~Temperature) + xlab(" ") + ylab("Growth (mm)") + 
  scale_fill_discrete(name = "Isolate name") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Swim_RF_figure.jpeg")

ggplot(data = subset(swim_swarm2, Motility_type == "swarm")) + 
  geom_boxplot(mapping = aes(x = CFS_number_new, y = Growth, fill = CFS_number_new, group = CFS_number_new)) + 
  facet_wrap(~Temperature) + xlab(" ") + ylab("Growth (mm)") + 
  scale_fill_discrete(name = "Isolate name") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Swarm_RF_figure.jpeg")

# Two-way ANOVA #### 

# swim - no interaction 
swim.aov2_int <- aov(Growth ~ Isolate_no * Temperature, data = subset(swim_swarm, Motility_type == "swim"))
summary(swim.aov2_int)

#                       Df Sum Sq Mean Sq F value   Pr(>F)    
#Isolate_no              7   5234     748   4.405  0.00161 ** 
#  Temperature             1  18624   18624 109.706 7.21e-12 ***
#  Isolate_no:Temperature  7   2603     372   2.191  0.06180 .  
#Residuals              32   5432     170 

# swarm - no interaction 
swarm.aov2_int <- aov(Growth ~ Isolate_no * Temperature, data = subset(swim_swarm, Motility_type == "swarm"))
summary(swarm.aov2_int)

#                       Df Sum Sq Mean Sq F value  Pr(>F)   
#Isolate_no              7  39.25   5.607   3.333 0.00884 **
#  Temperature             1   8.33   8.333   4.954 0.03321 * 
#  Isolate_no:Temperature  7  23.25   3.321   1.974 0.08995 . 
#Residuals              32  53.83   1.682 

# ** Proceed without interaction ####
# *** swim #### 
# Check assumptions 
# homogeneity of variances 
plot(swim.aov2, 1)

install.packages("car")
library(car)
leveneTest(Growth ~ Isolate_no + Temperature, data = subset(swim_swarm, Motility_type == "swim"))
# Package does not work !!!!

# Normality assumption - QQplot 
plot(swim.aov2, 2) 

# Extract the residuals
aov_swim_residuals <- residuals(object = swim.aov2)
# Run Shapiro-Wilk test
shapiro.test(x = aov_swim_residuals ) 
# W = 0.93385, p-value = 0.009478 

# Data is not normally distributed 

# Two way ANOVA
swim.aov2 <- aov(Growth ~ Isolate_no + Temperature, data = subset(swim_swarm, Motility_type == "swim"))
summary(swim.aov2)

#             Df Sum Sq Mean Sq F value   Pr(>F)    
#Isolate_no   7   5234     748   3.629  0.00417 ** 
#  Temperature  1  18624   18624  90.390 1.05e-11 ***
#  Residuals   39   8036     206

# *** swarm #### 
# Assumptions 
plot(swarm.aov2, 1)

install.packages("car")
library(car)
leveneTest(Growth ~ Isolate_no + Temperature, data = subset(swim_swarm, Motility_type == "swarm"))
# Package does not work !!!!

# Normality assumption - QQplot 
plot(swarm.aov2, 2) 

# Extract the residuals
aov_swarm_residuals <- residuals(object = swarm.aov2)
# Run Shapiro-Wilk test
shapiro.test(x = aov_swarm_residuals ) 
# W = 0.76577, p-value = 2.374e-07 

# Not normally distributed 

# Two way ANOVA 
swarm.aov2 <- aov(Growth ~ Isolate_no + Temperature, data = subset(swim_swarm, Motility_type == "swarm"))
summary(swarm.aov2)

#            Df Sum Sq Mean Sq F value Pr(>F)  
#Isolate_no   7  39.25   5.607   2.837 0.0172 *
#  Temperature  1   8.33   8.333   4.216 0.0468 *
#  Residuals   39  77.08   1.976  

#** Multiple comparisons #### 

# *** swim #### 
swim_isolate <- TukeyHSD(swim.aov2, which = "Isolate_no")

Tukey_swim <- as.data.frame(swim_isolate$Isolate_no)
write.csv(Tukey_swim, "Tukey_swim.csv")

library(lsmeans)
library(multcomp)
library(multcompView)

# For Dunnets, set Typhimurium as first as control 
swim_swarm$Isolate_no<-factor(swim_swarm$Isolate_no,levels=c("S. Typhimurium ST4/74","1587078","1443116","1451178","1475735","1475739", "1489873", "S. Mbandaka NCTC7892"),ordered=TRUE)

fit_swim <- lm(Growth ~ Isolate_no + Temperature + Isolate_no:Temperature, data = subset(swim_swarm, Motility_type == "swim"))

fit.tukey_swim <- lsmeans(fit_swim, pairwise ~ Isolate_no | Temperature)[[2]]
fit.tukey_swim
cld(fit.tukey_swim) 
tukey_swim_temp <- as.data.frame(cld(fit.tukey_swim))
tukey_swim_temp
write.csv(tukey_swim_temp, "Tukey_swim_sliced_by_temp.csv") 

# Each strain at two temperatures 
fit.tukey_swim2 <- lsmeans(fit_swim, pairwise ~ Temperature | Isolate_no)[[2]]
fit.tukey_swim2
cld(fit.tukey_swim2) 
tukey_swim_temp2 <- as.data.frame(cld(fit.tukey_swim2))
tukey_swim_temp2
write.csv(tukey_swim_temp2, "Tukey_swim_sliced_by_isolate.csv")

# Compare sliced least squares means via Dunnett's method.
fit.lsm_swim <- lsmeans(fit_swim, "Isolate_no", by=c("Temperature"))
contrast(fit.lsm_swim, "trt.vs.ctrl")

Dunnett_swim_temp <- as.data.frame(contrast(fit.lsm_swim, "trt.vs.ctrl"))
write.csv(Dunnett_swim_temp, "Dunnet_swim_temp.csv") 

# Pairwise t test - between temperatures? 
swim <- subset(swim_swarm, Motility_type == "swim")

# Don't use this! 
pairwise.t.test(swim$Growth, swim$Isolate_no, 
                p.adjust.method = "BH")

# *** swarm #### 
swarm_isolate <- TukeyHSD(swarm.aov2, which = "Isolate_no")

Tukey_swarm <- as.data.frame(swarm_isolate$Isolate_no)
write.csv(Tukey_swarm, "Tukey_swarm.csv")

# For Dunnets, set Typhimurium as first as control - this is done in swim, no need 
swim_swarm$Isolate_no<-factor(swim_swarm$Isolate_no,levels=c("S. Typhimurium ST4/74","1587078","1443116","1451178","1475735","1475739", "1489873", "S. Mbandaka NCTC7892"),ordered=TRUE)

fit_swarm <- lm(Growth ~ Isolate_no + Temperature + Isolate_no:Temperature, data = subset(swim_swarm, Motility_type == "swarm"))

fit.tukey_swarm <- lsmeans(fit_swarm, pairwise ~ Isolate_no | Temperature)[[2]]
fit.tukey_swarm
cld(fit.tukey_swarm) 
tukey_swarm_temp <- as.data.frame(cld(fit.tukey_swarm))
tukey_swarm_temp
write.csv(tukey_swarm_temp, "Tukey_swarm_sliced_by_temp.csv")

# Each strain at two temperatures 
fit.tukey_swarm2 <- lsmeans(fit_swarm, pairwise ~ Temperature | Isolate_no)[[2]]
fit.tukey_swarm2
cld(fit.tukey_swarm2) 
tukey_swarm_temp2 <- as.data.frame(cld(fit.tukey_swarm2))
tukey_swarm_temp2
write.csv(tukey_swarm_temp2, "Tukey_swarm_sliced_by_isolate.csv")

# Compare sliced least squares means via Dunnett's method.
fit.lsm_swarm <- lsmeans(fit_swarm, "Isolate_no", by=c("Temperature"))
contrast(fit.lsm_swarm, "trt.vs.ctrl")

Dunnett_swarm_temp <- as.data.frame(contrast(fit.lsm_swarm, "trt.vs.ctrl"))
write.csv(Dunnett_swarm_temp, "Dunnet_swarm_temp.csv")
