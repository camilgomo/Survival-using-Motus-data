# Translating survival results from jags models into daily survival and visualizations gor: Gómez et al. Estimating apparent survival along hemispheric migratory routes using Motus tracking. Movement Ecology.

library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(tidyverse)

dat <- read.csv("Results_Final_AllSpp_Scott.csv", h = T, sep = ";")
head(dat)
surv <- subset(dat, name == 'survival')
surv <- surv %>%
  mutate(Parameter = fct_relevel(Parameter, 
                                 "A-Colombia Coast (11 Lat)", "B-Gulf Coast            (29 Lat)", "C-Great Lakes         (39 Lat)", 
                                 "D-Boreal forest (>49 Lat )", "E-Returned"))

tot2 <- subset(dat, name == 'survival TOTAL 2')
tot1 <-subset(dat, name == 'survival TOTAL 1') 
total <- rbind(tot1, tot2)

det <- subset(dat, name == 'detection probability')

pd <- position_dodge(0.2)
pd2 <- position_dodge(0.5)

#Survival by Lat
windows(12,7)
ggplot(surv, aes(Parameter, Estimate, ymin = c(Estimate - SD), ymax = c(Estimate + SD), colour = Species)) +
  geom_point(position = pd, size=4) + scale_x_discrete(labels = label_wrap_gen(16))+
  geom_errorbar(aes(ymin = c(Estimate - SD), ymax = c(Estimate + SD)), width = 0.15, position = pd)+
  scale_color_manual(values = c("darkgoldenrod1", "cornsilk4", "chocolate4"))+
  labs(x = "",
       y = expression(paste("Survival by latitudinal band", ~ ~ phi)),
       title = "")+ ylim(0.5,1)+
  theme_bw(18) 

#Survival by day
windows(12,7)
ggplot(surv, aes(Parameter, Surv_day, ymin = Surv_day_low, ymax = Surv_day_up, colour = Species)) +
  geom_point(position = pd, size=4) + scale_x_discrete(labels = label_wrap_gen(16))+
  geom_errorbar(aes(ymin = Surv_day_low, ymax = Surv_day_up), width = 0.15, position = pd)+
  scale_color_manual(values = c("darkgoldenrod1", "cornsilk4", "chocolate4"))+
  labs(x = "",
       y = expression(paste("Daily survival", ~ ~ phi)),
       title = "")+ ylim(0.50,1.00)+
  theme_bw(18) 

#Overall Spring Survival 
windows(12,7)
ggplot(total, aes(Species, Estimate, ymin = P2.5., ymax = P97.5., colour = Species, group = name)) +
  geom_point(position = pd, size=4) + scale_x_discrete(labels = label_wrap_gen(21))+
  geom_errorbar(aes(ymin = c(Estimate - SD), ymax = c(Estimate + SD)), width = 0.15, position = pd)+
  scale_color_manual(values = c("darkgoldenrod1", "cornsilk4", "chocolate4"))+
  labs(x = "",
       y = expression(paste("Overall spring migration survival", ~ ~ phi)),
       title = "")+ ylim(0,1)+
  annotate(geom = "text", x = 0.95, y = 0.73, label = "1")+
  annotate(geom = "text", x = 1.05, y = 0.63, label = "2")+
  annotate(geom = "text", x = 1.95, y = 0.67, label = "1")+
  annotate(geom = "text", x = 2.05, y = 0.57, label = "2")+
  annotate(geom = "text", x = 2.95, y = 0.36, label = "1")+
  annotate(geom = "text", x = 3.05, y = 0.29, label = "2")+
  theme_bw(18) 

#Cummulative survival by lat
cum <- dat[dat$Species %in% c("BLPW", "SWTH", "GCTH"), ]
 cum <- cum %>% drop_na(CumSurv)
 
windows(12,7)
ggplot(cum, aes(LAT, CumSurv, ymin = 0, ymax = 1, colour = Species)) +
  geom_point(position = position_dodge(width = 1), size=4) + geom_line(position = position_dodge(width = 1)) +
  geom_errorbar(aes(ymin = CumSurv-CumSurv_SE, ymax = CumSurv+CumSurv_SE), width = 0.15, position = position_dodge(width = 1)) +
  scale_color_manual(values = c("darkgoldenrod1", "cornsilk4", "chocolate4"))+
  labs(x = "Latitude", y = expression(paste("Cummulative survival", ~ ~ phi)), title = "") +
  theme_bw(18) 


#Detection probability
det$Parameter <- recode(det$Parameter, 
                        'recap stage1' = 'Colombia Coast (11 Lat)',
                        'recap stage2' = 'Gulf Coast (29 Lat)',
                        'recap stage3' = 'Great Lakes (39 Lat)',
                        'recap stage4' = 'Breeding Grounds (>49 Lat)',
                        'recap stage5' = 'Returned')

det <- det %>%
  mutate(Parameter = fct_relevel(Parameter, 
                                 
                                 "Colombia Coast (11 Lat)", 
                                 "Gulf Coast (29 Lat)", 
                                 "Great Lakes (39 Lat)", 
                                 "Breeding Grounds (>49 Lat)",
                                 "Returned"))
windows(12,7)
ggplot(det, aes(Parameter, Estimate, ymin = P2.5., ymax = P97.5., colour = Species)) +
  geom_point(position = pd, size=4) + scale_x_discrete(labels = label_wrap_gen(16))+
  geom_errorbar(aes(ymin = P2.5., ymax = P97.5.), width = 0.15, position = pd)+
  scale_color_manual(values = c("darkgoldenrod1", "cornsilk4", "chocolate4"))+
  labs(x = "",
       y = "Latitudinal Detection Probability (p)",
       title = "") +
  theme_bw(15) 

#Rate of migration

dat1 <- dat[dat$ID %in% c("beta[1,5]", "beta[2,5]", "beta[3,5]","alpha[1,1]","alpha[2,1]","alpha[3,1]","alpha[1,2]","alpha[2,2]","alpha[3,2]","alpha[3,1]","alpha[1,3]","alpha[2,3]","alpha[3,3]", "alpha[1,4]","alpha[2,4]","alpha[3,4]"), ]

windows(12,7)
ggplot(dat1, aes(LAT, Mean_Days, colour = Species)) +
  geom_point(position = pd, size=4) +
  geom_line(position = position_dodge(width = 1)) +
  scale_color_manual(values = c("darkgoldenrod1", "cornsilk4", "chocolate4"))+
  ylim(0,20)+
  labs(x = "Latitude (°)",
       y = str_wrap("Mean No.Days taken to fly between latitudinal bands",35),
       title = "")+
  theme_bw(20) 
