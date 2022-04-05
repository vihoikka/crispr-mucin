install.packages("glmmTMB", repos="https://glmmTMB.github.io/glmmTMB/repos", type="binary")
install.packages("glmm")
install.packages("TMB")
install.packages("lme4")
install.packages("lmerTest")
library(ggplot2)
library(tidyr)
library(glmm)
library(lme4)
library(lmerTest)

library(glmmTMB)
        
dada <- read.csv('../aeromonas_competition_data.csv', header = TRUE, sep = ";", dec=",")
dada$condRep <- paste(dada$condition, dada$replicate,sep="")
dada$condRep <- as.factor(dada$condRep)
dada_long <- gather(data = dada, key = species, value = titer, columnare:aeromonas)
dada_long$titer <- as.numeric(dada_long$titer)
dada_long_plot <- subset(dada_long, condition != 6)

dada_long$cond_name <- as.factor(dada_long$cond_name)


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#000000", "#CC79A7","#0072B2","#F0E442", "#D55E00")

labels <- c(`Aero + col R` = "Aeromonas sp. + F. col R",
            `Aero + col Z` = "Aeromonas sp. + F. col Z",
            `Aeromonas` = "Aeromonas sp.",
            `columnare R` = "F. columnare R",
            `columnare Z` = "F. columnare Z")

#plot aeromonas titers
p1 <- ggplot(data = dada_long_plot, aes(x = day, y = titer)) +
  geom_point(aes(color = species)) +
  labs(color = "Species") +
  geom_smooth(aes(color = species), method = "loess") +
  facet_grid(~cond_name, labeller = as_labeller(labels)) +
  scale_y_continuous(trans='log10') +
  theme_bw() +
  theme(legend.text = element_text(face = "italic")) +
  scale_color_manual(values=cbPalette,labels = c("Aeromonas sp.", "F. columnare")) +
  ylab("Bacterial titer (cfu/ml)") +
  xlab("Day")
  

p1
ggsave(p1, filename="competition_experiment.png", height = 3, width = 10)

dada_long$cond_name <- relevel(dada_long$cond_name,ref = "Aeromonas") # Declare LW_M as reference
dada_long$condRep <- as.factor(dada_long$condRep)

#using gamma distribution for the model with replicate as possible random effect
# "is the titer of aeromonas different in the presence or absence of columnare (both Z and R), on day 3?

#without random effect the AIC is >2 units smaller
m1 <- glm(titer ~ cond_name,
          family = Gamma("log"),
          data = subset(dada_long, day > 0 & species == "aeromonas" & cond_name != "columnare R" & cond_name != "columnare Z" & cond_name != "neg ctrl" ))
summary(m1)
AIC(m1) #994.1869

#alternative model with random effect
m2 <- glmer(titer ~ cond_name + (1|condRep),
              family = Gamma("log"),
              data = subset(dada_long, day > 0 & species == "aeromonas" & cond_name != "columnare R" & cond_name != "columnare Z" & cond_name != "neg ctrl" ))
summary(m2)
AIC(m2) #998.1869
