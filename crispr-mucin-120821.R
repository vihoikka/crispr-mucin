# Scripts for analysing and plotting OD data from unpublished long term coevolution study (Almeida 2022). Code by Ville Hoikkala.

#### Load data and libraries ####
library(data.table)
dada <- fread('longerm_stats_vh_220620_updated.csv') #read as data table

dada_21 <- fread('longtermdata_spacers_270421.csv') #read as data table
data_21_titers <- fread('aeromonas_titers_020621.csv',na.strings=c("","NA")) #read as data table

#data_21 modification
dada_21$totalspacers <- dada_21$spacers_C1 + dada_21$spacers_C2
dada_21$rep <- substr(dada_21$colony,1,2)
dada_21$booleanNewSpacersVerbal <- ifelse(dada_21$totalspacers == 0, "=0", "More than 0")
dada_21$booleanNewSpacers <- ifelse(dada_21$totalspacers == 0, FALSE, TRUE)

dada_21_subset <- subset(dada_21, Condition == "columnare_aeromonas_phage" | Condition == "columnare_phage")

dada$sample <-factor(dada$sample, ordered = FALSE )
dada$plate_phage <-factor(dada$plate.phage)
dada$plate_nophage <-factor(dada$plate.nophage)

#data_21 titers modification
data_21_titers_long <- melt(data_21_titers, 
                            measure.vars=c("titer_rough","titer_rhizoid","titer_total"),
                            variable.name="titer_type"
)

data_21_titers_long <- na.omit(data_21_titers_long)

#phage titers
titers <- fread('phage_titers_280520.csv') #read as data table
titers$z_prop <- titers$bac_z/titers$bac_tot
titers$r_prop <- titers$bac_r/titers$bac_tot
titers$bac_phage_difference <- titers$phage_titer - titers$bac_tot

#spacer accumulation
spacerTime <- fread('spacer_accumulation.csv')

#bacterial isolate times
columns <- c("sample","")

#mutation file
mutDada <- fread('mutCounts_070520.csv') #read as data table

#Library preparations
install.packages("ggplot2")
install.packages("glmmTMB", repos="https://glmmTMB.github.io/glmmTMB/repos", type="binary")
install.packages("splines")
install.packages("reshape2")
install.packages("tidyr")
install.packages("dplyr")
install.packages("patchwork")
install.packages("ggrepel")
install.packages("ggpubr")
install.packages("fitdistrplus")
install.packages("gamlss")
install.packages("AICcmodavg")

library(ggplot2)
theme_update(plot.title = element_text(hjust = 0.5)) #make all headings centerized
library(glmmTMB)
library(splines)
library(reshape2)
library(tidyr)
library(dplyr)
library(patchwork)
library(ggrepel)
library(ggpubr)
library(fitdistrplus)
library(gamlss)          # defines pdf, cdf of ZIP
library(AICcmodavg)
library(dplyr)



#### Arrange data ####
#Remove NA's, transform to long format, calculate time relative to ancestor, create datasets for spacer analysis, subset in different samples

# Only use isolates with OD data
dadaOD <- na.omit(dada, cols=c("max_od.phage")) #NOTE: na.omit with cols requires a data table, not a data frame
dadaOD$OD_diff_10.phage <- dadaOD$OD_diff.phage + 10 #add 10 to OD difference values
dadaOD$OD_diff_10.nophage <- dadaOD$OD_diff.nophage + 10

dadaOD$sample <- relevel(dadaOD$sample,ref = "LW_M") # Declare LW_M as reference
dadaOD$sample <- relevel(dadaOD$sample,ref = "B245") # Declare B245 as reference


#Long format transformation
dadaOD_long <- dadaOD %>%
  gather (key, value, 'plate.phage','plate.nophage','max_od.phage','max_od.nophage','OD_diff_10.phage', 'OD_diff_10.nophage', 'time.phage', 'time.nophage', 'timeDiff.phage', 'timeDiff.nophage') %>%
  separate(key, into = c("measurement","phage"), sep = "\\.") %>%
  spread(measurement, value)

#change "phage" and "nophage" to "yes" and "no"
dadaOD_long$phage[dadaOD_long$phage == "phage"] <- "yes"
dadaOD_long$phage[dadaOD_long$phage == "nophage"] <- "no"
dadaOD_long$plate <- as.factor(dadaOD_long$plate)

#Add a relative time column (compare time to ancestor). B245 with phage 3064, without phage 3024
dadaOD_long$timeDiff10k <- dadaOD_long$timeDiff+10000

#Add extraction time column (information extracted from sample name using regexpr [everything after dot character])
dadaOD_long$extractTime <- as.integer(regmatches(dadaOD_long$id, regexpr("[^.]*$", dadaOD_long$id, perl=TRUE)))

# Trim for tests that use spacers
na.omit(subset(dadaOD_long, sample =="B245"), cols = "id")
test <- na.omit(dadaOD_long, cols = "id")

dadaODSpacers1long <- drop_na(dadaOD_long,C1_spacers)
dadaODSpacersAlllong <- drop_na(dadaODSpacers1long,C2_spacers)
dadaODSpacersAlllong$totalSpacers <- dadaODSpacersAlllong$C1_spacers + dadaODSpacersAlllong$C2_spacers

# Subset into different treatments
LW <- subset(dadaOD_long, sample == "LW")
LW_M <- subset(dadaOD_long, sample == "LW_M")
SHIEH <- subset(dadaOD_long, sample == "Shieh")
SHIEH_M <- subset(dadaOD_long, sample == "Shieh_M")

dadaMucin <- subset(dadaOD_long, sample == "Shieh_M" | sample == "LW_M")
dadaNoShieh <- subset(dadaOD_long, sample == "Shieh_M" | sample == "LW_M" | sample == "LW")


samplesNames <- c(
  `LW` = "Lake water",
  `LW_M` = "Lake water + mucin",
  `Shieh` = "Shieh",
  `Shieh_M` = "Shieh + mucin"
)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#000000", "#CC79A7","#0072B2","#F0E442", "#D55E00")
linetypes <- c("solid","dashed","twodash","dotted")

#### Titers (Fig 2) ####

phage_titers <- ggplot(data = titers[!is.na(titers$phage_titer)], aes(x=week,y=phage_titer)) +
  geom_line(aes(color = rep, linetype = rep), size = 0.8, alpha = 0.9) +
  #  geom_line(aes(x=week, y=bac_tot, color = rep, linetype = rep, group = rep), size = 1, alpha = 0.3) +
  scale_y_continuous(trans='log10') +
  ylab("Phage (pfu/ml)") +
  facet_grid(~sample,  
             labeller = as_labeller(samplesNames)) +
  scale_colour_manual(values=cbPalette) +
  theme_bw() +
  xlab("")

bac_titers <- ggplot(data = titers[!is.na(titers$bac_tot), ], aes(x=week,y=bac_tot)) +
  geom_line(aes(color = rep, linetype = rep, group = rep), size = 0.8, alpha = 0.9) +
  scale_y_continuous(trans='log10') +
  ylab("Bacteria (cfu/ml)") +
  facet_grid(~sample,
             labeller = as_labeller(samplesNames)) +
  scale_colour_manual(values=cbPalette) +
  scale_linetype_manual(values=linetypes) +
  theme_bw()

titerplot <- phage_titers / bac_titers + plot_annotation(tag_levels = 'A')
ggsave(titerplot, filename = "manuscript_figs/fig_2_titerplot.png", width = 8.5, height = 5)

# Phage titer models LW and LWM
titer_stats_LW_all <- subset(titers, (sample == "LW" | sample == "LW_M") & week > 1 & ctrl == "no")
titer_stats_LW <- subset(titers, (sample == "LW" | sample == "LW_M") & week >= 9 & ctrl == "no")
library(lme4)
library(glmmTMB)

model_phagetiter_LW_LWM_wk9_16 <- glmmTMB(data = titer_stats_LW, phage_titer ~ sample, family = Gamma("log"))
summary(model_phagetiter_LW_LWM_wk9_16)
AIC(model_phagetiter_LW_LWM_wk9_16)

#Visualize results
SuppFig_LW_phagetiters <- ggplot(data = titer_stats_LW_all, aes(x = sample, y = phage_titer)) +
  geom_boxplot() +
  geom_point(shape = 21, size = 2) +
  facet_grid(~week) +
  scale_y_continuous(trans='log10') +
  ylab("Phage titer (pfu/ml)") +
  xlab("Week") +
  theme(axis.text.x = element_text(angle = 45))

#Bacterial titer models in LW and LWM
model_bacterial_titer_LW_LWM_wk9_16 <- glm(data = titer_stats_LW_all, bac_tot ~ sample)
summary(model_bacterial_titer_LW_LWM_wk9_16)
AIC(model_bacterial_titer_LW_LWM_wk9_16)

#Visualize results
SuppFig_LW_bactiters <- ggplot(data = titer_stats_LW_all, aes(x = sample, y = bac_tot)) +
  geom_boxplot() +
  geom_point(shape = 21, size = 2) +
  facet_grid(~week) +
  scale_y_continuous(trans='log10') +
  ylab("Bacterial titer (cfu/ml)") +
  xlab("Week") +
  theme(axis.text.x = element_text(angle = 45))

#Save as supplementary figure
SuppFig_LW_titers <- SuppFig_LW_phagetiters / SuppFig_LW_bactiters + plot_annotation(tag_levels = "A")
ggsave(SuppFig_LW_titers, file = "manuscript_figs/SI/SuppFig_2_LW_titers.png", height = 8, width = 10.9)

#### Spacers and morphologies (Fig 3) ####
spacerTime$noC1spacers <- spacerTime$screened - spacerTime$C1_spacers
spacerTime$noC2spacers <- spacerTime$screened - spacerTime$C2_spacers

spacerTimeLongC1 <- melt(spacerTime, 
                         measure.vars=c("C1_spacers","noC1spacers"),
                         variable.name="locus",
                         value.name="spacers"
)

spacerTimeLongC1$noSpacers <- spacerTimeLongC1$screened - spacerTimeLongC1$spacers

spacerTimeLongC2 <- melt(spacerTime, 
                         measure.vars=c("C2_spacers","noC2spacers"),
                         variable.name="locus",
                         value.name="spacers"
)

spacerTimeLongC2$noSpacers <- spacerTimeLongC2$screened - spacerTimeLongC2$spacers

spacerTimeLongC1Per <- spacerTimeLongC1  %>%
  group_by(week, sample, locus) %>%
  summarise(n = sum(spacers)) %>%
  mutate(percentage = n / sum(n))

spacerTimeLongC2Per <- spacerTimeLongC2  %>%
  group_by(week, sample, locus) %>%
  summarise(n = sum(spacers)) %>%
  mutate(percentage = n / sum(n))

#Only screened ones
spacerTimeLongC1_onlyScreened <- subset(spacerTimeLongC1, screened != 0)
spacerTimeLongC2_onlyScreened <- subset(spacerTimeLongC2, screened != 0)

spacerTimeLongC1Per_screened <- spacerTimeLongC1_onlyScreened  %>%
  group_by(week, sample, locus) %>%
  summarise(n = sum(spacers)) %>%
  mutate(percentage = n / sum(n))

spacerTimeLongC2Per_screened <- spacerTimeLongC2_onlyScreened  %>%
  group_by(week, sample, locus) %>%
  summarise(n = sum(spacers)) %>%
  mutate(percentage = n / sum(n))


ylim.prim <- c(0, 1)
ylim.sec <- c(0, 30)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

#Combined C1 and C2 spacers over time
spacers_time_plot3 <- ggplot(data = spacerTimeLongC1Per_screened, aes(x = week, y = percentage, fill = locus)) +
  geom_area(data = subset(spacerTimeLongC1Per_screened, locus == "C1_spacers"), aes(x = week, y = percentage), alpha=0.6 , size=0.5, colour="black", position = position_stack(reverse = F)) +
  geom_area(data = subset(spacerTimeLongC2Per_screened, locus == "C2_spacers"), aes(x = week, y = percentage), alpha=0.6 , size=0.5, colour="black", position = position_stack(reverse = F)) +
  facet_grid(~sample,  
             labeller = as_labeller(samplesNames)) +
  scale_fill_manual(values = c("#E69F00","#56B4E9","white"),name = "", labels = c("II-C spacers","VI-B spacers","No spacers"),guide = guide_legend(reverse=TRUE)) +
  xlab("Weeks") +
  theme_bw() +
  scale_y_continuous("Proportion of colonies with \n spacer acquisition (area graph)", sec.axis = sec_axis(~ (. - a)/b, name = "Number of \n colonies screened (bars)")) +
  geom_col(data = spacerTimeLongC1_onlyScreened, aes(y=a+screened*b), alpha = 0.08, fill = "red",position = "dodge") +
  theme( axis.line.y.right = element_line(color = "#D55E00"), 
         axis.ticks.y.right = element_line(color = "#D55E00"),
         axis.text.y.right = element_text(color = "#D55E00"),
         axis.title.y.right = element_text(color = "#D55E00"),
         axis.text.x = element_text(angle = 0))

#Relevel to put LW first
dadaODSpacersAlllongSpacerPanel <- dadaODSpacersAlllong
dadaODSpacersAlllongSpacerPanel$sample <- relevel(dadaODSpacersAlllongSpacerPanel$sample,ref = "LW") # Declare B245 as reference

plot_spacersSamplePanel <- ggplot(data = subset(dadaODSpacersAlllongSpacerPanel, ctrl == "no" & ancestral == "no"), aes(x=sample, y=totalSpacers)) +
  geom_jitter(width=0.2, alpha = 0.5, height = 0, size = 1.5, shape = 21, aes(fill = Morphology)) +
  xlab("Treatment") + ylab("Total spacers \n (II-C + VI-B)") +
  scale_y_continuous(breaks = seq(0, 10, len = 11)) +
  scale_fill_manual(values=cbPalette) +
  theme(panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=c("LW","LW+M","Shieh","Shieh+M"))

#Morphologies
dadaOD_long$sample <- relevel(dadaOD_long$sample,ref = "LW")
percentMorphotypes <- subset(dadaOD_long, ctrl == "no" & ancestral == "no") %>% group_by(sample) %>% count(Morphology) %>%
  mutate(ratio=scales::percent(n/sum(n)))
plot_morphotypeProportions <- ggplot(data = subset(dadaOD_long, ctrl == "no" & ancestral == "no"), aes(x=sample,fill=Morphology)) + 
  geom_bar(position = "fill") +
  ggtitle("Morphotypes") +
  scale_fill_manual(values=cbPalette) +
  ylab("Proportion") +
  xlab("Treatment") +
  theme_bw()

#Make mutations long (syn vs nonsyn)
mutData_long <- melt(mutDada, 
                     measure.vars=c("base_substitution.synonymous","base_substitution.nonsynonymous"),
                     variable.name="mutation_type",
                     value.name="count"
)


plot_mutTypes <- ggplot(data = mutData_long, aes(x=Condition, y=count, group = sample, fill = mutation_type)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 45)) +
  ylab("# of mutations") +
  facet_grid(~Ctrl) +
  scale_fill_manual(values=cbPalette) +
  scale_x_discrete(labels=c("LW","LW+M","Shieh","Shieh+M"))
plot_mutTypes

plot_mutTypes <- ggplot(data = mutData_long, aes(x=Condition, y=count, group = sample, fill = Ctrl)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge2(width=.5,preserve = "single"),color="black") +
  theme(axis.text.x = element_text(angle = 45)) +
  ylab("# of mutations") +
  xlab("Treatment") +
  ggtitle("Mutations") +
  facet_grid(~mutation_type,
             labeller = as_labeller(c("base_substitution.synonymous" = "Synonymous", "base_substitution.nonsynonymous" = "Non-synonymous"))) +
  scale_fill_manual(values=cbPalette, name = "Control") +
  scale_x_discrete(labels=c("LW","LW+M","Shieh","Shieh+M"))

plot_mutTypes

fig3a <- spacers_time_plot3
fig3bc <- (plot_spacersSamplePanel + plot_morphotypeProportions) + plot_spacer() + plot_layout(widths = c(2,2,1))
fig3 <- fig3a / fig3bc + plot_annotation(tag_levels = 'A')
ggsave(fig3, filename = "manuscript_figs/fig3.png", height = 5.5, width = 11)

#### Follow-up MAX OD (Fig 4) ####
dadaOD_long$sample <- relevel(dadaOD_long$sample,ref = "B245") # Declare B245 as reference

#Plots
plot_timeVsPhage <- ggplot(data = dadaOD_long, aes(x=sample, y=time)) +
  ggtitle("Time to reach max OD") +
  geom_boxplot(width=0.5, outlier.alpha = 0) +
  geom_jitter(aes(colour = ctrl), alpha = 0.5, width = 0.1, height = 0) +
  xlab("Phage") + ylab("Time relative to ancestor (min)") +
  facet_grid(~phage) +
  theme(legend.position = "none") 

fig_4A <- ggplot(data = subset(dadaOD_long, ctrl == "no"), aes(x=sample, y=max_od)) +
  geom_boxplot(width=0.5, outlier.alpha = 0) +
  geom_jitter(data = subset(dadaOD_long, ctrl == "no"), aes(fill = Morphology), alpha = 0.8, width = 0.1, height = 0, shape = 21) +
  labs(shape="") +
  xlab("Treatment") + ylab("ODmax") +
  ggtitle("Phage resistance") +
  theme(axis.title.x = element_blank()) +
  facet_grid(~phage,
             labeller = as_labeller(c(
               "yes" = "Phage",
               "no" = "No phage")
             )) +
  scale_x_discrete(labels=c("Ancestor","LW","LW+M","Shieh","Shieh+M")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=cbPalette)

fig_4A_time <- ggplot(data = subset(dadaOD_long, ctrl == "no"), aes(x=sample, y=time)) +
  geom_boxplot(width=0.5, outlier.alpha = 0) +
  geom_jitter(data = subset(dadaOD_long, ctrl == "no"), aes(fill = Morphology), alpha = 0.8, width = 0.1, height = 0, shape = 21) +
  labs(shape="") +
  xlab("Treatment") + ylab("Time to ODmax") +
  ggtitle("Phage resistance") +
  theme(axis.title.x = element_blank()) +
  facet_grid(~phage,
             labeller = as_labeller(c(
               "yes" = "Phage",
               "no" = "No phage")
             )) +
  scale_x_discrete(labels=c("Ancestor","LW","LW+M","Shieh","Shieh+M")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=cbPalette)

#Controls
fig_4B <- ggplot(data = subset(dadaOD_long, sample == "B245" | ctrl == "yes"), aes(x=phage, y=max_od)) +
  geom_boxplot(width=0.5, outlier.alpha = 0) +
  geom_line(aes(group = sample, linetype = sample, color = sample), size = 0.8, alpha = 0.8, stat  = "summary") +
  geom_jitter(aes(fill = sample), alpha = 0.7, width = 0.1, height = 0, size = 3, shape = 21) +
  labs(shape="") +
  ggtitle("Control samples") +
  xlab("") + ylab("ODmax") +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels=c("No phage","Phage")) +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  theme_bw()

fig_4B_time <- ggplot(data = subset(dadaOD_long, sample == "B245" | ctrl == "yes"), aes(x=phage, y=time)) +
  geom_boxplot(width=0.5, outlier.alpha = 0) +
  geom_line(aes(group = sample, linetype = sample, color = sample), size = 0.8, alpha = 0.8, stat  = "summary") +
  geom_jitter(aes(fill = sample), alpha = 0.7, width = 0.1, height = 0, size = 3, shape = 21) +
  labs(shape="") +
  ggtitle("Control samples") +
  xlab("") + ylab("Time to ODmax") +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels=c("No phage","Phage")) +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  theme_bw()

fig4 <- fig_4A + fig_4B +plot_annotation(tag_levels = "A")
ggsave(fig4, filename = "manuscript_figs/fig4.png", width = 9, height = 4)

fig4_time <- fig_4A + fig_4A_time + fig_4B + fig_4B_time + plot_annotation(tag_levels = "A")
ggsave(fig4_time, filename = "manuscript_figs/fig4_time.png", width = 9, height = 6)


#### Spacers / morphology (Fig 5) and associated models ####

phageFacetNames <- c(
  `yes` = "With phage",
  `no` = "No phage"
)

hello <- subset(dadaODSpacersAlllong, sample == "LW_M")

booleanSpacerPlot_maxOD_LWM <- ggplot(data = subset(dadaODSpacersAlllong, sample == "LW_M"), 
                                      aes(x = booleanNewSpacers, y = max_od, fill = Morphology, label = id)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(height = 0, width = 0.05, aes(fill = Morphology), size = 2, alpha = 0.5, shape = 21) +
  scale_shape_manual(values=c(21,22)) +
  scale_fill_manual(values=cbPalette) +
  labs(
    x = "New spacers",
    y = "Max OD"
  ) +
  ylim(c(0.13,0.75)) +
  facet_grid(~phage,
             labeller = as_labeller(phageFacetNames)) +
  ggtitle("Lake water + mucin") +
  xlab("") +
  scale_x_discrete(labels=c("No new\nspacers","One or more\nnew spacers"))

booleanSpacerPlot_maxOD_LWM_time <- ggplot(data = subset(dadaODSpacersAlllong, sample == "LW_M"), 
                                           aes(x = booleanNewSpacers, y = time, fill = Morphology, label = id)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(height = 0, width = 0.05, aes(fill = Morphology), size = 2, alpha = 0.5, shape = 21) +
  scale_shape_manual(values=c(21,22)) +
  scale_fill_manual(values=cbPalette) +
  labs(
    x = "New spacers",
    y = "TIme to max OD"
  ) +
  facet_grid(~phage,
             labeller = as_labeller(phageFacetNames)) +
  ggtitle("Lake water + mucin") +
  xlab("") +
  scale_x_discrete(labels=c("No new\nspacers","One or more\nnew spacers"))


booleanSpacerPlot_maxOD_ShiehM <- ggplot(data = subset(dadaODSpacersAlllong, sample == "Shieh_M"), 
                                         aes(x = booleanNewSpacers, y = max_od, fill = Morphology, label = id)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(height = 0, width = 0.05, aes(fill = Morphology), size = 2, alpha = 0.5, shape = 21) +
  scale_shape_manual(values=c(21,22)) +
  scale_fill_manual(values=cbPalette) +
  labs(
    x = "",
    y = "Max OD"
  ) +
  ylim(c(0.13,0.75)) +
  facet_grid(~phage,
             labeller = as_labeller(phageFacetNames)) + 
  ggtitle("Shieh + mucin") +
  scale_x_discrete(labels=c("No new\nspacers","One or more\nnew spacers"))

booleanSpacerPlot_maxOD_ShiehM_time <- ggplot(data = subset(dadaODSpacersAlllong, sample == "Shieh_M"), 
                                              aes(x = booleanNewSpacers, y = time, fill = Morphology, label = id)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(height = 0, width = 0.05, aes(fill = Morphology), size = 2, alpha = 0.5, shape = 21) +
  scale_shape_manual(values=c(21,22)) +
  scale_fill_manual(values=cbPalette) +
  labs(
    x = "",
    y = "Time to max OD"
  ) +
  facet_grid(~phage,
             labeller = as_labeller(phageFacetNames)) + 
  ggtitle("Shieh + mucin") +
  scale_x_discrete(labels=c("No new\nspacers","One or more\nnew spacers"))


spacerPlots <- booleanSpacerPlot_maxOD_LWM / booleanSpacerPlot_maxOD_ShiehM + plot_annotation(tag_levels = "A")
ggsave(spacerPlots, filename = "manuscript_figs/fig_5.png", width = 7.5, height = 6)

spacerPlots_time <- (booleanSpacerPlot_maxOD_LWM | booleanSpacerPlot_maxOD_LWM_time) / (booleanSpacerPlot_maxOD_ShiehM | booleanSpacerPlot_maxOD_ShiehM_time) + plot_annotation(tag_levels = "A")
ggsave(spacerPlots_time, filename = "manuscript_figs/fig_5_time.png", width = 10, height = 7)

#Models for fig 5
#LWM with phage
spacers_morphotype_LWM_phage_model <- glmmTMB(max_od ~ booleanNewSpacers * Morphology+ (1|plate) + (1|rep), 
                                              family = Gamma("log"),
                                              data = subset(dadaODSpacersAlllong, sample == "LW_M" & phage == "yes"))
summary(spacers_morphotype_LWM_phage_model)
AIC(spacers_morphotype_LWM_phage_model)

#LWM without phage
spacers_morphotype_LWM_noPhage_model <- glmmTMB(max_od ~ booleanNewSpacers * Morphology + (1|plate) + (1|rep), 
                                                family = Gamma("log"), 
                                                data = subset(dadaODSpacersAlllong, sample == "LW_M" & phage == "no"))
summary(spacers_morphotype_LWM_noPhage_model)
AIC(spacers_morphotype_LWM_noPhage_model)

hist(subset(dadaODSpacersAlllong, sample == "LW_M" & phage == "no")$max_od)

#ShiehM with phage
spacers_morphotype_ShiehM_phage_model <- glmmTMB(max_od ~ booleanNewSpacers * Morphology + (1|plate) + (1|rep), 
                                                 family = Gamma("log"),  
                                                 data = subset(dadaODSpacersAlllong, sample == "Shieh_M" & phage == "yes"))
summary(spacers_morphotype_ShiehM_phage_model)
AIC(spacers_morphotype_ShiehM_phage_model)

#ShiehM without phage
spacers_morphotype_ShiehM_noPhage_model <- glmmTMB(max_od ~ booleanNewSpacers * Morphology + (1|plate) + (1|rep), 
                                                   family = Gamma("log"), 
                                                   data = subset(dadaODSpacersAlllong, sample == "Shieh_M" & phage == "no"))
summary(spacers_morphotype_ShiehM_noPhage_model)
AIC(spacers_morphotype_ShiehM_noPhage_model)

#ShiehM with and without phagae combined
spacers_morphotype_ShiehM_pooled_model <- glmmTMB(max_od ~ booleanNewSpacers * Morphology * phage+ (1|plate) + (1|rep) + (1|id), 
                                                  family = Gamma("log"), 
                                                  data = subset(dadaODSpacersAlllong, sample == "Shieh_M"))
summary(spacers_morphotype_ShiehM_pooled_model)
AIC(spacers_morphotype_ShiehM_pooled_model)


#### Aeromonas experiment (Fig 6) and associated models ####
spacerdata_21_totalspacers_plot <- ggplot(data = dada_21_subset, aes(x=Condition, y=totalspacers)) +
  geom_boxplot() +
  geom_jitter(aes(colour = morphology), alpha = 0.5, width = 0.1, height = 0) +
  ylab("New spacers in a colony") +
  scale_x_discrete(labels=c("Aeromonas\n(+)","Aeromonas\n(-)")) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  ggtitle("Total spacers (II-C + VI-B)") +
  ylim(0,11) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
spacerdata_21_totalspacers_plot

spacerdata_21_totalspacers_plot_C1 <- ggplot(data = dada_21_subset, aes(x=Condition, y=spacers_C1)) +
  geom_boxplot() +
  geom_jitter(aes(colour = morphology), alpha = 0.5, width = 0.2, height = 0) +
  ylab("II-C spacers") +
  scale_x_discrete(labels=c("Aeromonas\n(+)","Aeromonas\n(-)")) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  ylim(0,11) +
  ggtitle("II-C spacers") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
spacerdata_21_totalspacers_plot_C1

spacerdata_21_totalspacers_plot_C2 <- ggplot(data = dada_21_subset, aes(x=Condition, y=spacers_C2)) +
  geom_boxplot() +
  geom_jitter(aes(colour = morphology), alpha = 0.5, width = 0.2, height = 0) +
  ylab("VI-B spacers") +
  ggtitle("VI-B spacers") +
  scale_x_discrete(labels=c("Aeromonas\n(+)","Aeromonas\n(-)")) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  ylim(0,11) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
spacerdata_21_totalspacers_plot_C2

fig_6_both_spacers_aeromonas <- spacerdata_21_totalspacers_plot + spacerdata_21_totalspacers_plot_C1 + spacerdata_21_totalspacers_plot_C2

ggsave(fig_6_both_spacers_aeromonas, filename = "manuscript_figs/Fig_6_both_spacers_aeromonas.png", width = 8, height = 3)

#MODELS

dada_21$Condition <- relevel(as.factor(dada_21$Condition),ref = "columnare_phage") # Declare columnare with phage (no aeromonas) as reference

#model 1: total number of spacers explained by condition
dada_21_model1 <- glmmTMB(totalspacers ~ Condition + (1|rep) + (1|day), 
                          family = nbinom2,
                          ziformula=~0,
                          data = dada_21_subset)
summary(dada_21_model1)
AIC(dada_21_model1)

#model 2: probability of having aq explained by condition

dada_21_model2 <- glmmTMB(as.numeric(booleanNewSpacers) ~ Condition + (1|rep) + (1|day), 
                          family = binomial,
                          data = dada_21)

summary(dada_21_model2)
AIC(dada_21_model2)
#### Supplementary plots and associated models ####

#Effect of time on OD in follow-up experiment
plot_extractTimeAndTimetoMaxODP <- ggplot(data = subset(dadaOD_long, phage == "yes" & ancestral == "no"), aes(x = extractTime, y = time)) +
  geom_point(size = 1.5, alpha = point_alpha, shape = 21) +
  facet_grid(~sample) +
  ggtitle("With phage") +
  xlab("") + ylab("Time to max OD (min)")

plot_extractTimeAndTimetoMaxODNP <- ggplot(data = subset(dadaOD_long, phage == "no" & ancestral == "no"), aes(x = extractTime, y = time)) +
  geom_point(size = 1.5, alpha = point_alpha, shape = 21) +
  facet_grid(~sample) +
  ggtitle("No phage") +
  xlab("") + ylab("")

plot_extractTimeAndMaxODP <- ggplot(data = subset(dadaOD_long, phage == "yes" & ancestral == "no"), aes(x = extractTime, y = max_od)) +
  geom_point(size = 1.5, alpha = point_alpha, shape = 21) +
  facet_grid(~sample) +
  ggtitle("") +
  xlab("Extraction time") + ylab("Max OD")

plot_extractTimeAndMaxODNP <- ggplot(data = subset(dadaOD_long, phage == "no" & ancestral == "no"), aes(x = extractTime, y = max_od)) +
  geom_point(size = 1.5, alpha = point_alpha, shape = 21) +
  facet_grid(~sample) +
  ggtitle("") +
  xlab("Extraction time") + ylab("")

extractionTimePlots <- (plot_extractTimeAndTimetoMaxODP | plot_extractTimeAndTimetoMaxODNP) / (plot_extractTimeAndMaxODP | plot_extractTimeAndMaxODNP) +
  plot_annotation("Effect of extraction time on growth")
ggsave(extractionTimePlots, filename = "manuscript_figs/SI/extractTimePlots.png", height = 4, width = 8)


titers_21_plot_all <- ggplot(data = subset(data_21_titers_long, titer_type == "titer_total"), 
                             aes(x=as.factor(sampling_day),y=value), group = species) +
  stat_summary(aes(y = value, group = species, color= species), fun=mean, geom="line") +
  geom_point(aes(color = species), alpha = 0.5) +
  ylab("Titer") +
  xlab("Sampling day") +
  scale_y_continuous(trans='log10') +
  facet_grid(~condition_alt) +
  theme_minimal() +
  scale_colour_manual(values=cbPalette)
titers_21_plot_all

ggsave(titers_21_plot_all, filename = "manuscript_figs/SI/SI_fig_4_titers_aeromonasExp_all.png", width = 16, height = 3)


library(lme4)
#Modeling extraction time on ODMAX
summary(lmer(max_od ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "no" & ancestral == "no" & sample == "LW")))
summary(lmer(max_od ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "no" & ancestral == "no" & sample == "LW_M")))
summary(lmer(max_od ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "no" & ancestral == "no" & sample == "Shieh")))
summary(lmer(max_od ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "no" & ancestral == "no" & sample == "Shieh_M")))

summary(lmer(max_od ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "yes" & ancestral == "no" & sample == "LW")))
summary(lmer(max_od ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "yes" & ancestral == "no" & sample == "LW_M")))
summary(lmer(max_od ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "yes" & ancestral == "no" & sample == "Shieh")))
summary(lmer(max_od ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "yes" & ancestral == "no" & sample == "Shieh_M")))
#... and on time-to-ODMAX
summary(lmer(time ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "no" & ancestral == "no" & sample == "LW")))
summary(lmer(time ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "no" & ancestral == "no" & sample == "LW_M")))
summary(lmer(time ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "no" & ancestral == "no" & sample == "Shieh")))
summary(lmer(time ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "no" & ancestral == "no" & sample == "Shieh_M")))

summary(lmer(time ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "yes" & ancestral == "no" & sample == "LW")))
summary(lmer(time ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "yes" & ancestral == "no" & sample == "LW_M")))
summary(lmer(time ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "yes" & ancestral == "no" & sample == "Shieh")))
summary(lmer(time ~ extractTime + (1|rep), data = subset(dadaOD_long, phage == "yes" & ancestral == "no" & sample == "Shieh_M")))


#### Models ####

# Follow-up growth experiment overall (fig 4A)
dadaOD_long$sample <- relevel(dadaOD_long$sample,ref = "B245") # Declare B245 as reference
hist(subset(dadaOD_long, ctrl =="no")$max_od)
treatODphage_nonrelative <- glmmTMB(max_od ~ phage * sample + (1|rep) + (1|id) + (1|time), 
                                    family = gaussian, 
                                    data = subset(dadaOD_long, ctrl =="no"))
summary(treatODphage_nonrelative)
AIC(treatODphage_nonrelative)

treatODphage_nonrelative_phage <- glmmTMB(max_od ~ sample + (1|rep) + (1|id) + (1|time), 
                                          family = gaussian, 
                                          data = subset(dadaOD_long, ctrl =="no" & phage =="yes"))
summary(treatODphage_nonrelative_phage)
AIC(treatODphage_nonrelative_phage)

#the same for time-to-ODmax
hist(subset(dadaOD_long, ctrl =="no")$time)
treatODphage_nonrelative_time <- glmmTMB(time ~ phage * sample + (1|rep), 
                                         family = gaussian, 
                                         data = subset(dadaOD_long, ctrl =="no"))
summary(treatODphage_nonrelative_time)
AIC(treatODphage_nonrelative_time)

#controls (growth experiment, fig 4B)
summary(lm(data = subset(dadaOD_long, ctrl == "yes" | ancestral == "yes"), max_od ~ phage*sample))


### 16-week experiment
#LW and LW+M phage
summary(lm(data = subset(titers, week > 8 & (sample == "LW" | sample == "LW_M")), phage_titer ~ sample))
#LW and LW+M bacteria
summary(lm(data = subset(titers, week > 1 & (sample == "LW" | sample == "LW_M")), bac_tot ~ sample))

#Spacer accumulation in different settings
dadaODSpacersAlllong$sample <- relevel(dadaODSpacersAlllong$sample,ref = "LW") # Declare LW_M as reference

spacersTestLW_LWM_LW <- subset(dadaODSpacersAlllong, ancestral == "no" & ctrl == "no" & phage == "yes" & (sample == "LW" | sample =="LW_M"))
spacersTestLW_LWM_ShiehMucin <- subset(dadaODSpacersAlllong, ancestral == "no" & ctrl == "no" & phage == "yes" & (sample == "LW_M" | sample =="Shieh_M"))

spacersByTreatment_model_LW <- glmmTMB(totalSpacers ~ sample + (1|rep), family=nbinom2, spacersTestLW_LWM_LW)
summary(spacersByTreatment_model)
AIC(spacersByTreatment_model)

spacersByTreatment_model_LW_ShiehM <- glmmTMB(totalSpacers ~ sample + (1|rep), family=nbinom2, spacersTestLW_LWM_ShiehMucin)
summary(spacersByTreatment_model_LW_ShiehM)
AIC(spacersByTreatment_model)

#### Other ####
#table for showing OD-max and time-to-OD-max
od_table <- subset(dadaOD_long, ancestral == "no" & ctrl == "no")
od_table <- select(od_table, sample,max_od,time,id,phage)
od_table <- od_table[order(as.character(od_table$id)),]
rownames(od_table) <- NULL
write.csv(od_table, "od_metrics.csv")

#correlation between OD-max and time-to-OD-max
od_table_NA <- na.omit(od_table)

cor(od_table_NA$max_od, od_table_NA$time)

plot(od_table_NA$max_od, od_table_NA$time)

cor.test(od_table_NA$max_od, od_table_NA$time)

od_table_NA_phage <- subset(od_table_NA, phage == "yes")
cor.test(od_table_NA_phage$max_od, od_table_NA_phage$time)

od_table_NA_nophage <- subset(od_table_NA, phage == "no")
cor.test(od_table_NA_nophage$max_od, od_table_NA_nophage$time)
