library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

gencog = read.csv('final_gen_cognitive.csv')
gencog$derived = replace(gencog$derived, gencog$derived=="YES", "AA, CA, A")
gencog$derived = replace(gencog$derived, gencog$derived=="NO", "CC, C")
colnames(gencog)[9] = "Genotypes"


#for all
#2480 x 1395 jpeg

stat.test.sdr <- gencog %>%
  t_test(SDR.GA ~ Genotypes) %>%
  add_significance() %>%
  add_xy_position(x="Genotypes")

stat.test.ldr <- gencog %>%
  t_test(LDR.GA ~ Genotypes) %>%
  add_significance() %>%
  add_xy_position(x="Genotypes")

stat.test.london <- gencog %>%
  t_test(London_GeneralAccuracy ~ Genotypes) %>%
  add_significance() %>%
  add_xy_position(x="Genotypes")

##standalone sdr_all
# 
# sdr_all = ggboxplot(gencog, x="Genotypes", y="SDR.GA", fill="Genotypes",width=0.5, outlier.shape=NA) +
#   geom_jitter() + stat_pvalue_manual(stat.test.sdr,
#                                      label = "Welch's t-test one-tailed p-value = {p/2}", size = 15) +
#   xlab('rs1343879 Genotypes') +
#   ylab("Short-Term Word Recall Score") +
#   scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, by=0.2)) +
#   # ggtitle("Short-Term Word Recall Scores of Individuals with and without the Nonsense Allele (A)") +
#   theme(legend.position="none", axis.title=element_text(size=30),
#         axis.text.x=element_text(size=30), axis.text.y=element_text(size=30), title=element_text(size=25))
# 
# png('gencog_sdr.png', width=2480, height=1395)
# sdr_all
# dev.off()

sdr_all = ggboxplot(gencog, x="Genotypes", y="SDR.GA", fill="Genotypes",width=0.5, outlier.shape=NA) + 
  geom_jitter() + stat_pvalue_manual(stat.test.sdr, 
                                     label = "One-tailed p-value = {p/2}", size = 10) + 
  xlab('rs1343879 Genotypes') +
  ylab("Short-Term Word Recall Score") + 
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, by=0.2)) +
  # ggtitle("Short-Term Word Recall Scores of Individuals with and without the Nonsense Allele (A)") + 
  theme(legend.position="none", axis.text=element_text(size=25), title=element_text(size=25))

ldr_all = ggboxplot(gencog, x="Genotypes", y="LDR.GA", fill="Genotypes",width=0.5, outlier.shape=NA) + 
  geom_jitter() + stat_pvalue_manual(stat.test.ldr, 
                                     label = "One-tailed p-value = {p/2}", size = 10) + 
  xlab('rs1343879 Genotypes') +
  ylab("Long-Term Word Recall Score") + 
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, by=0.2)) +
  # ggtitle("Long-Term Word Recall Scores of Individuals with and without the Nonsense Allele (A)") + 
  theme(legend.position="none", axis.text=element_text(size=25), title=element_text(size=25))

london_all = ggboxplot(gencog, x="Genotypes", y="London_GeneralAccuracy", 
                       fill="Genotypes",width=0.5, outlier.shape=NA) + 
  geom_jitter() + stat_pvalue_manual(stat.test.london, 
                                     label = "One-tailed p-value = {p/2}", size = 10) + 
  xlab('rs1343879 Genotypes') +
  ylab("Tower of London Accuracy Score") + 
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, by=0.2)) +
  # ggtitle("Tower of London Accuracy Scores of Individuals with and without the Nonsense Allele (A)") + 
  theme(legend.position="none", axis.text=element_text(size=25), title=element_text(size=25))

figure_all <- ggarrange(sdr_all, ldr_all, london_all, ncol=3, nrow=1)
figure_all

male_gencog = gencog[which((gencog['Gen']== "C") | (gencog['Gen'] == "A")), ]

stat.test.sdr.m <- male_gencog %>%
  t_test(SDR.GA ~ Genotypes) %>%
  add_significance() %>%
  add_xy_position(x="Genotypes")

stat.test.ldr.m <- male_gencog %>%
  t_test(LDR.GA ~ Genotypes) %>%
  add_significance() %>%
  add_xy_position(x="Genotypes")

stat.test.london.m <- male_gencog %>%
  t_test(London_GeneralAccuracy ~ Genotypes) %>%
  add_significance() %>%
  add_xy_position(x="Genotypes")


sdr_male = ggboxplot(male_gencog, x="Genotypes", y="SDR.GA", fill="Genotypes",width=0.5, outlier.shape=NA) + 
  geom_jitter() + stat_pvalue_manual(stat.test.sdr.m, 
                                     label = "One-tailed p-value = {p/2}", size = 10) + 
  xlab('rs1343879 Genotypes') +
  ylab("Short-Term Word Recall Score") + 
  # ggtitle("Short-Term Word Recall Scores of Male Individuals with and without the Nonsense Allele (A)") + 
  theme(legend.position="none", axis.text=element_text(size=25), title=element_text(size=25))

ldr_male = ggboxplot(male_gencog, x="Genotypes", y="LDR.GA", fill="Genotypes",width=0.5, outlier.shape=NA) + 
  geom_jitter() + stat_pvalue_manual(stat.test.ldr.m, 
                                     label = "One-tailed p-value = {p/2}", size = 10) + 
  xlab('rs1343879 Genotypes') +
  ylab("Long-Term Word Recall Score") + 
  # ggtitle("Long-Term Word Recall Scores of Male Individuals with and without the Nonsense Allele (A)") + 
  theme(legend.position="none", axis.text=element_text(size=25), title=element_text(size=25))

london_male = ggboxplot(male_gencog, x="Genotypes", y="London_GeneralAccuracy", 
                       fill="Genotypes",width=0.5, outlier.shape=NA) + 
  geom_jitter() + stat_pvalue_manual(stat.test.london.m, 
                                     label = "One-tailed p-value = {p/2}", size = 10) + 
  xlab('rs1343879 Genotypes') +
  ylab("Tower of London Accuracy Score") + 
  # ggtitle("Tower of London Accuracy Scores of Male Individuals with and without the Nonsense Allele (A)") + 
  theme(legend.position="none", axis.text=element_text(size=25), title=element_text(size=25))

figure_male <- ggarrange(sdr_male, ldr_male, london_male, ncol=3, nrow=1)
figure_male


female_gencog = gencog[which((gencog['Gen']== "CC") | (gencog['Gen'] == "AA") | (gencog['Gen'] == "CA")), ]

stat.test.sdr.f <- female_gencog %>%
  t_test(SDR.GA ~ Genotypes) %>%
  add_significance() %>%
  add_xy_position(x="Genotypes")

stat.test.ldr.f <- female_gencog %>%
  t_test(LDR.GA ~ Genotypes) %>%
  add_significance() %>%
  add_xy_position(x="Genotypes")

stat.test.london.f <- female_gencog %>%
  t_test(London_GeneralAccuracy ~ Genotypes) %>%
  add_significance() %>%
  add_xy_position(x="Genotypes")


sdr_female = ggboxplot(female_gencog, x="Genotypes", y="SDR.GA", fill="Genotypes",width=0.5, outlier.shape=NA) + 
  geom_jitter() + stat_pvalue_manual(stat.test.sdr.f, 
                                     label = "One-tailed p-value = {p/2}", size = 10) + 
  xlab('rs1343879 Genotypes') +
  ylab("Short-Term Word Recall Score") + 
  # ggtitle("Short-Term Word Recall Scores of Female Individuals with and without the Nonsense Allele (A)") + 
  theme(legend.position="none", axis.text=element_text(size=25), title=element_text(size=25))

ldr_female = ggboxplot(female_gencog, x="Genotypes", y="LDR.GA", fill="Genotypes",width=0.5, outlier.shape=NA) + 
  geom_jitter() + stat_pvalue_manual(stat.test.ldr.f, 
                                     label = "One-tailed p-value = {p/2}", size = 10) + 
  xlab('rs1343879 Genotypes') +
  ylab("Long-Term Word Recall Score") + 
  # ggtitle("Long-Term Word Recall Scores of Female Individuals with and without the Nonsense Allele (A)") + 
  theme(legend.position="none", axis.text=element_text(size=25), title=element_text(size=25))

london_female = ggboxplot(female_gencog, x="Genotypes", y="London_GeneralAccuracy", 
                        fill="Genotypes",width=0.5, outlier.shape=NA) + 
  geom_jitter() + stat_pvalue_manual(stat.test.london.f, 
                                     label = "One-tailed p-value = {p/2}", size = 10) + 
  xlab('rs1343879 Genotypes') +
  ylab("Tower of London Accuracy Score") + 
  # ggtitle("Tower of London Accuracy Scores of Female Individuals with and without the Nonsense Allele (A)") + 
  theme(legend.position="none", axis.text=element_text(size=25), title=element_text(size=25))

figure_female <- ggarrange(sdr_female, ldr_female, london_female, ncol=3, nrow=1)
figure_female