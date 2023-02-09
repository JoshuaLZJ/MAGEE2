library(dplyr)
library(ggplot2)
library(Seurat)
library(tidyr)
library(tibble)

ss.seurat <- readRDS('Seurat_ss_magee2.rda')
# 
# 
print('getting some barplots')
png('barplot_region_magee2_proportion.png', width = 1240, height = 1395)
par(mar=c(8,4.1,4.1,2.1))
pos_table = table(ss.seurat@meta.data$region_label[which(ss.seurat@meta.data$magee2 == 'Pos')])
full_table = table(ss.seurat@meta.data$region_label)
for (name in names(full_table)) {
  if (!name %in% names(pos_table)) {
    pos_table[name] = 0
  }
}
pos_table = pos_table[order(names(pos_table))]
full_table = subset(full_table, subset = names(full_table) != 'CLA')
pos_table = subset(pos_table, subset = names(pos_table) != 'CLA')
as_tibble(cbind(as.data.frame(( pos_table / full_table) * 100),
                paste0(as.character(pos_table),"/", as.character(full_table)))) %>%
  `colnames<-`(c("Var1", "Freq", "Cell_counts")) %>%
  ggplot(aes(x=reorder(Var1, desc(Freq)), y=Freq)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=Cell_counts), hjust = -0.05, size=12, angle=90) +
  # ggtitle("Plot of Magee2 Cell Proportions and cell counts (on top) for each brain region") +
  xlab("Brain Region") + ylab("Proportion (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=35),
        axis.text.y = element_text(size=35),
        axis.title = element_text(size=40),
        plot.margin = margin(t=100, l=10, unit="pt"),
        axis.line=element_line(linewidth=2)) +
  coord_cartesian(clip = "off")
dev.off()
# 
# png('outputs/barplot_neighborhood_magee2_proportion.png', width = 2480, height = 1395)
# par(mar=c(8,4.1,4.1,2.1))
# pos_table = table(ss.seurat@meta.data$neighborhood_label[which(ss.seurat@meta.data$magee2 == 'Pos')])
# full_table = table(ss.seurat@meta.data$neighborhood_label)
# for (name in names(full_table)) {
#   if (!name %in% names(pos_table)) {
#     pos_table[name] = 0
#   }
# }
# pos_table = pos_table[order(names(pos_table))]
# as_tibble(cbind(as.data.frame(( pos_table / full_table) * 100), paste0(as.character(pos_table),"/", as.character(full_table)))) %>% `colnames<-`(c("Var1", "Freq", "Cell_counts")) %>% ggplot(aes(x=Var1, y=Freq)) + geom_bar(stat="identity") + geom_text(aes(label=Cell_counts), vjust = -0.5) + ggtitle("Plot of Magee2 Cell Proportions and cell counts (on top) for each neighbourhood of cells") + xlab("Neighbourhood") + ylab("Proportion (%)") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10), title = element_text(size = 15))
# dev.off()

png('barplot_subclass_magee2_proportion.png', width = 2480, height = 1395)
par(mar=c(8,4.1,4.1,2.1))
pos_table = table(ss.seurat@meta.data$subclass_label[which(ss.seurat@meta.data$magee2 == 'Pos')])
full_table = table(ss.seurat@meta.data$subclass_label)
for (name in names(full_table)) {
  if (!name %in% names(pos_table)) {
    pos_table[name] = 0
  }
}
pos_table = pos_table[order(names(pos_table))]
as_tibble(cbind(as.data.frame(( pos_table / full_table) * 100), 
                paste0(as.character(pos_table),"/", as.character(full_table)))) %>% 
  `colnames<-`(c("Var1", "Freq", "Cell_counts")) %>% 
  ggplot(aes(x=reorder(Var1, desc(Freq)), y=Freq)) + geom_bar(stat="identity") + 
  geom_text(aes(label=Cell_counts), hjust = -0.1, size=12.5, angle=90) + 
  # ggtitle("Plot of Magee2 Cell Proportions and cell counts (on top) for each subclass of cells") + 
  xlab("Subclass") + ylab("Proportion (%)") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=35), 
        axis.text.y = element_text(size=35),
        axis.title = element_text(size=40),
        plot.margin = margin(t=100, l=10, unit="pt"),
        axis.line=element_line(linewidth=2)) +
  coord_cartesian(clip = "off")
dev.off()

png('barplot_class_magee2_proportion.png', width = 1240, height = 1395)
par(mar=c(8,4.1,4.1,2.1))
pos_table = table(ss.seurat@meta.data$class_label[which(ss.seurat@meta.data$magee2 == 'Pos')])
full_table = table(ss.seurat@meta.data$class_label)
for (name in names(full_table)) {
  if (!name %in% names(pos_table)) {
    pos_table[name] = 0
  }
}
pos_table = pos_table[order(names(pos_table))]
as_tibble(cbind(as.data.frame(( pos_table / full_table) * 100),
                paste0(as.character(pos_table),"/", as.character(full_table)))) %>%
  `colnames<-`(c("Var1", "Freq", "Cell_counts")) %>%
  ggplot(aes(x=reorder(Var1, desc(Freq)), y=Freq)) + geom_bar(stat="identity", width=0.5) +
  geom_text(aes(label=Cell_counts), vjust = -0.5, size=15) +
  # ggtitle("Plot of Magee2 Cell Proportions and cell counts (on top) for each class of cells") +
  xlab("Class") + ylab("Proportion (%)") +
  theme_classic() +
  theme(axis.text.x=element_text(size=35) , axis.text.y=element_text(size=35),
        axis.title = element_text(size=40),
        axis.line=element_line(linewidth=2))
dev.off()


#looking at sexual differences in region

male_region <- as_tibble_col(ss.seurat$region_label[which(ss.seurat$donor_sex_label == 'M')], column_name='region') %>% 
  add_column(magee2 = ss.seurat$magee2[which(ss.seurat$donor_sex_label == 'M')])
male_region_prop <- male_region %>% group_by(region, .drop = FALSE) %>% 
  count(magee2) %>% mutate(freq = (n/sum(n))*100)
male_region_prop <- male_region_prop[which(male_region_prop$region != 'CLA'), ]
male_region_copy <- male_region_prop
male_region_prop <- male_region_prop[which(male_region_prop$magee2 == 'Pos'), ]
male_region_total_cell_counts <- male_region_copy %>% group_by(region) %>% 
  summarise(sum = sum(n))
male_region_prop$total_cell_counts <- male_region_total_cell_counts$sum

female_region <- as_tibble_col(ss.seurat$region_label[which(ss.seurat$donor_sex_label == 'F')], column_name='region') %>% 
  add_column(magee2 = ss.seurat$magee2[which(ss.seurat$donor_sex_label == 'F')])
female_region_prop <- female_region %>% group_by(region, .drop = FALSE) %>% 
  count(magee2) %>% mutate(freq = (n/sum(n))*100)
female_region_prop <- female_region_prop[which(female_region_prop$region != 'CLA'), ]
female_region_copy <- female_region_prop
female_region_prop <- female_region_prop[which(female_region_prop$magee2 == 'Pos'), ]
female_region_total_cell_counts <- female_region_copy %>% group_by(region) %>% 
  summarise(sum = sum(n))
female_region_prop$total_cell_counts <- female_region_total_cell_counts$sum

combined_region_diff <- as_tibble_col(male_region_prop$region, column_name='region') %>% 
  add_column(difference = male_region_prop$freq - female_region_prop$freq)
combined_region_diff <- combined_region_diff %>% arrange(desc(difference))
combined_region_diff$sex <- NA
combined_region_diff[which(combined_region_diff$difference > 0), ]$sex <- 'M'
combined_region_diff[which(combined_region_diff$difference < 0), ]$sex <- 'F'

combined_region_diff$cell_count <- NA
combined_region_diff$male_magee2_count <- NA
combined_region_diff$female_magee2_count <- NA
for (reg in combined_region_diff[which(combined_region_diff$difference > 0), ]$region) {
  combined_region_diff[which(combined_region_diff$region == reg), ]$cell_count <- paste0(
    male_region_prop[which(male_region_prop$region == reg), ]$n, 
    '/', male_region_prop[which(male_region_prop$region == reg), ]
    $total_cell_counts)
}

for (reg in combined_region_diff[which(combined_region_diff$difference < 0), ]$region) {
  combined_region_diff[which(combined_region_diff$region == reg), ]$cell_count <- paste0(
    female_region_prop[which(female_region_prop$region == reg), ]$n, 
    '/', female_region_prop[which(female_region_prop$region == reg), ]
    $total_cell_counts)
}

for (reg in combined_region_diff$region) {
  combined_region_diff[which(combined_region_diff$region == reg), ]$male_magee2_count <- male_region_prop[which(male_region_prop$region == reg), ]$n
  combined_region_diff[which(combined_region_diff$region == reg), ]$female_magee2_count <- female_region_prop[which(female_region_prop$region == reg), ]$n
}


combined_region_diff <- combined_region_diff[which(combined_region_diff$male_magee2_count > 100 & 
                                                     combined_region_diff$female_magee2_count > 100), ]

combined_region_diff_top <- combined_region_diff[which(combined_region_diff$difference > 10 
                                                       | combined_region_diff$difference < -10), ]

#plotting sexual differences in region

png('sexual_diff_region.png', width=2480, height=1395)
ggplot(combined_region_diff_top,aes(x=reorder(region, difference) ,
                                 y=difference,fill=sex))+
  geom_bar(stat="identity", width=0.8)+
  coord_flip()+
  ylab('Difference in percentage (%) of Magee2-expressing cells')+
  xlab('Brain Region')+
  # geom_text(aes(label=cell_count), size=15, 
  #           position = position_dodge2(width=0.5, padding=0.75))+
  ggtitle("Plot of Sexual Differences in Magee2 Cell Proportions and 
          cell counts for each region of cells") + 
  scale_fill_discrete(name = "Gender") +
  theme(aspect.ratio = 1/2.5, 
        title = element_text(size = 50),
        plot.title = element_text(size = 50, hjust=0.5),
        legend.text=element_text(size=50),
        axis.text=element_text(size=40),
        axis.line=element_line(linewidth=2),
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(2, "cm"))
dev.off()

#skipping subclass because no subclass in either male or female mice have magee2-expressing cells > 100 and percentage diff > 10%

# #looking at sexual differences in subclass
# 
# male_subclass <- as_tibble_col(ss.seurat$subclass_label[which(ss.seurat$donor_sex_label == 'M')], column_name='subclass') %>% 
#   add_column(magee2 = ss.seurat$magee2[which(ss.seurat$donor_sex_label == 'M')])
# male_subclass_prop <- male_subclass %>% group_by(subclass, .drop = FALSE) %>% 
#   count(magee2) %>% mutate(freq = (n/sum(n))*100)
# male_subclass_prop <- male_subclass_prop[which(male_subclass_prop$subclass != 'NA'), ]
# male_subclass_copy <- male_subclass_prop
# male_subclass_prop <- male_subclass_prop[which(male_subclass_prop$magee2 == 'Pos'), ]
# male_subclass_total_cell_counts <- male_subclass_copy %>% group_by(subclass) %>% 
#   summarise(sum = sum(n))
# male_subclass_prop$total_cell_counts <- male_subclass_total_cell_counts$sum
# 
# female_subclass <- as_tibble_col(ss.seurat$subclass_label[which(ss.seurat$donor_sex_label == 'F')], column_name='subclass') %>% 
#   add_column(magee2 = ss.seurat$magee2[which(ss.seurat$donor_sex_label == 'F')])
# female_subclass_prop <- female_subclass %>% group_by(subclass, .drop = FALSE) %>% 
#   count(magee2) %>% mutate(freq = (n/sum(n))*100)
# female_subclass_prop <- female_subclass_prop[which(female_subclass_prop$subclass != 'NA'), ]
# female_subclass_copy <- female_subclass_prop
# female_subclass_prop <- female_subclass_prop[which(female_subclass_prop$magee2 == 'Pos'), ]
# female_subclass_total_cell_counts <- female_subclass_copy %>% group_by(subclass) %>% 
#   summarise(sum = sum(n))
# female_subclass_prop$total_cell_counts <- female_subclass_total_cell_counts$sum
# 
# combined_subclass_diff <- as_tibble_col(male_subclass_prop$subclass, column_name='subclass') %>% 
#   add_column(difference = male_subclass_prop$freq - female_subclass_prop$freq)
# combined_subclass_diff <- combined_subclass_diff %>% arrange(desc(difference))
# combined_subclass_diff$sex <- NA
# combined_subclass_diff[which(combined_subclass_diff$difference > 0), ]$sex <- 'M'
# combined_subclass_diff[which(combined_subclass_diff$difference < 0), ]$sex <- 'F'
# 
# combined_subclass_diff$cell_count <- NA
# combined_subclass_diff$male_magee2_count <- NA
# combined_subclass_diff$female_magee2_count <- NA
# for (reg in combined_subclass_diff[which(combined_subclass_diff$difference > 0), ]$subclass) {
#   combined_subclass_diff[which(combined_subclass_diff$subclass == reg), ]$cell_count <- paste0(
#     male_subclass_prop[which(male_subclass_prop$subclass == reg), ]$n, 
#     '/', male_subclass_prop[which(male_subclass_prop$subclass == reg), ]
#     $total_cell_counts)
# }
# 
# for (reg in combined_subclass_diff[which(combined_subclass_diff$difference < 0), ]$subclass) {
#   combined_subclass_diff[which(combined_subclass_diff$subclass == reg), ]$cell_count <- paste0(
#     female_subclass_prop[which(female_subclass_prop$subclass == reg), ]$n, 
#     '/', female_subclass_prop[which(female_subclass_prop$subclass == reg), ]
#     $total_cell_counts)
# }
# 
# for (reg in combined_subclass_diff$subclass) {
#   combined_subclass_diff[which(combined_subclass_diff$subclass == reg), ]$male_magee2_count <- male_subclass_prop[which(male_subclass_prop$subclass == reg), ]$n
#   combined_subclass_diff[which(combined_subclass_diff$subclass == reg), ]$female_magee2_count <- female_subclass_prop[which(female_subclass_prop$subclass == reg), ]$n
# }
# 
# 
# combined_subclass_diff <- combined_subclass_diff[which(combined_subclass_diff$male_magee2_count > 100 & 
#                                                          combined_subclass_diff$female_magee2_count > 100), ]
# 
# combined_subclass_diff_top <- combined_subclass_diff[which(combined_subclass_diff$difference > 10 
#                                                            | combined_subclass_diff$difference < -10), ]
# 
# #plotting sexual differences in subclass
# 
# ggplot(combined_subclass_diff_top,aes(x=reorder(subclass, difference) ,
#                                       y=difference,fill=sex))+geom_bar(stat="identity", width=0.8)+
#   coord_flip()+ylab('Difference in percentage (%) of Magee2-expressing cells')+
#   xlab('Brain subclass')+geom_text(aes(label=cell_count), 
#                                    position = position_dodge2(width=.75, padding=0.5))+
#   ggtitle("Plot of Sexual Differences in Magee2 Cell Proportions and 
#           cell counts for each subclass of cells") + 
#   theme(text=element_text(size=15), aspect.ratio = 1/2.5, 
#         title = element_text(size = 30))




