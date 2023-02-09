library(mixOmics)
library(popgen.tools)
library(dplyr)

CHB_matrix = read.table('/home/zhulaw/Documents/Honours/R_main/discoal_plsda/CHB_50KB_MAGEE2.csv', sep = ',', header=FALSE)
CHB_matrix <- mutate_all(CHB_matrix, function(x) as.integer(as.character(x)))
CHB_matrix <- as.matrix(CHB_matrix)
is_genome_matrix(CHB_matrix)

pos_lst = read.table('/home/zhulaw/Documents/Honours/R_main/discoal_plsda/CHB_MAGEE2_50KB_pos_lst.csv', sep=',', head=FALSE)
pos_lst = as.matrix(pos_lst)

CHB_sim = sim_obj(cmd = NA, seeds=NA, , segsites=length(CHB_matrix[1, ]), positions=pos_lst, sweep='CHB', select_coeff=NA, genome_matrix = CHB_matrix)
sum_stats_CHB <- sum_stats(CHB_sim, nwins=2, ID=0, snp=1000)

neutral_df = do.call(rbind.data.frame, sum_stats_neutral)
neutral_df_sum = neutral_df[9:28]
neutral_df_sum <- cbind (neutral_df[2], neutral_df_sum)
ssv_df = do.call(rbind.data.frame, sum_stats_ssv)
ssv_df_sum = ssv_df[9:28]
ssv_df_sum <- cbind (ssv_df[2], ssv_df_sum)
sdn_df = do.call(rbind.data.frame, sum_stats_sdn)
sdn_df_sum = sdn_df[9:28]
sdn_df_sum <- cbind (sdn_df[2], sdn_df_sum)

CHB_df = data.frame(sum_stats_CHB)
CHB_df_sum = CHB_df[9:28]
CHB_df_sum <- cbind (CHB_df[2], CHB_df_sum)


combined_df = rbind(ssv_df_sum, neutral_df_sum, sdn_df_sum)
combined_plsda = plsda(combined_df[c(2:7, 12:21)], combined_df$sweep)

test_predict = predict(combined_plsda, CHB_df_sum[c(2:7, 12:21)], method = "centroids.dist")

plotIndiv(combined_plsda, group = combined_df$sweep, style='graphics', pch=21, title='PLS-DA of Summary Statistics from Coalescent Simulations and 1000Genomes Chinese Han Beijing (CHB)', legend = TRUE, size.legend = 0.8, cex = 1, size.title = 0.8)
points(test_predict$variates, pch = 21, bg = "red")
text(test_predict$variates, labels = c('CHB'), adj = c(0, 1))