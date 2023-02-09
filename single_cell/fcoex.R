install.packages('Seurat')
install.packages('dqrng')
BiocManager::install("fcoex")
BiocManager::install("CEMiTool")
install.packages('dplyr')
install.packages('ggplot2')
install.packages('corrr')

library(Seurat)
library(dqrng)
library(fcoex)
library(corrr)
library(dplyr)
library(stringr)
library(schex)

dir.create('outputs')

ss.seurat <- readRDS('Seurat_ss_magee2_20000.rda')

mito.genes <- grep(pattern = "^mt-", x = rownames(x = ss.seurat@assays$RNA@data), value = TRUE)
percent.mito <- colSums(ss.seurat@assays$RNA@data[mito.genes, ])/colSums(ss.seurat@assays$RNA@data)
ss.seurat <- AddMetaData(object = ss.seurat, metadata = percent.mito, col.name = "percent.mito")
hkgenes.found <- which(toupper(rownames(ss.seurat@assays$RNA@data)) %in% hkgenes)
hkgenes.sum <- colSums(ss.seurat@assays$RNA@data[hkgenes.found, ] > 0)
ss.seurat <- AddMetaData(object = ss.seurat, metadata = hkgenes.sum, col.name = "hkgenes.sum")
ss.seurat <- subset(ss.seurat, subset = nFeature_RNA > 800 & nFeature_RNA < 15000 & percent.mito < 0.05 & hkgenes.sum > 10 & nCount_RNA > 5000 )
cat('Total number of cells after filtering: ', ss.seurat@assays$RNA@data@Dim[2], '\n')

ss.seurat <- NormalizeData(ss.seurat, scale.factor = 1e6)

#Female
ss.seurat_female <- subset(ss.seurat, subset = donor_sex_label == 'F')
FM_ss_F <- FindMarkers(ss.seurat_female, ident.1 = "Pos", ident.2 = "Neg", only.pos = TRUE)
write.csv(FM_ss_F, 'outputs/FM_ss_F.csv')


#some correlation functions
correlation_matrix <- function(sample, gene) {
  matrix_mod <- as.matrix(sample@assays$RNA@data)
  
  gene <- as.numeric(matrix_mod[gene,])
  correlations <- apply(matrix_mod,1,function(x){cor(gene,x)})
  correlations <- correlations[order(correlations, decreasing = TRUE)]
  correlations <- correlations[which(is.na(correlations) == FALSE)]
  return(correlations)
}

simple_correlation_network <- function(sample, gene) {
  correlations <- correlation_matrix(sample, gene)
  top_genes <- append(names(correlations[2:10]), 
                      names(tail(correlations, 10)))
  print(top_genes)
  network_genes <- append("Magee2", top_genes)
  for (top_gene in top_genes){
    secondary_correlations <- correlation_matrix(sample, top_gene)
    upper_correlations <- secondary_correlations[
      which(secondary_correlations > 0.5)]
    lower_correlations <- secondary_correlations[
      which(secondary_correlations < -0.25)]
    if (length(lower_correlations) > 4){
      network_genes <- append(network_genes, 
                              names(upper_correlations[1:15]))
      network_genes <- append(network_genes, names(lower_correlations[1:5]))
    }
    else {
      network_genes <- append(network_genes, 
                              names(upper_correlations[1:20]))
    }
  }
  network_genes <- network_genes[which(is.na(network_genes) == FALSE)]
  rm(secondary_correlations, lower_correlations)
  expr_table <- as_tibble(data.frame
                          (matrix(nrow=sample@assays$RNA@data@Dim[2],
                                  ncol=0)))
  # colnames(expr_table) <- network_genes
  for (network_gene in network_genes) {
    gene_idx <- which(rownames(ss.seurat@assays$RNA@data) == network_gene)
    expr_table[, network_gene] <- sample@assays$RNA@data[gene_idx, ]
  }
  res.cor <- correlate(expr_table)
  return(res.cor)
}


print('computing correlation matrices..')
dir.create('outputs/pearson_corr')
correlations <- correlation_matrix(ss.seurat, "Magee2")
write.csv(correlations, 'outputs/SS_full_corr.csv')

#simple correlation network
simple_corr_matrix <- simple_correlation_network(ss.seurat, 'Magee2')
png('outputs/pearson_corr/simple_correlation_network.png', width = 960, height = 540)
network_plot(simple_corr_matrix, min_cor = 0.1, repel = TRUE)
dev.off()

#for each class_label
dir.create('outputs/pearson_corr/class_label')
labels = unique(ss.seurat@meta.data$class_label)
for (label in labels[1:length(labels)-1]) {
  label_subset <- subset(ss.seurat, subset = class_label == label)
  print(label_subset@assays$RNA@data@Dim[[2]])
  label_corr <- correlation_matrix(label_subset, 'Magee2')
  label_filename = paste0('outputs/pearson_corr/class_label/', label, '.csv')
  write.csv(label_corr, file=label_filename, row.names=T)
}

#for each subclass_label
dir.create('outputs/pearson_corr/subclass_label')
labels = unique(ss.seurat@meta.data$subclass_label)
for (label in labels[1:length(labels)-1]) {
  label_subset <- subset(ss.seurat, subset = subclass_label == label)
  print(label_subset@assays$RNA@data@Dim[[2]])
  label_corr <- correlation_matrix(label_subset, 'Magee2')
  label_filename = paste0('outputs/pearson_corr/subclass_label/', label, '.csv')
  write.csv(label_corr, file=label_filename, row.names=T)
}

#for each region_label
dir.create('outputs/pearson_corr/region_label')
labels = unique(ss.seurat@meta.data$region_label)
for (label in labels[1:length(labels)-1]) {
  label_subset <- subset(ss.seurat, subset = region_label == label)
  print(label_subset@assays$RNA@data@Dim[[2]])
  label_corr <- correlation_matrix(label_subset, 'Magee2')
  label_filename = paste0('outputs/pearson_corr/region_label/', label, '.csv')
  write.csv(label_corr, file=label_filename, row.names=T)
}

#for each neighborhood_label
dir.create('outputs/pearson_corr/neighborhood_label')
labels = unique(ss.seurat@meta.data$neighborhood_label)
for (label in labels[1:length(labels)-1]) {
  label_subset <- subset(ss.seurat, subset = neighborhood_label == label)
  print(label_subset@assays$RNA@data@Dim[[2]])
  label_corr <- correlation_matrix(label_subset, 'Magee2')
  label_filename = paste0('outputs/pearson_corr/neighborhood_label/', label, '.csv')
  write.csv(label_corr, file=label_filename, row.names=T)
}


#standard clustering using PCA and tsne reduction
ss.seurat <- FindVariableFeatures(object = ss.seurat, nfeatures = 3000)
ss.seurat <- ScaleData(object = ss.seurat)
ss.seurat <- RunPCA(object = ss.seurat)
ss.seurat <- FindNeighbors(object = ss.seurat)
ss.seurat <- FindClusters(object = ss.seurat)
ss.seurat <- RunTSNE(object = ss.seurat, check_duplicates = FALSE)

print('getting some tsne reduction plots')
png('outputs/Magee2_pca_tsne_plots.png', width = 960, height = 540)
FeaturePlot(object = ss.seurat, features = c("Magee2"), pt.size = 0.5, max.cutoff = 'q95', ncol = 1)
dev.off()

png('outputs/subclass_label_pca_tsne_plots.png', width = 960, height = 540)
DimPlot(ss.seurat, reduction='tsne', group.by = 'subclass_label')
dev.off()

png('outputs/class_label_pca_tsne_plots.png', width = 960, height = 540)
DimPlot(ss.seurat, reduction='tsne', group.by = 'class_label')
dev.off()

png('outputs/neighborhood_label_pca_tsne_plots.png', width = 960, height = 540)
DimPlot(ss.seurat, reduction='tsne', group.by = 'neighborhood_label')
dev.off()



#fcoex workflow
target <- ss.seurat@meta.data$magee2

# Get normalized table from the pre-processing
exprs <- as.data.frame(ss.seurat@assays$RNA@data)

# Create fcoex object
fc <- new_fcoex(data.frame(exprs),target)

rm(exprs, target)

fc <- discretize(fc, number_of_bins = 8)

fc <- find_cbf_modules(fc,n_genes_selected_in_first_step = 200, verbose = FALSE, is_parallel = FALSE)

fc <- get_nets(fc)

network_plots <- show_net(fc)

print('saving network plots')
for (name in names(network_plots)){
  pdf(str_glue('outputs/',name,'_network_plot.pdf'))
  network_plots[[name]]
  dev.off()
}

fc <- recluster(fc)

for (name in names(fc@module_list)){
  mod_name = str_glue('mod_',name)
  ss.seurat@meta.data[[mod_name]] <- idents(fc)[[name]]
}

ss.seurat <- make_hexbin(ss.seurat, nbins = 40, dimension_reduction = "tsne")

print('saving hexplots for each module')
for (name in names(fc@module_list)){
  pdf(str_glue('outputs/',name,'_hexplot.pdf'))
  plot_hexbin_meta(ss.seurat, col=str_glue('mod_',name), action="majority")
  dev.off()
}


#gene module enrichment
fc@module_list$Magee2 <- toupper(fc@module_list$Magee2)

gmt_filename <- system.file("extdata", "pathways.gmt", package = "CEMiTool")

if (gmt_filename == "")
{
  print("You likely need to install CEMiTool")
} else {
  gmt_in <- pathwayPCA::read_gmt(gmt_filename,  description = TRUE)
  
}

fc <- mod_ora(fc, gmt_in)
fc <- plot_ora(fc)

print('saving ora barplots')
for (name in names(fc@module_list)){
  pdf(str_glue('outputs/',name,'_fcoex_ora.pdf'))
  fc@barplot_ora[[name]]
  dev.off()
}


print('complete!')