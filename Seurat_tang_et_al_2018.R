# ==============================================================================
# Author : Paul Gueguen
# Email : paul.gueguen@curie.fr
# Date : 12/12/2017  
# Modified : 22/05/2017
# ==============================================================================
library(Seurat)
library(gplots)
library(dplyr)
library(Matrix)
library(clustree)
library(MAST)

### SIGNATURES
sign_moDC <- read.csv("~/Bureau/Projet/Mono_singlecell/svn/analyse/script/sign_mo-DC.txt", h = F)[,1]
sign_moMAC <- read.csv("~/Bureau/Projet/Mono_singlecell/svn/analyse/script/sign_mo-MAC.txt", h = F)[,1]
sign_activated_DC <- toupper(read.csv("~/Bureau/Signatures/Signatures Papier Tsing Lee/Activated_DC.txt", h = F)[,1])
sign_blood_cDC2 <- toupper(read.csv("~/Bureau/Signatures/Signatures Papier Tsing Lee/blood_cDC2.txt", h = F)[,1])
sign_blood_mono <- toupper(read.csv("~/Bureau/Signatures/Signatures Papier Tsing Lee/blood_mono.txt", h = F)[,1])
sign_blood_DC1 <- toupper(read.csv("~/Bureau/Signatures/Signatures Papier Tsing Lee/DC1.txt", h = F)[,1])
sign_cDC <- toupper(read.csv("~/Bureau/Signatures/Signatures Papier Tsing Lee/genes cDC signature.txt", h = F)[,1])
sign_DC2 <- toupper(read.csv("~/Bureau/Signatures/Signatures Papier Tsing Lee/DC2.txt", h = F)[,1])
sign_in_vitro_moDC <- toupper(read.csv("~/Bureau/Signatures/Signatures Papier Tsing Lee/in_vitro_moDC.txt", h = F)[,1])

### INITIATION
mono_mac_dc_a_dc_am.data <- Read10X("~/Bureau/Projet/Mono_singlecell/svn/project/script/Tsing_Lee/aggregate_monocytes_macro_dc_dc_amygdale/outs/filtered_gene_bc_matrices_mex/GRCh38/")
mono_mac_dc_a_dc_am <- CreateSeuratObject(raw.data = mono_mac_dc_a_dc_am.data, min.cells = 3, min.genes = 200, project = "10X_mono_mac_dc_a_dc_am")

### QC and selecting cells for further analysis
dim(mono_mac_dc_a_dc_am@data)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = mono_mac_dc_a_dc_am@data), value = TRUE)
ribo.genes <- grep(pattern = "^RPL|^RPS", x = rownames(x = mono_mac_dc_a_dc_am@data), value = TRUE)
percent.mito <- Matrix::colSums(mono_mac_dc_a_dc_am@raw.data[mito.genes, ])/Matrix::colSums(mono_mac_dc_a_dc_am@raw.data)
percent.ribo <- Matrix::colSums(mono_mac_dc_a_dc_am@raw.data[ribo.genes, ])/Matrix::colSums(mono_mac_dc_a_dc_am@raw.data)
mono_mac_dc_a_dc_am <- AddMetaData(mono_mac_dc_a_dc_am, percent.mito, "percent.mito")
mono_mac_dc_a_dc_am <- AddMetaData(mono_mac_dc_a_dc_am, percent.ribo, "percent.ribo")
VlnPlot(mono_mac_dc_a_dc_am, features.plot = c("nUMI","nGene","percent.mito","percent.ribo"), nCol = 4)
mono_mac_dc_a_dc_am <- FilterCells(object = mono_mac_dc_a_dc_am, subset.names = c("nGene", "percent.mito", "percent.ribo"), low.thresholds = c(400, -Inf, -Inf), high.thresholds = c(5000, 0.06, 0.5))
dim(mono_mac_dc_a_dc_am@data)

### Orig.ident Grepping
aggr_d1 <-grep ("*-1$",names(mono_mac_dc_a_dc_am@ident), value = T)
aggr_d2 <-grep ("*-2",names(mono_mac_dc_a_dc_am@ident), value = T)
aggr_d3 <-grep ("*-3",names(mono_mac_dc_a_dc_am@ident), value = T)
aggr_d4 <-grep ("*-4",names(mono_mac_dc_a_dc_am@ident), value = T)
aggr_d5 <-grep ("*-5",names(mono_mac_dc_a_dc_am@ident), value = T)

### Setup orig.sample
orig.sample<-numeric()
for (i in 1:length(mono_mac_dc_a_dc_am@cell.names)){
  if(mono_mac_dc_a_dc_am@cell.names[i]%in%aggr_d1){
    orig.sample[i] <- "Monocyte 1"
  }
  if(mono_mac_dc_a_dc_am@cell.names[i]%in%aggr_d2){
    orig.sample[i] <- "Monocyte 2"
  }
  if(mono_mac_dc_a_dc_am@cell.names[i]%in%aggr_d3){
    orig.sample[i] <- "Ascites Macrophages"
  }
  if(mono_mac_dc_a_dc_am@cell.names[i]%in%aggr_d4){
    orig.sample[i] <- "Ascites DCs"
  }
  if(mono_mac_dc_a_dc_am@cell.names[i]%in%aggr_d5){
    orig.sample[i] <- "Tonsil DCs"
  }
}
names(orig.sample) <- WhichCells(mono_mac_dc_a_dc_am)
mono_mac_dc_a_dc_am <- AddMetaData(mono_mac_dc_a_dc_am, orig.sample, "orig.sample")

### Normalizing the data
mono_mac_dc_a_dc_am <- NormalizeData(object = mono_mac_dc_a_dc_am, normalization.method = "LogNormalize", scale.factor = 1e4)

### Detection of variable genes across the single cells
mono_mac_dc_a_dc_am <- FindVariableGenes(object = mono_mac_dc_a_dc_am, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0, y.cutoff = 0.5, do.contour = F)
length(x = mono_mac_dc_a_dc_am@var.genes)

### Scaling the data and removing unwanted sources of variation
mono_mac_dc_a_dc_am <- ScaleData(object = mono_mac_dc_a_dc_am, vars.to.regress = c("nUMI", "percent.mito"))

### Perform linear dimensional reduction
mono_mac_dc_a_dc_am <- RunPCA(object = mono_mac_dc_a_dc_am, pc.genes = mono_mac_dc_a_dc_am@var.genes)
dim(mono_mac_dc_a_dc_am@data)

### Determine statistically significant principal components
PCElbowPlot(mono_mac_dc_a_dc_am, num.pc = 20)

### Cluster the cells
mono_mac_dc_a_dc_am <- FindClusters(mono_mac_dc_a_dc_am, dims.use = 1:19, save.SNN = T, resolution = 0.8)
mono_mac_dc_a_dc_am_sub <- SubsetData(mono_mac_dc_a_dc_am, ident.use= c(0:12)) 

### Run Non-linear dimensional reduction (tSNE) (fig 1)
mono_mac_dc_a_dc_am <- RunTSNE(mono_mac_dc_a_dc_am, dims.use = 1:19, do.fast = T, seed.use = 40)
TSNEPlot(mono_mac_dc_a_dc_am)

p1 <- DimPlot(mono_mac_dc_a_dc_am_sub, reduction.use = "tsne", do.label = T, pt.size=1.5, do.return = T, 
              cols.use = c("red4","purple","red2","deepskyblue","purple3","green","orchid4","darkturquoise","indianred3","deepskyblue3","orchid3","magenta4","mediumpurple4"))
p2 <- DimPlot(mono_mac_dc_a_dc_am_sub, reduction.use = "tsne", do.label = F, pt.size=1.5, group.by = "orig.sample", do.return = T,
              cols.use = c("deepskyblue","red2","green","green4","purple3"))
plot_grid(p1,p2)

### Find Differentially Expressed Genes - find all markers of all clusters (fig 1)
mono_mac_dc_a_dc_am.markers_sub <- FindAllMarkers(object = mono_mac_dc_a_dc_am_sub, only.pos = TRUE, test.use = "MAST")
mono_mac_dc_a_dc_am.markers %>% group_by(cluster) %>% top_n(20, avg_logFC) -> top20
top20_ordered <- top20[which(top20$cluster==5),]
top20_ordered <- rbind(top20_ordered, top20[which(top20$cluster==0),])
top20_ordered <- rbind(top20_ordered, top20[which(top20$cluster==2),])
top20_ordered <- rbind(top20_ordered, top20[which(top20$cluster==8),])
top20_ordered <- rbind(top20_ordered, top20[which(top20$cluster==3),])
top20_ordered <- rbind(top20_ordered, top20[which(top20$cluster==7),])
top20_ordered <- rbind(top20_ordered, top20[which(top20$cluster==9),])
top20_ordered <- rbind(top20_ordered, top20[which(top20$cluster==12),])
top20_ordered <- rbind(top20_ordered, top20[which(top20$cluster==1),])
top20_ordered <- rbind(top20_ordered, top20[which(top20$cluster==10),])
top20_ordered <- rbind(top20_ordered, top20[which(top20$cluster==6),])
top20_ordered <- rbind(top20_ordered, top20[which(top20$cluster==4),])
top20_ordered <- rbind(top20_ordered, top20[which(top20$cluster==11),])

pdf("~/Bureau/heatmap_mono_mac_dc_a_dc_am_sub_good_order_top20.pdf", height = 6, width = 8)
DoHeatmap(object = mono_mac_dc_a_dc_am_sub, genes.use = top20_ordered$gene, slim.col.label = TRUE, remove.key = TRUE, cex.row = 2, group.cex = 5, group.order = c(5,0,2,8,3,7,9,12,1,10,6,4,11))
dev.off()
 
### Computing the DEG between activated mo_DC and activated DC (cluster 9 and 12)
mono_mac_dc_a_dc_am_act <- SubsetData(mono_mac_dc_a_dc_am_sub, ident.use= c(9,12))
mono_mac_dc_a_dc_am_act.markers <- FindAllMarkers(object = mono_mac_dc_a_dc_am_act, only.pos = T, test.use = "MAST")
mono_mac_dc_a_dc_am_act.markers %>% group_by(cluster) %>% top_n(100, avg_logFC) -> top100
pdf("heatmap_mo_dc_activated_dc_activated.pdf", height = 8, width = 6)
DoHeatmap(object = mono_mac_dc_a_dc_am_sub, genes.use = top100$gene, slim.col.label = TRUE, remove.key = F, cex.row = 3, group.cex = 5, disp.min = -5, disp.max = 5, rotate.key = T)
dev.off()

### Signature scores (replace names with the desired signatures) (fig 2)
mono_mac_dc_a_dc_am_sub <- AddModuleScore(mono_mac_dc_a_dc_am_sub, genes.list = list(sign_moDC), enrich.name = "sign_moDC")
pdf("~/Bureau/sign_moDC_red.pdf", width = 9)
FeaturePlot(mono_mac_dc_a_dc_am_sub, features.plot = "sign_moDC1",  cols.use = c("red","grey"), min.cutoff = "-Inf", max.cutoff = "0", no.axes = F, no.legend = F)
dev.off()

### Violinplots top5 per cluster (fig S1)
mono_mac_dc_a_dc_am.markers %>% group_by(cluster) %>% top_n(5, avg_logFC) -> top5
top5_ordered <- top5[which(top5$cluster==5),]
top5_ordered <- rbind(top5_ordered, top5[which(top5$cluster==0),])
top5_ordered <- rbind(top5_ordered, top5[which(top5$cluster==2),])
top5_ordered <- rbind(top5_ordered, top5[which(top5$cluster==8),])
top5_ordered <- rbind(top5_ordered, top5[which(top5$cluster==3),])
top5_ordered <- rbind(top5_ordered, top5[which(top5$cluster==7),])
top5_ordered <- rbind(top5_ordered, top5[which(top5$cluster==9),])
top5_ordered <- rbind(top5_ordered, top5[which(top5$cluster==12),])
top5_ordered <- rbind(top5_ordered, top5[which(top5$cluster==1),])
top5_ordered <- rbind(top5_ordered, top5[which(top5$cluster==10),])
top5_ordered <- rbind(top5_ordered, top5[which(top5$cluster==6),])
top5_ordered <- rbind(top5_ordered, top5[which(top5$cluster==4),])
top5_ordered <- rbind(top5_ordered, top5[which(top5$cluster==11),])
top5_ordered <- top5_ordered[1:65,]
pdf("~/Bureau/Vlnplots_top5_cluster_Tang_et_al.pdf")
for(i in 1:length(top5$gene)){
  print(VlnPlot2(mono_mac_dc_a_dc_am_sub, top5_ordered$gene[i], group.by = "res.0.8", c.levels = c(5,0,2,8,3,7,9,12,1,10,6,4,11), cols.use = type_cols))
}
dev.off()

VlnPlot2(mono_mac_dc_a_dc_am_sub, top5$gene, group.by = "res.0.8", c.levels = c(5,0,2,8,3,7,9,12,1,10,6,4,11), cols.use = type_cols, nCol = 5, single.legend = T)

### Violin plots for DC / MAC genes (Fig S2)
tonsil_dc_genes <- c("YBX1","COTL1","RELB","FAM60A","IER2","TNFAIP2","SPI1","PTP4A2")
ascites_dc_genes <- c("TYROBP","CST7","TNFSF13B","NMT1","CD1E","MNDA")
pdf("~/Bureau/Vlnplots_ascites_dc_genes_Tang_et_al.pdf")
for(i in 1:length(mac_dc_genes)){
  print(VlnPlot2(mono_mac_dc_a_dc_am_sub, ascites_dc_genes[i], group.by = "res.0.8", c.levels = c(5,0,2,8,3,7,9,12,1,10,6,4,11), cols.use = type_cols, point.size.use = 0.5))
}
dev.off()

### Reordering violinplots (you must compute the function to be able to reorder the clusters using c.levels)
VlnPlot2 <- function (object, features.plot, ident.include = NULL, c.levels = NULL, nCol = NULL, 
                      do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16, 
                      size.y.use = 16, size.title.use = 20, adjust.use = 1, point.size.use = 1, 
                      cols.use = NULL, group.by = NULL, y.log = FALSE, x.lab.rot = FALSE, 
                      y.lab.rot = FALSE, legend.position = "right", single.legend = TRUE, 
                      remove.legend = FALSE, do.return = FALSE, return.plotlist = FALSE, 
                      ...) 
{
  if (is.null(x = nCol)) {
    if (length(x = features.plot) > 9) {
      nCol <- 4
    }
    else {
      nCol <- min(length(x = features.plot), 3)
    }
  }
  data.use <- data.frame(FetchData(object = object, vars.all = features.plot, 
                                   ...), check.names = F)
  if (is.null(x = ident.include)) {
    cells.to.include <- object@cell.names
  }
  else {
    cells.to.include <- WhichCells(object = object, ident = ident.include)
  }
  data.use <- data.use[cells.to.include, , drop = FALSE]
  if (!is.null(x = group.by)) {
    ident.use <- factor(x = FetchData(object = object, 
                                      vars.all = group.by)[cells.to.include, 1], levels = c.levels)
  }
  else {
    ident.use <- object@ident[cells.to.include]
  }
  gene.names <- colnames(x = data.use)[colnames(x = data.use) %in% 
                                         rownames(x = object@data)]
  if (single.legend) {
    remove.legend <- TRUE
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data.use)
  }
  plots <- lapply(X = features.plot, FUN = function(x) {
    return(SingleVlnPlot(feature = x, data = data.use[, x, 
                                                      drop = FALSE], cell.ident = ident.use, do.sort = do.sort, 
                         y.max = y.max, size.x.use = size.x.use, size.y.use = size.y.use, 
                         size.title.use = size.title.use, adjust.use = adjust.use, 
                         point.size.use = point.size.use, cols.use = cols.use, 
                         gene.names = gene.names, y.log = y.log, x.lab.rot = x.lab.rot, 
                         y.lab.rot = y.lab.rot, legend.position = legend.position, 
                         remove.legend = remove.legend))
  })
  if (length(x = features.plot) > 1) {
    plots.combined <- plot_grid(plotlist = plots, ncol = nCol)
    if (single.legend && !remove.legend) {
      legend <- get_legend(plot = plots[[1]] + theme(legend.position = legend.position))
      if (legend.position == "bottom") {
        plots.combined <- plot_grid(plots.combined, legend, 
                                    ncol = 1, rel_heights = c(1, 0.2))
      }
      else if (legend.position == "right") {
        plots.combined <- plot_grid(plots.combined, legend, 
                                    rel_widths = c(3, 0.3))
      }
      else {
        warning("Shared legends must be at the bottom or right of the plot")
      }
    }
  }
  else {
    plots.combined <- plots[[1]]
  }
  if (do.return) {
    if (return.plotlist) {
      return(plots)
    }
    else {
      return(plots.combined)
    }
  }
  else {
    if (length(x = plots.combined) > 1) {
      plots.combined
    }
    else {
      invisible(x = lapply(X = plots.combined, FUN = print))
    }
  }
}