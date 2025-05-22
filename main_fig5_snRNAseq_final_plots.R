#_--------------------2023-10-02 updated, Jun-Liszt Li; @CIBR, @PKU


#--------------------For main fig.5d--------------------------------#
library(Seurat)

PID15a1 <- readRDS("E:/NGS_data/bAVM_snRNAseq/lijun_analysis/PID15_after_umap_1st.rds")
DimPlot(PID15a1, reduction = "umap", label = T)
DimPlot(PID15c1, reduction = "umap", label = T)



PID27a1 <- readRDS("E:/NGS_data/bAVM_snRNAseq/lijun_analysis/PID27_after_umap_1st.rds")
DimPlot(PID27a1, reduction = "umap", label = T)
DimPlot(PID27c1, reduction = "umap", label = T)



###2023-09-29 re-analysis;
load("E:/NGS_data/bAVM_snRNAseq/lijun_analysis/2023_0409_umap_dim20_PID15_a1_c1_all_merge1.RData")

###2023-09-30 re-analysis;
load("E:/NGS_data/bAVM_snRNAseq/lijun_analysis/2023_0409_umap_dim20_PID27_a1_c1_all_merge1.RData")


#---------------------------------Count number of cells per cluster, PID15--------------------------------------#
# Store cluster identities in object@meta.data$my.clusters
object <- Seurat::StashIdent(object = ECs_PID15.combined, save.name = "my.clusters")

# Get number of cells per cluster and per sample of origin
table(object@meta.data$my.clusters, object@meta.data$orig.ident)
#Cluster  number of cells
#0    120
#1     35
Idents(object = ECs_PID15.combined)


#---------------------------------Count number of cells per cluster, PID27--------------------------------------#
# Store cluster identities in object@meta.data$my.clusters
object <- Seurat::StashIdent(object = ECs_PID27.combined, save.name = "my.clusters")

# Get number of cells per cluster and per sample of origin
table(object@meta.data$my.clusters, object@meta.data$orig.ident)
#Cluster  number of cells
#0    120
#1     35
Idents(object = ECs_PID27.combined)



PID15_EC_meta <- ECs_PID15.combined@meta.data
x1_PID15_EC_meta <- ECs_PID15.combined@meta.data
#> rownames(x1_PID15_EC_meta)[117]
#[1] "c1_Ctrl_GTTCGCTGCTGCCTTCT"
#> rownames(x1_PID15_EC_meta)[116]
#[1] "a1_Exp_ACGGTCGTACCCTCTTC"

PID15_EC_meta$seurat_clusters[1:116] = 0
PID15_EC_meta$seurat_clusters[117:nrow(x1_PID15_EC_meta)] = 1




#--------------------------------------------------PID15--------------------------------------------------------------#
#Using DEPSeq2 to conduct diff expression analysis between ECs_PID27a1 and ECs_PID27c1;

ECs_PID15.combined <- merge(ECs_PID15a1, y = ECs_PID15c1, add.cell.ids = c("a1_Exp", "c1_Ctrl"), project = "ECs_PID15")
ECs_PID15.combined
#An object of class Seurat 
#24565 features across 155 samples within 1 assay 
#Active assay: RNA (24565 features, 0 variable features)
saveRDS(ECs_PID15.combined, file = "E:/NGS_data/bAVM_snRNAseq/lijun_analysis/ECs_PID15.combined.rds")




#------------------------------------Normalizing the data--------------------------------------------#
ECs_PID15.combined <- NormalizeData(ECs_PID15.combined, normalization.method = "LogNormalize", scale.factor = 10000)

#-------------------Identification of highly variable features (feature selection)--------------------#
ECs_PID15.combined <- FindVariableFeatures(ECs_PID15.combined, selection.method = "vst", nfeatures = 2000)


#-------------------------------Scaling the data----------------------------------------#
all.genes <- rownames(ECs_PID15.combined)
ECs_PID15.combined <- ScaleData(ECs_PID15.combined, features = all.genes)
#Centering and scaling data matrix
#|==================================================================================| 100%


#------------Just normalized and scale the data, Do not run non-linear dimensional reduction-------------------#

#---------------------Perform linear dimensional reduction------------------------------#
#Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.
ECs_PID15.combined <- RunPCA(ECs_PID15.combined, features = VariableFeatures(object = ECs_PID15.combined))

#Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()
# Examine and visualize PCA results a few different ways
print(ECs_PID15.combined[["pca"]], dims = 1:5, nfeatures = 5)




#------------------------------Cluster the cells----------------------------------#
ECs_PID15.combined <- FindNeighbors(ECs_PID15.combined, dims = 1:20)
#Computing nearest neighbor graph
#Computing SNN
########################do not need to run the find cluster and choose the resolution; use previous one;

#---------------Run non-linear dimensional reduction (UMAP/tSNE)------------------#
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
ECs_PID15.combined <- RunUMAP(ECs_PID15.combined, dims = 1:20)

###Color: "darkturquoise", "lightcoral"
DimPlot(ECs_PID15.combined, reduction = "umap", label = T, cols = c("darkturquoise", "lightcoral"))




library(dplyr)
# find markers for every cluster compared to all remaining cells, report only the positive ones
ECs_PID15.markers <- FindAllMarkers(ECs_PID15.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ECs_PID15.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


ECs_PID15.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(ECs_PID15.combined, features = top10$gene) + NoLegend()






#----------------------Finding differentially expressed features (cluster biomarkers)------------------------------#
###Function description: Finds markers (differentially expressed genes) for identity classes
# find all markers of cluster 15(bAVM-ECs); Compare Ctrl-EC and bAVM-ECs;
Ecs_PID15a1_exp_diffs <- FindMarkers(ECs_PID15.combined, ident.1 = 15, min.pct = 0.1)
#nrow(Ecs_PID15a1_exp_diffs); 1316 (min.pct = 0.25)
#nrow(Ecs_PID15a1_exp_diffs) ; 3257 (min.pct = 0.1)


head(Ecs_PID15a1_exp_diffs, n = 10)
#p_val avg_log2FC pct.1 pct.2    p_val_adj
#Slco1a4 1.541619e-12  -2.838547 0.103 0.667 3.786987e-08
#Slco2a1 1.597115e-11   3.254060 0.724 0.051 3.923314e-07
#Lgmn    6.319206e-11   2.395185 0.750 0.077 1.552313e-06
#Cxcl12  1.234326e-10  -3.809444 0.026 0.410 3.032121e-06


####Volcano plot
##data$significant <- as.factor(data$PValue<0.05 & abs(data$logFC) > 1)
##data$gene <- data$GENE_NAME
Ecs_PID15a1_exp_diffs$SYMBOL <- rownames(Ecs_PID15a1_exp_diffs)
colnames(Ecs_PID15a1_exp_diffs)[5] <- "padj"
colnames(Ecs_PID15a1_exp_diffs)[2] <- "log2FoldChange"


data <- Ecs_PID15a1_exp_diffs
data$significant <- as.factor(data$padj<0.05 & abs(data$log2FoldChange) > 1)
#########store gene names
data$gene <- Ecs_PID15a1_exp_diffs$SYMBOL

library(ggrepel)
p <- ggplot(data=data, aes(x=log2FoldChange, y =-log10(padj), color=significant)) +
  geom_point(alpha=0.8, size=1.2, col="grey50")+
  geom_point(data=subset(data, data$padj<0.05 & data$log2FoldChange > 1),alpha=0.8, size=1.2,col="#cf558e")+
  geom_point(data=subset(data, data$padj<0.05 & data$log2FoldChange < -1),alpha=0.8, size=1.2,col="#558acf")+
  labs(x="log2(Fold Change)",y="-log10 (padj)", title="Ctrl vs bAVM(PID15)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  theme(legend.position='none')+
  theme_bw()+
  ####remove the border
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_point(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange > 2),alpha=0.8, size=3,col="red") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Rps21"),alpha=0.8, size=3,col="Red") +
  geom_point(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange < -2),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Slc39a10"),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Spock2"),alpha=0.8, size=3,col="blue") +
  geom_text_repel(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange > 2), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange < -2), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Rps21"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Slc39a10"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Spock2"), aes(label=gene), col="black", size = 5.5, alpha = 0.8)

p+ theme(
  title = element_text(color="black", size=19, face="bold"),
  axis.title.x = element_text(color="black", size=19, face="bold"),
  axis.title.y = element_text(color="black", size=19, face="bold"),
  axis.text=element_text(size=18,face = "bold"),
  line = element_line(size =1)) + scale_x_continuous(name="Log2(Fold change)", limits=c(-5, 5)) 


VlnPlot(obj = ECs_PID15.combined, features = c("Flt1", "Erg", "Cldn5", "Slco1a4", "Ackr1", "Tll1", "Slc39a10", "Spock2", "Postn"), stack = TRUE, sort = TRUE, flip = F) +
  theme(legend.position = "none")

VlnPlot(obj = ECs_PID15.combined, features = c("Flt1", "Erg", "Cldn5"), stack = TRUE, sort = TRUE, flip = F) +
  theme(legend.position = "none")




#---------------------------------------------------------Final fig.5k------------------------------------------------------@
#######################Proliferative cells
###Cell cycle analysis
###Warning: All cells have the same value of Ccne1
############
###The following requested variables were not found: Tp53
###The following requested variables were not found: Cdc25
VlnPlot(obj = ECs_PID15.combined, features = c("Ackr1", "Tll1", "Gse1", "Ccdc141", "Nr2f2", "Klf2", "Nrp2", "Mki67", "Cdk1", "Cdk6", "Ccnd3", "Ephb4"), stack = TRUE, sort = T, flip = T) +
  theme(legend.position = "none")



#---------------------------------------------------------Final fig.5q------------------------------------------------------@
#####################Tip cell analysis
##########Common tip cell markers: Npr2; Npr1; Igfbp3; Flt4; Dll4; Robo4; Kdr;
###Stalk cell marker: Jag1; Rhoa; Notch1; Ackr1; Tll1
#----------------------------------------------------------------------------#
###The following requested variables were not found: Vegfr2
############Tip: "Esm1"
###Warning: All cells have the same value of Esm1
##Kdr = Vegfr2
###in primary cell line, "Igf1r", "Igf2",
###The following requested variables were not found: Kdrl
VlnPlot(obj = ECs_PID15.combined, features = c("Ackr1", "Tll1", "Flt4", "Robo4", "Igfbp3", "Dll4", "Kdr", "Nrp1", "Nrp2", "Notch1", "Rhoa", "Jag1"), stack = TRUE, sort = T, flip = TRUE) +
  theme(legend.position = "none")


VlnPlot(obj = ECs_PID15.combined, features = c("Flt4", "Robo4", "Igfbp3", "Dll4", "Kdr", "Nrp1", "Nrp2", "Notch1", "Rhoa", "Jag1"), stack = TRUE, sort = T, flip = TRUE) +
  theme(legend.position = "none")





#----------------------------------------------For main fig.5g-------------------------------------------------------#
#--------------------------------------------------PID27--------------------------------------------------------------#
#Using DEPSeq2 to conduct diff expression analysis between ECs_PID27a1 and ECs_PID27c1;

ECs_PID27.combined <- merge(ECs_PID27a1, y = ECs_PID27c1, add.cell.ids = c("a1_Exp", "c1_Ctrl"), project = "ECs_PID27")
ECs_PID27.combined
#An object of class Seurat 
#24211 features across 522 samples within 1 assay 
#Active assay: RNA (24211 features, 0 variable features)
saveRDS(ECs_PID27.combined, file = "E:/NGS_data/bAVM_snRNAseq/lijun_analysis/ECs_PID27.combined.rds")


#------------------------------------Normalizing the data--------------------------------------------#
ECs_PID27.combined <- NormalizeData(ECs_PID27.combined, normalization.method = "LogNormalize", scale.factor = 10000)

#-------------------Identification of highly variable features (feature selection)--------------------#
ECs_PID27.combined <- FindVariableFeatures(ECs_PID27.combined, selection.method = "vst", nfeatures = 2000)


#-------------------------------Scaling the data----------------------------------------#
all.genes <- rownames(ECs_PID27.combined)
ECs_PID27.combined <- ScaleData(ECs_PID27.combined, features = all.genes)
#Centering and scaling data matrix


#------------Just normalized and scale the data, Do not run non-linear dimensional reduction-------------------#

#---------------------Perform linear dimensional reduction------------------------------#
#Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.
ECs_PID27.combined <- RunPCA(ECs_PID27.combined, features = VariableFeatures(object = ECs_PID27.combined))

#Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()
# Examine and visualize PCA results a few different ways
print(ECs_PID27.combined[["pca"]], dims = 1:5, nfeatures = 5)



#------------------------------Cluster the cells----------------------------------#
ECs_PID27.combined <- FindNeighbors(ECs_PID27.combined, dims = 1:20)
#Computing nearest neighbor graph
#Computing SNN
########################do not need to run the find cluster and choose the resolution; use previous one;

#---------------Run non-linear dimensional reduction (UMAP/tSNE)------------------#
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
ECs_PID27.combined <- RunUMAP(ECs_PID27.combined, dims = 1:20)

DimPlot(ECs_PID27.combined, reduction = "umap", label = T)




#---------------------------------Count number of cells per cluster--------------------------------------#
# Store cluster identities in object@meta.data$my.clusters
object <- Seurat::StashIdent(object = ECs_PID27.combined, save.name = "my.clusters")

# Get number of cells per cluster and per sample of origin
table(object@meta.data$my.clusters, object@meta.data$orig.ident)
#bAVM27
#24     33
#8     489



###########################################################################


PID27_EC_meta <- ECs_PID27.combined@meta.data
x1_PID27_EC_meta <- ECs_PID27.combined@meta.data

PID27_EC_meta$seurat_clusters[1:490] = 0
PID27_EC_meta$seurat_clusters[490:522] = 1


PID27_EC_meta$RNA_snn_res.0.4[1:490] = 0
PID27_EC_meta$RNA_snn_res.0.4[490:522] = 1



#----------------------Finding differentially expressed features (cluster biomarkers)------------------------------#
# find all markers of cluster 4; Compare Ctrl-EC(ident = 24) and bAVM-ECs(ident = 8);
Ecs_PID27a1_exp_diffs <- FindMarkers(ECs_PID27.combined, ident.1 = c(8), min.pct = 0.1)
nrow(Ecs_PID27a1_exp_diffs)
###min.pct = 0.25: [1] 1149
head(Ecs_PID27a1_exp_diffs, n = 10)




####Volcano plot
##data$significant <- as.factor(data$PValue<0.05 & abs(data$logFC) > 1)
##data$gene <- data$GENE_NAME
Ecs_PID27a1_exp_diffs$SYMBOL <- rownames(Ecs_PID27a1_exp_diffs)
colnames(Ecs_PID27a1_exp_diffs)[5] <- "padj"
colnames(Ecs_PID27a1_exp_diffs)[2] <- "log2FoldChange"


data <- Ecs_PID27a1_exp_diffs
data$significant <- as.factor(data$padj<0.05 & abs(data$log2FoldChange) > 1)
#########store gene names
data$gene <- Ecs_PID27a1_exp_diffs$SYMBOL

library(ggrepel)
p <- ggplot(data=data, aes(x=log2FoldChange, y =-log10(padj), color=significant)) +
  geom_point(alpha=0.8, size=1.2, col="grey50")+
  geom_point(data=subset(data, data$padj<0.05 & data$log2FoldChange > 1),alpha=0.8, size=1.2,col="#cf558e")+
  geom_point(data=subset(data, data$padj<0.05 & data$log2FoldChange < -1),alpha=0.8, size=1.2,col="#558acf")+
  labs(x="log2(Fold Change)",y="-log10 (padj)", title="Ctrl vs bAVM(PID27)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  theme(legend.position='none')+
  theme_bw()+
  ####remove the border
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_point(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange > 1.8),alpha=0.8, size=3,col="red") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Tll1"),alpha=0.8, size=3,col="Red") +
  #geom_point(data=dplyr::filter(data, data$SYMBOL == "Slco2a1"),alpha=0.8, size=3,col="Red") +
  geom_point(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange < -2),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Slc39a10"),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Spock2"),alpha=0.8, size=3,col="blue") +
  geom_text_repel(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange > 1.8), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange < -2), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Tll1"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  #geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Slco2a1"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Slc39a10"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Spock2"), aes(label=gene), col="black", size = 5.5, alpha = 0.8)

p+ theme(
  title = element_text(color="black", size=19, face="bold"),
  axis.title.x = element_text(color="black", size=19, face="bold"),
  axis.title.y = element_text(color="black", size=19, face="bold"),
  axis.text=element_text(size=18,face = "bold"),
  line = element_line(size =1)) + scale_x_continuous(name="Log2(Fold change)", limits=c(-5, 5)) 


VlnPlot(obj = ECs_PID27.combined, features = c("Tll1", "Ackr1", "Postn", "Vcam1", "Spock2", "Slc39a10", "Slco1a4"), stack = TRUE, sort = F, flip = T) +
  theme(legend.position = "none")

VlnPlot(obj = ECs_PID27.combined, features = c("Flt1", "Erg", "Cldn5"), stack = TRUE, sort = F, flip = F) +
  theme(legend.position = "none")


#######################Proliferative cells
###Cell cycle analysis
###Warning: All cells have the same value of Ccne1
############
###The following requested variables were not found: Tp53
###The following requested variables were not found: Cdc25

###Do not sort;
VlnPlot(obj = ECs_PID27.combined, features = c("Ephb4", "Cdk6", "Tll1", "Ackr1", "Gse1", "Nrp2", "Mki67", "Cdk1", "Nr2f2", "Klf2", "Ccnd3", "Ccdc141"), stack = TRUE, sort = F, flip = T) +
  theme(legend.position = "none")



#####################Tip cell analysis
##########Common tip cell markers:
###Stalk cell marker: Jag1
###The following requested variables were not found: Vegfr2
############Tip: "Esm1"
###Warning: All cells have the same value of Esm1
##Kdr = Vegfr2
###in primary cell line, "Igf1r", "Igf2",
###The following requested variables were not found: Kdrl
###The arrangement of features is set; please do not disturb it;
###do not sort;
VlnPlot(obj = ECs_PID27.combined, features = c("Nrp2", "Nrp1", "Tll1", "Ackr1", "Igfbp3", "Flt4", "Notch1", "Dll4", "Robo4", "Kdr", "Rhoa", "Jag1"), stack = TRUE, sort = F, flip = TRUE) +
  theme(legend.position = "none")


###
VlnPlot(obj = ECs_PID15.combined, features = c("Nrp2", "Nrp1", "Igfbp3", "Flt4", "Notch1", "Dll4", "Robo4", "Kdr", "Rhoa", "Jag1"), stack = TRUE, sort = T, flip = TRUE) +
  theme(legend.position = "none")


###Compare with human bAVM snRNAseq data; 
###Nidus ECs diff genes;
###再mapping下Mfsd2a, Slc16a1, Slc38a5.
###Upregulated:Itga6, Ccl14, Pgf, Stc1, Plvap, Angpt2.
###Downregulated: Timp1

###PID27
VlnPlot(obj = ECs_PID27.combined, features = c("Tll1", "Ackr1", "Postn", "Vcam1", "Ccl14", "Pgf", "Stc1", "Plvap", "Angpt2"), stack = TRUE, sort = F, flip = T) +
  theme(legend.position = "none")

###PID15
VlnPlot(obj = ECs_PID15.combined, features = c("Itga6", "Tll1", "Ackr1", "Postn", "Vcam1", "Ccl14", "Pgf", "Plvap", "Angpt2"), stack = TRUE, sort = T, flip = T) +
  theme(legend.position = "none")

###
VlnPlot(obj = ECs_PID15.combined, features = c("Timp1", "Tll1", "Ackr1", "Postn", "Vcam1", "Rgcc", "Pgf", "Plvap", "Angpt2"), stack = TRUE, sort = F, flip = T) +
  +   theme(legend.position = "none")





###Capillary marker;
VlnPlot(obj = ECs_PID15.combined, features = c("Slco1c1", "Slc7a5", "Mfsd2a"), stack = TRUE, sort = T, flip = T) +
  theme(legend.position = "none")
VlnPlot(obj = ECs_PID27.combined, features = c("Slco1c1", "Slc7a5", "Mfsd2a"), stack = TRUE, sort = F, flip = T) +
  theme(legend.position = "none")





#------------------------------------------------------fig.5f-------------------------------------------------------#
###2023-10-07updated
data <- Ecs_PID15a1_exp_diffs
data$significant <- as.factor(data$padj<0.05 & abs(data$log2FoldChange) > 1)
#########store gene names
data$gene <- Ecs_PID15a1_exp_diffs$SYMBOL

#############################padj as y axis;
p <- ggplot(data=data, aes(x=log2FoldChange, y =-log10(padj), color=significant)) +
  geom_point(alpha=0.8, size=1.2, col="grey50")+
  geom_point(data=subset(data, data$padj<0.05 & data$log2FoldChange > 1),alpha=0.8, size=1.2,col="#cf558e")+
  geom_point(data=subset(data, data$padj<0.05 & data$log2FoldChange < -1),alpha=0.8, size=1.2,col="#558acf")+
  labs(x="log2(Fold Change)",y="-log10 (padj)", title="Ctrl vs bAVM(PID15)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  theme(legend.position='none')+
  theme_bw()+
  ####remove the border
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  #geom_point(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange > 2),alpha=0.8, size=3,col="red") +
  geom_point(data=dplyr::filter(data, data$SYMBOL  == "Rhoj"),alpha=0.8, size=3,col="Red") +
  geom_point(data=dplyr::filter(data, data$SYMBOL  == "Slco2a1"),alpha=0.8, size=3,col="Red") +
  geom_point(data=dplyr::filter(data, data$SYMBOL  == "Tacr1"),alpha=0.8, size=3,col="Red") +
  geom_point(data=dplyr::filter(data, data$SYMBOL  == "Cmip"),alpha=0.8, size=3,col="Red") +
  #geom_point(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange < -2),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Slco1a4"),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Atp10a"),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Adipor2"),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Cxcl12"),alpha=0.8, size=3,col="blue") +
  #geom_text_repel(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange > 1.8), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  #geom_text_repel(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange < -2), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  #geom_text_repel(data=dplyr::filter(data, data$SYMBOL %in% Pdgf_list_gene_sig$SYMBOL), aes(label=gene), col="black", size = 5.5, alpha = 0.8) + 
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Rhoj"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) + 
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Slco2a1"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Tacr1"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Cmip"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  ###
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Slco1a4"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Atp10a"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) + 
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Adipor2"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Cxcl12"), aes(label=gene), col="black", size = 5.5, alpha = 0.8)

#geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Slc39a10"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
#geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Spock2"), aes(label=gene), col="black", size = 5.5, alpha = 0.8)

p+ theme(
  title = element_text(color="black", size=19, face="bold"),
  axis.title.x = element_text(color="black", size=19, face="bold"),
  axis.title.y = element_text(color="black", size=19, face="bold"),
  axis.text=element_text(size=18,face = "bold"),
  line = element_line(size =1)) + scale_x_continuous(name="Log2(Fold change)", limits=c(-4, 4)) 




#---------------------------------------------------PID27-----------------------------------------------------
#--------------------------------------------------fig.5i----------------------------------------------------------#
###2023-10-07updated
data <- Ecs_PID27a1_exp_diffs
data$significant <- as.factor(data$padj<0.05 & abs(data$log2FoldChange) > 1)
#########store gene names
data$gene <- Ecs_PID27a1_exp_diffs$SYMBOL
#############################padj as y axis;
p <- ggplot(data=data, aes(x=log2FoldChange, y =-log10(padj), color=significant)) +
  geom_point(alpha=0.8, size=1.2, col="grey50")+
  geom_point(data=subset(data, data$padj<0.05 & data$log2FoldChange > 1),alpha=0.8, size=1.2,col="#cf558e")+
  geom_point(data=subset(data, data$padj<0.05 & data$log2FoldChange < -1),alpha=0.8, size=1.2,col="#558acf")+
  labs(x="log2(Fold Change)",y="-log10 (padj)", title="Ctrl vs bAVM(PID27)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  theme(legend.position='none')+
  theme_bw()+
  ####remove the border
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  #geom_point(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange > 2),alpha=0.8, size=3,col="red") +
  geom_point(data=dplyr::filter(data, data$SYMBOL  == "Nav3"),alpha=0.8, size=3,col="Red") +
  geom_point(data=dplyr::filter(data, data$SYMBOL  == "Slco2a1"),alpha=0.8, size=3,col="Red") +
  geom_point(data=dplyr::filter(data, data$SYMBOL  == "Tacr1"),alpha=0.8, size=3,col="Red") +
  #geom_point(data=dplyr::filter(data, data$SYMBOL  == "Smad1"),alpha=0.8, size=3,col="Red") +
  #geom_point(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange < -2),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Slco1a4"),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Slc39a10"),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Slc2a1"),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Bsg"),alpha=0.8, size=3,col="blue") +
  #geom_text_repel(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange > 1.8), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  #geom_text_repel(data=subset(data, data$padj < 1*10^(-3) & data$log2FoldChange < -2), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  #geom_text_repel(data=dplyr::filter(data, data$SYMBOL %in% Pdgf_list_gene_sig$SYMBOL), aes(label=gene), col="black", size = 5.5, alpha = 0.8) + 
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Nav3"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) + 
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Slco2a1"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Tacr1"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  #geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Smad1"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  ###
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Slco1a4"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Slc39a10"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) + 
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Slc2a1"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Bsg"), aes(label=gene), col="black", size = 5.5, alpha = 0.8)

#geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Slc39a10"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
#geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Spock2"), aes(label=gene), col="black", size = 5.5, alpha = 0.8)

p+ theme(
  title = element_text(color="black", size=19, face="bold"),
  axis.title.x = element_text(color="black", size=19, face="bold"),
  axis.title.y = element_text(color="black", size=19, face="bold"),
  axis.text=element_text(size=18,face = "bold"),
  line = element_line(size =1)) + scale_x_continuous(name="Log2(Fold change)", limits=c(-5, 5)) 





#----------------------------------------------main fig.5e and 5h---------------------------------------------------------#
VlnPlot(obj = ECs_PID15.combined, features = c("Erg", "Cldn5", "Flt1"), stack = TRUE, sort = TRUE, flip = F) +
  theme(legend.position = "none")

VlnPlot(obj = ECs_PID27.combined, features = c("Erg", "Cldn5", "Flt1"), stack = TRUE, sort = F, flip = F) +
  theme(legend.position = "none")











#---------------------------------------------For Supplementary Figure 40a---------------------------------------------#

#------------------------------------------------For GSEA analysis----------------------------------------------#
###To generate a pre-ranked gene list(for RNAseq, go from raw counts to DESeq2):
library(DESeq2)
library(readr)
## I have 20 samples for each condition: knockout gene X, and WT gene X
## KO vs WT
#dds$condition <- factor(dds$condition, levels = c("KO","WT"))
#res <- results(dds, contrast=c("condition", "KO", "WT"))

#########################
data$ENSEMBL<- toupper(data$SYMBOL)      ##########covert to uppercase using toupper function;


res_pre_rank_EC <- sign(data$log2FoldChange) * -log10(data$p_val)
res_pre_rank_EC[is.na(res_pre_rank_EC)] = 0        ######convert na into "0"
# you will need the gene symbol, and suffix the file with rnk for GSEA to recognize
write_tsv(data.frame(Name = data$ENSEMBL, metric = res_pre_rank_EC), "L:/ECs_PID15_bAVM_prerank_2023_10_03updated.rnk")


#------------------------------------------------GSEA analysis----------------------------------------------#
###To generate a pre-ranked gene list(for RNAseq, go from raw counts to DESeq2):
library(DESeq2)
library(readr)
## I have 20 samples for each condition: knockout gene X, and WT gene X
## KO vs WT
#dds$condition <- factor(dds$condition, levels = c("KO","WT"))
#res <- results(dds, contrast=c("condition", "KO", "WT"))

#########################
data$ENSEMBL<- data$gene      ##########covert to uppercase using toupper function;


res_pre_rank_EC <- sign(data$log2FoldChange) * -log10(data$p_val)
res_pre_rank_EC[is.na(res_pre_rank_EC)] = 0        ######convert na into "0"
# you will need the gene symbol, and suffix the file with rnk for GSEA to recognize
write_tsv(data.frame(Name = data$ENSEMBL, metric = res_pre_rank_EC), "K:/lijun_bAVM_snRNAseq/GSEA_analysis/human_bAVM_EC_prerank_2023_10_03updated.rnk")






#-----------------------------------------For main fig.5j and supplementary fig.40a----------------------------------------#


#----------------------------------------For comparison of human bAVM-ECs dataset------------------------------------------#
#-----------------------------------------------Pre-process-----------------------------------------------------#
{
###########################################load human bAVM data;
require(Seurat)
require(data.table)
#setwd("adultPancreas")
mat <- fread("L:/NGS_data_all/ref_2022_science_human_bAVM_snRNAseq/exprMatrix.tsv.gz")
meta <- read.table("L:/NGS_data_all/ref_2022_science_human_bAVM_snRNAseq/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
h_bAVM <- CreateSeuratObject(counts = mat, project = "human_bAVM", meta.data=meta)
#An object of class Seurat 
#33495 features across 9541 samples within 1 assay 
#Active assay: RNA (33495 features, 0 variable features)


#------------------------------------Normalizing the data--------------------------------------------#
h_bAVM <- NormalizeData(h_bAVM, normalization.method = "LogNormalize", scale.factor = 10000)

#-------------------Identification of highly variable features (feature selection)--------------------#
h_bAVM <- FindVariableFeatures(h_bAVM, selection.method = "vst", nfeatures = 2000)


#-------------------------------Scaling the data----------------------------------------#
all.genes <- rownames(h_bAVM)
h_bAVM <- ScaleData(h_bAVM, features = all.genes)
#Centering and scaling data matrix
#|==================================================================================| 100%

#------------Just normalized and scale the data, Do not run non-linear dimensional reduction-------------------#



#---------------------Perform linear dimensional reduction------------------------------#
#Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.
h_bAVM <- RunPCA(h_bAVM, features = VariableFeatures(object = h_bAVM))

#Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()
# Examine and visualize PCA results a few different ways
print(h_bAVM[["pca"]], dims = 1:5, nfeatures = 5)



VizDimLoadings(h_bAVM, dims = 1:2, reduction = "pca")
DimPlot(h_bAVM, reduction = "pca")



DimHeatmap(h_bAVM, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(h_bAVM, dims = 2, cells = 500, balanced = TRUE)

DimHeatmap(h_bAVM, dims = 1:9, cells = 500, balanced = TRUE)


#----------------Determine the ‘dimensionality’ of the dataset---------------------#
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
h_bAVM <- JackStraw(h_bAVM, num.replicate = 100)
h_bAVM <- ScoreJackStraw(h_bAVM, dims = 1:20)

JackStrawPlot(h_bAVM, dims = 1:20)
#Warning message:
#  Removed 28000 rows containing missing values (geom_point). 

ElbowPlot(h_bAVM)


#------------------------------Cluster the cells----------------------------------#
h_bAVM <- FindNeighbors(h_bAVM, dims = 1:20)
#Computing nearest neighbor graph
#Computing SNN
########################do not need to run the find cluster and choose the resolution; use previous one;


meta_h_bAVM <- h_bAVM@meta.data
###Ctrl:[1:4940]; AVM:[4941:9541]


#h_bAVM <- FindClusters(h_bAVM, resolution = 0.5)


# Look at cluster IDs of the first 5 cells
head(Idents(h_bAVM), 5)



#---------------Run non-linear dimensional reduction (UMAP/tSNE)------------------#
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
h_bAVM <- RunUMAP(h_bAVM, dims = 1:20)

DimPlot(h_bAVM, reduction = "umap", label = T)



#---------------------------------Count number of cells per cluster--------------------------------------#
# Store cluster identities in object@meta.data$my.clusters
object <- Seurat::StashIdent(object = h_bAVM, save.name = "my.clusters")

# Get number of cells per cluster and per sample of origin
table(object@meta.data$my.clusters, object@meta.data$orig.ident)
#8 AVM patients, 7 control tissues;

colnames(table(object@meta.data$my.clusters, object@meta.data$orig.ident))



# find markers for every cluster compared to all remaining cells, report only the positive
# ones
h_bAVM.markers <- FindAllMarkers(h_bAVM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
h_bAVM.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#
h_bAVM.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(h_bAVM, features = top5$gene) + NoLegend()


############Brain endothelial cells
FeaturePlot(object = h_bAVM, features = c("CLDN5", "SLC2A1", "CDH5"), cols = c("grey", "red"), raster=FALSE, label = T) 
StackedVlnPlot(obj = h_bAVM, features = c("Cldn5", "Slc2a1", "Cdh5"))


FeaturePlot(object = h_bAVM, features = c("SLCO1C1"), cols = c("grey", "red"), raster=FALSE, label = T) 


VlnPlot(obj = PID15, features = c("Slc1a3", "Slco1c1"), stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")

############################2023-04-09updated
VlnPlot(obj = h_bAVM, features = c("CLDN5", "GFAP", "SLC38A5", "SLC2A1", "CDH5", "SLCO1C1", "SLCO1A4", "FLT1", "ERG", "PECAM1", "PDGFRB", "CSPG4", "KCNJ8", "ACTA2","PDLIM3", "AQP4"), stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none")


h_bAVM_meta_data <- h_bAVM@meta.data



#-------------------------------------------------------For main fig.5j-----------------------------------------------------#
hbAVM_ECs_rds <- readRDS("L:/NGS_data_all/ref_2022_science_human_bAVM_snRNAseq/final_merged_ec_cb.rds")
#View(hbAVM_ECs_rds)
DimPlot(hbAVM_ECs_rds, reduction = "umap", label = T)


#---------------------------------Count number of cells per cluster------------------------------------------------#
# Store cluster identities in object@meta.data$my.clusters
object <- Seurat::StashIdent(object = hbAVM_ECs_rds, save.name = "my.clusters")

# Get number of cells per cluster and per sample of origin
table(object@meta.data$my.clusters, object@meta.data$orig.ident)



h_bAVM_meta <- hbAVM_ECs_rds@meta.data
nidus_meta <- dplyr::filter(h_bAVM_meta, h_bAVM_meta$clusters == "Nidus")


####################find nidus specific differential expressed genes;
#----------------------Finding differentially expressed features (cluster biomarkers)------------------------------#
# find all markers of cluster 4; Compare Ctrl-EC(ident = 24) and bAVM-ECs(ident = 8);
hbAVM_Nidus_exp_diffs <- FindMarkers(hbAVM_ECs_rds, ident.1 = c("Nidus"), min.pct = 0.25)
head(hbAVM_Nidus_exp_diffs, n = 10)


####Volcano plot
##data$significant <- as.factor(data$PValue<0.05 & abs(data$logFC) > 1)
##data$gene <- data$GENE_NAME
hbAVM_Nidus_exp_diffs$SYMBOL <- rownames(hbAVM_Nidus_exp_diffs)
colnames(hbAVM_Nidus_exp_diffs)[5] <- "padj"
colnames(hbAVM_Nidus_exp_diffs)[2] <- "log2FoldChange"


data <- hbAVM_Nidus_exp_diffs
data$significant <- as.factor(data$padj<0.05 & abs(data$log2FoldChange) > 1)
#########store gene names
data$gene <- hbAVM_Nidus_exp_diffs$SYMBOL

library(ggrepel)
p <- ggplot(data=data, aes(x=log2FoldChange, y =-log10(padj), color=significant)) +
  geom_point(alpha=0.8, size=1.2, col="grey50")+
  geom_point(data=subset(data, data$padj<0.05 & data$log2FoldChange > 1),alpha=0.8, size=1.2,col="#cf558e")+
  geom_point(data=subset(data, data$padj<0.05 & data$log2FoldChange < -1),alpha=0.8, size=1.2,col="#558acf")+
  labs(x="log2(Fold Change)",y="-log10 (padj)", title="Others vs Nidus(h-bAVM)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  theme(legend.position='none')+
  theme_bw()+
  ####remove the border
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_point(data=subset(data, data$padj < 1*10^(-2) & data$log2FoldChange > 1),alpha=0.8, size=3,col="red") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Tll1"),alpha=0.8, size=3,col="Red") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Ackr1"),alpha=0.8, size=3,col="Red") +
  geom_point(data=subset(data, data$padj < 1*10^(-2) & data$log2FoldChange < -1),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Slc39a10"),alpha=0.8, size=3,col="blue") +
  geom_point(data=dplyr::filter(data, data$SYMBOL == "Spock2"),alpha=0.8, size=3,col="blue") +
  geom_text_repel(data=subset(data, data$padj < 1*10^(-2) & data$log2FoldChange > 1), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=subset(data, data$padj < 1*10^(-2) & data$log2FoldChange < -1), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Tll1"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Ackr1"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Slc39a10"), aes(label=gene), col="black", size = 5.5, alpha = 0.8) +
  geom_text_repel(data=dplyr::filter(data, data$SYMBOL == "Spock2"), aes(label=gene), col="black", size = 5.5, alpha = 0.8)

p+ theme(
  title = element_text(color="black", size=19, face="bold"),
  axis.title.x = element_text(color="black", size=19, face="bold"),
  axis.title.y = element_text(color="black", size=19, face="bold"),
  axis.text=element_text(size=18,face = "bold"),
  line = element_line(size =1)) + scale_x_continuous(name="Log2(Fold change)", limits=c(-2.5, 2.5)) 



VlnPlot(obj = hbAVM_ECs_rds, features = c("TLL1", "ACKR1", "POSTN", "VCAM1"), stack = TRUE, sort = F, flip = T) +
  theme(legend.position = "none")

VlnPlot(obj = hbAVM_ECs_rds, features = c("TLL1", "PLVAP"), stack = TRUE, sort = F, flip = T) +
  theme(legend.position = "none")




AVM_ctrl_separate_rds <- hbAVM_ECs_rds

#########################re-cluster in to two;
AVM_Ctrl_separ_meta <- AVM_ctrl_separate_rds@meta.data
AVM_ctrl_separate_rds@meta.data$clusters[1:4940] <- "Ctrl"
AVM_ctrl_separate_rds@meta.data$clusters[4940:nrow(AVM_Ctrl_separ_meta)] <- "AVM"



DimPlot(AVM_ctrl_separate_rds, reduction = "umap", label = T)


#----------------------Finding differentially expressed features (cluster biomarkers)------------------------------#
# find all markers of cluster 4; Compare Ctrl-EC(ident = 24) and bAVM-ECs(ident = 8);
hbAVMxx_exp_diffs <- FindMarkers(AVM_ctrl_separate_rds, ident.1 = c("AVM"), min.pct = 0.25)
head(hbAVM_Nidus_exp_diffs, n = 10)


human_bAVM_EC_diffs <- read.csv("L:/NGS_data_all/ref_2022_science_human_bAVM_snRNAseq/hbAVM_S6_EC_diffs_AVM_ctrl.csv")



##############################################
colnames(human_bAVM_EC_diffs)[6] <- "padj"
colnames(human_bAVM_EC_diffs)[3] <- "log2FoldChange"


data <- human_bAVM_EC_diffs
data$significant <- as.factor(data$padj<0.05 & abs(data$log2FoldChange) > 1)
#########store gene names
data$gene <- human_bAVM_EC_diffs$gene



####################################convert to mouse homolog gene
library(Hmisc)
human_bAVM_EC_diffs$mouse_homolog <- capitalize(tolower(human_bAVM_EC_diffs$gene))



#################Screen out the significant genes
human_bAVM_DEGs_sig <- dplyr::filter(human_bAVM_EC_diffs, human_bAVM_EC_diffs$p_val_adj<0.1 & abs(human_bAVM_EC_diffs$avg_log2FC) > 0.5)
nrow(human_bAVM_DEGs_sig)
#[1]  2665  (abs(avg_log2FC) > 0.75))
#[1] 5140    (abs(avg_log2FC) > 0.5))

####Human up and down
human_bAVM_EC_sigs_up <- dplyr::filter(human_bAVM_DEGs_sig, human_bAVM_DEGs_sig$avg_log2FC>0)
human_bAVM_EC_sigs_dw <- dplyr::filter(human_bAVM_DEGs_sig, human_bAVM_DEGs_sig$avg_log2FC<0)
genes_bAVM_EC_sigs_up <- human_bAVM_EC_sigs_up$mouse_homolog
length(genes_bAVM_EC_sigs_up)
###[1] 2660
genes_bAVM_EC_sigs_dw <- human_bAVM_EC_sigs_dw$mouse_homolog
length(genes_bAVM_EC_sigs_dw)
###[1] 2660



#------------------------------------------------------------PID15---------------------------------------------------------------#
library(dplyr)
#################Screen out the significant genes
PID15_DEGs_sig <- dplyr::filter(Ecs_PID15a1_exp_diffs, Ecs_PID15a1_exp_diffs$p_val_adj<0.1 & abs(Ecs_PID15a1_exp_diffs$avg_log2FC) > 0.5)
nrow(PID15_DEGs_sig)
#[1] 1708


####Mouse PID15 up and down
Ecs_PID15a1_exp_diffs$SYMBOL <- rownames(Ecs_PID15a1_exp_diffs)
Ecs_PID15a1_exp_diffs_up <- dplyr::filter(PID15_DEGs_sig, PID15_DEGs_sig$avg_log2FC>0)
Ecs_PID15a1_exp_diffs_dw <- dplyr::filter(PID15_DEGs_sig, PID15_DEGs_sig$avg_log2FC<0)
genes_mus_PID15_sig_up <- Ecs_PID15a1_exp_diffs_up$SYMBOL
length(genes_mus_PID15_sig_up)
#[1] 1012
genes_mus_PID15_sig_dw <- Ecs_PID15a1_exp_diffs_dw$SYMBOL
length(genes_mus_PID15_sig_dw)
#[1] 696




###identify the overlaped gene sets between human bAVM-EC and mouse bAVM-ECs (both PID15 and PID17)
PID15_overlap_up_genes <- intersect(human_bAVM_EC_sigs_up$mouse_homolog, genes_mus_PID15_sig_up)
length(PID15_overlap_up_genes)
###95/413
###229/1012
#############[1] "Rhoj"    "Magi1"   "Ext1"    "Slco2a1" "Myof"
###Up regulated overlapping gene number: n = 249 

PID15_overlap_dw_genes <- intersect(human_bAVM_EC_sigs_dw$mouse_homolog, genes_mus_PID15_sig_dw)
length(PID15_overlap_dw_genes)
#84/324
##165/696
###[1] 3; "Slc2a1" "Atp10a" "Cxcl12"
###Down regulated overlapping gene number: n = 156 




#------------------------------------------------------------PID27---------------------------------------------------------------#
#################Screen out the significant genes
PID27_DEGs_sig <- dplyr::filter(Ecs_PID27a1_exp_diffs, Ecs_PID27a1_exp_diffs$p_val_adj<5 & abs(Ecs_PID27a1_exp_diffs$avg_log2FC) > 0.5)
nrow(PID27_DEGs_sig)
#[1] 1672
####Mouse PID27 up and down
Ecs_PID27a1_exp_diffs$SYMBOL <- rownames(Ecs_PID27a1_exp_diffs)
Ecs_PID27a1_exp_diffs_up <- dplyr::filter(PID27_DEGs_sig, PID27_DEGs_sig$avg_log2FC>0)
Ecs_PID27a1_exp_diffs_dw <- dplyr::filter(PID27_DEGs_sig, PID27_DEGs_sig$avg_log2FC<0)
genes_mus_PID27_sig_up <- Ecs_PID27a1_exp_diffs_up$SYMBOL
length(genes_mus_PID27_sig_up)
#900
genes_mus_PID27_sig_dw <- Ecs_PID27a1_exp_diffs_dw$SYMBOL
length(genes_mus_PID27_sig_dw)
#772


PID27_overlap_up_genes <- intersect(human_bAVM_EC_sigs_up$mouse_homolog, genes_mus_PID27_sig_up)
length(PID27_overlap_up_genes)
###[1] 42
###118/900

PID27_overlap_dw_genes <- intersect(human_bAVM_EC_sigs_dw$mouse_homolog, genes_mus_PID27_sig_dw)
length(PID27_overlap_dw_genes)
###[1] 80
###151/772






#-------------------------------------------------------For main fig.5j-----------------------------------------------------#
########PID15
#-----------------------------------------------GO enrichment analysis-----------------------------------------------------#
#----------Extract the data table for GO enrichment analysis;
Ecs_PID15a1_human_overlap_up <- dplyr::filter(Ecs_PID15a1_exp_diffs_up, Ecs_PID15a1_exp_diffs_up$SYMBOL %in% PID15_overlap_up_genes)


library(org.Mm.eg.db)
library(clusterProfiler) 
gene_ids <- bitr(Ecs_PID15a1_human_overlap_up$gene, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Mm.eg.db")
# I used ENSEMBL gene id for the matrix, and annotate with gene symbols.
#res_anno<- as.data.frame(res_3v_dw_up) %>% left_join(gene_symbols, c("ENSEMBL" = "gene_id"))
Ecs_PID15a1_human_overlap_up <- left_join(Ecs_PID15a1_human_overlap_up, gene_ids, c("gene" = "SYMBOL"))

########ID choose 0.001; 4 fold
res_de_bAVM_upx_id = data.frame(ID=Ecs_PID15a1_human_overlap_up$ENSEMBL) 


#------------------------------------------up regulated genes enrichment analysis-----------------------------------------------#
library(DOSE) 
library(annotate)
eg = bitr(de_bAVM_x_id$ID, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
############Cellular function;
ego_cc<-enrichGO(gene = res_de_bAVM_upx_id$ID,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'ENSEMBL',
                 ont        = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.01,
                 readable      = TRUE) #TRUE 
ego_cc_results<-summary(ego_cc)
barplot(ego_cc,showCategory = 15, title="The GO_CC of DEGs(up_0-01_4fc)")
dotplot(ego_cc, showCategory = 15, title="The GO_CC of DEGs(up_0-01_4fc)")
#######dotplot(ego_cc, showCategory = 20, split=".sign", title="The GO_CC enrichment analysis of all DEGs") + facet_grid(.~.sign)
#  scale_size(range=c(2, 12))+
#  scale_x_discrete(labels=function(ego_cc) str_wrap(ego_cc,width = 25))
library(topGO)
library(GSEABase)
##########################################
plotGOgraph(ego_cc)
##igraph format of plot GOgraph
goplot(ego_cc)
##go term overlapping relationship
emapplot(ego_cc, showCategory = 30)
##gene and go term
cnetplot(ego_cc, showCategory = 5)
##enrichmap
library(enrichplot)
require(doseplot)
library(igraph)
###enrichMap(ego_cc)
####molecular function


#####################biological process;
ego_bp<-enrichGO(gene       = res_de_bAVM_upx_id$ID,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'ENSEMBL',
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05,
                 readable      = TRUE) #TRUE shows SYMBOL，FALSE shows original ID)
ego_bp_results<-summary(ego_bp)
barplot(ego_bp,showCategory = 15,title="The GO_BP of DEGs(up_0-01_4fc)")
##scale_size(range=c(2, 12))+
##scale_x_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 40))
dotplot(ego_bp, showCategory = 15, title="The GO_BP of DEGs(up_0-01_4fc)")
#+  scale_size(range=c(2, 12))


######################molecular function;
ego_mf_up<- enrichGO(gene = res_de_bAVM_upx_id$ID,
                     OrgDb      = org.Mm.eg.db,
                     keyType    = 'ENSEMBL',
                     ont        = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable      = TRUE) #TRUE shows SYMBOL，FALSE shows original ID)
ego_mf_up_results<-summary(ego_mf_up)
barplot(ego_mf_up,showCategory = 15,title="The GO_MF of DEGs(up_0-05_2fc)")
dotplot(ego_mf_up, showCategory = 15, title="The GO_MF of DEGs(up_0-05_2fc)")
#  scale_size(range=c(2, 12))+
#  scale_x_discrete(labels=function(ego_mf) str_wrap(ego_mf,width = 25))
cnetplot(ego_mf_dw, showCategory = 5)






#------------------------------------------------------------PID27--------------------------------------------------------------#
#-----------------------------------------------GO enrichment analysis-----------------------------------------------------#
#----------Extract the data table for GO enrichment analysis;
Ecs_PID27a1_human_overlap_up <- dplyr::filter(Ecs_PID27a1_exp_diffs_up, Ecs_PID27a1_exp_diffs_up$SYMBOL %in% PID27_overlap_up_genes)


library(org.Mm.eg.db)
library(clusterProfiler) 
gene_ids <- bitr(Ecs_PID27a1_human_overlap_up$gene, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Mm.eg.db")
# I used ENSEMBL gene id for the matrix, and annotate with gene symbols.
#res_anno<- as.data.frame(res_3v_dw_up) %>% left_join(gene_symbols, c("ENSEMBL" = "gene_id"))
Ecs_PID27a1_human_overlap_up <- left_join(Ecs_PID27a1_human_overlap_up, gene_ids, c("gene" = "SYMBOL"))

########ID choose 0.001; 4 fold
res_de_bAVM_upx_id = data.frame(ID=Ecs_PID27a1_human_overlap_up$ENSEMBL) 


#------------------------------------------up regulated genes enrichment analysis-----------------------------------------------#
library(DOSE) 
library(annotate)
eg = bitr(de_bAVM_x_id$ID, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
############Cellular function;
ego_cc<-enrichGO(gene = res_de_bAVM_upx_id$ID,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'ENSEMBL',
                 ont        = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.01,
                 readable      = TRUE) #TRUE 
ego_cc_results<-summary(ego_cc)
barplot(ego_cc,showCategory = 15, title="The GO_CC of DEGs(up_0-01_4fc)")
dotplot(ego_cc, showCategory = 15, title="The GO_CC of DEGs(up_0-01_4fc)")
#######dotplot(ego_cc, showCategory = 20, split=".sign", title="The GO_CC enrichment analysis of all DEGs") + facet_grid(.~.sign)
#  scale_size(range=c(2, 12))+
#  scale_x_discrete(labels=function(ego_cc) str_wrap(ego_cc,width = 25))
library(topGO)
library(GSEABase)
##########################################
plotGOgraph(ego_cc)
##igraph format of plot GOgraph
goplot(ego_cc)
##go term overlapping relationship
emapplot(ego_cc, showCategory = 30)
##gene and go term
cnetplot(ego_cc, showCategory = 5)
##enrichmap
library(enrichplot)
require(doseplot)
library(igraph)
###enrichMap(ego_cc)
####molecular function


#####################biological process;
ego_bp<-enrichGO(gene       = res_de_bAVM_upx_id$ID,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'ENSEMBL',
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05,
                 readable      = TRUE) #TRUE shows SYMBOL，FALSE shows original ID)
ego_bp_results<-summary(ego_bp)
barplot(ego_bp,showCategory = 15,title="The GO_BP of DEGs(up_0-01_4fc)")
##scale_size(range=c(2, 12))+
##scale_x_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 40))
dotplot(ego_bp, showCategory = 15, title="The GO_BP of DEGs(up_0-01_4fc)")
#+  scale_size(range=c(2, 12))


######################molecular function;
ego_mf_up<- enrichGO(gene = res_de_bAVM_upx_id$ID,
                     OrgDb      = org.Mm.eg.db,
                     keyType    = 'ENSEMBL',
                     ont        = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable      = TRUE) #TRUE shows SYMBOL，FALSE shows original ID)
ego_mf_up_results<-summary(ego_mf_up)
barplot(ego_mf_up,showCategory = 15,title="The GO_MF of DEGs(up_0-05_2fc)")
dotplot(ego_mf_up, showCategory = 15, title="The GO_MF of DEGs(up_0-05_2fc)")
#  scale_size(range=c(2, 12))+
#  scale_x_discrete(labels=function(ego_mf) str_wrap(ego_mf,width = 25))
cnetplot(ego_mf_dw, showCategory = 5)








