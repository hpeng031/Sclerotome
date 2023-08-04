##scRNAseq project for tenocyte regeneration paper
##nkx3.1:Gal4; UAS:NTR-mCherry  + FACS for mCherry+
##52 hpf trunk only
##library created with 10X single cell 3' v3 chemistry


#loading necessary packages
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)
library(ggplot2)
library(ggprism)

#Reading in 10X data and creating Seurat object
counts_52hpf <- Read10X_h5(filename = "raw/filtered_feature_bc_matrix_52hpf.h5", use.names = TRUE, unique.features = TRUE)
unfiltered_52hpf <- CreateSeuratObject(counts = counts_52hpf, min.cells = 3, min.features = 200)
unfiltered_52hpf@meta.data$sample <- "52hpf"

##Adding QC data
unfiltered_52hpf[["percent.mt"]] <- PercentageFeatureSet(object = unfiltered_52hpf, pattern = "^mt-")
unfiltered_52hpf[["percent.rb"]] <- PercentageFeatureSet(object = unfiltered_52hpf, pattern = "^rp")
unfiltered_52hpf[["GenesPerUMI"]] <- unfiltered_52hpf$nFeature_RNA/unfiltered_52hpf$nCount_RNA
unfiltered_52hpf[["log10GenesPerUMI"]] <- log10(unfiltered_52hpf$nFeature_RNA) / log10(unfiltered_52hpf$nCount_RNA)
saveRDS(unfiltered_52hpf, file = "raw/unfiltered_52hpf.rds")

##filtering
filtered_52hpf <- subset(unfiltered_52hpf, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 10 & log10GenesPerUMI > 0.8)
#filtering genes not expressed, keeping genes expressed in >=10 cells
counts <- GetAssayData(object = filtered_52hpf, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_52hpf <- CreateSeuratObject(filtered_counts, meta.data = filtered_52hpf@meta.data)
saveRDS(filtered_52hpf, file = "raw/filtered_52hpf.rds")

Idents(unfiltered_52hpf)<-"52hpf"

#Visualize QC metrics
VlnPlot(unfiltered_52hpf, features = "nFeature_RNA")+
theme_prism()+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", color = "black", size = 11))+
  theme(axis.text = element_text(color = "black", size = 10))+
  theme(legend.position = "none")+
  guides(y=guide_prism_minor())+
  geom_hline(yintercept = 300, linetype = "dashed", size = 0.5)+
  geom_hline(yintercept = 5000, linetype = "dashed", size = 0.5)+
  ylab("Number of Genes")+
  ggtitle(element_blank())

VlnPlot(unfiltered_52hpf, features = "nCount_RNA")+
  theme_prism()+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", color = "black", size = 11))+
  theme(axis.text = element_text(color = "black", size = 10))+
  theme(legend.position = "none")+
  guides(y=guide_prism_minor())+
  ylab("Number of Transcripts")+
  ggtitle(element_blank())

VlnPlot(unfiltered_52hpf, features = "percent.mt")+
  theme_prism()+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", color = "black", size = 11))+
  theme(axis.text = element_text(color = "black", size = 10))+
  theme(legend.position = "none")+
  guides(y=guide_prism_minor())+
  geom_hline(yintercept = 10, linetype = "dashed", size = 0.5)+
  ylab("Proportion of Mitochondrial Genes")+
  ggtitle(element_blank())


#Visualize nUMI vs nGenes percent.mt
unfiltered_metadata <- unfiltered_52hpf@meta.data
filtered_metadata <- filtered_52hpf@meta.data

p1 <- unfiltered_metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm, formula = y~x, size = 1) +
  scale_x_log10(limits=c(300,100000)) + 
  scale_y_log10(limits=c(200,10000)) + 
  theme_prism() +
  theme(axis.title = element_text(face = "bold", color = "black", size = 11))+
  theme(axis.text = element_text(color = "black", size = 10))+
  theme(legend.text = element_text(color = "black", face = "bold", size = 10))+
  theme(legend.title = element_text(size = 10))+
  ylab("Number of Genes")+
  xlab("Number of Transcripts")+
  labs(color = "% Mitochondrial")+
  geom_hline(yintercept = 300, linetype = "dashed", size = 0.5) +
  geom_hline (yintercept = 5000, linetype = "dashed", size = 0.5)

p2 <- filtered_metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm, formula = y~x, size = 1) +
  scale_x_log10(limits=c(300,100000)) + 
  scale_y_log10(limits=c(200,10000)) + 
  theme_prism() +
  theme(axis.title = element_text(face = "bold", color = "black", size = 11))+
  theme(axis.text = element_text(color = "black", size = 10))+
  theme(legend.text = element_text(color = "black", face = "bold", size = 10))+
  theme(legend.title = element_text(size = 10))+
  ylab("Number of Genes")+
  xlab("Number of Transcripts")+
  labs(color = "% Mitochondrial")+
  geom_hline(yintercept = 300, linetype = "dashed", size = 0.5) +
  geom_hline (yintercept = 5000, linetype = "dashed", size = 0.5)

p <- p1 +p2
p

##R squared=0.82


#look at cell complexity  
unfiltered_metadata %>%
   ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
   geom_density(alpha = 0.2) +
   theme_prism() +
  theme(axis.title = element_text(face = "bold", color = "black", size = 11))+
  theme(axis.text = element_text(color = "black", size = 10))+
  theme(legend.position = "none")+
  xlab("Log10 Genes/Transcripts")+
  ylab("Cell Density")+
   geom_vline(xintercept = 0.8, linetype = "dashed", size = 0.5)


##To compare cell nos. post filtering
ncol(unfiltered_52hpf@assays[["RNA"]]@counts)
#3167
ncol(filtered_52hpf@assays[["RNA"]]@counts)
#2705

#############################################################################
##Dimension reduction
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(ggprism)
library(AnnotationHub)

##Cell cycle scoring
g2m_genes <- readRDS("C:/Users/arshe/OneDrive/Desktop/Working data/scRNAseq/nkx3.1_48hpf/raw_data/g2mgenesannotated.rds")
s_genes <- readRDS("C:/Users/arshe/OneDrive/Desktop/Working data/scRNAseq/nkx3.1_48hpf/raw_data/sgenesannotated.rds")
filtered_52hpf <- CellCycleScoring(filtered_52hpf, 
                                    g2m.features = g2m_genes, 
                                    s.features = s_genes)
filtered_52hpf$CC.difference <- filtered_52hpf$S.Score - filtered_52hpf$G2M.Score
filtered_52hpf <- SCTransform(filtered_52hpf, vars.to.regress = "percent.mt")
filtered_52hpf <- RunPCA(filtered_52hpf)

ElbowPlot(filtered_52hpf, ndims = 50)+
  theme_prism()+
  theme(axis.title = element_text(face = "bold", color = "black", size = 11))+
  theme(axis.text = element_text(color = "black", size = 10))+
  geom_vline(xintercept = 40, linetype = "dashed", size = 0.5)


filtered_52hpf <- RunUMAP(filtered_52hpf, dims = 1:40)
DimPlot(filtered_52hpf)
saveRDS(filtered_52hpf, file = "IP/SCTransformed_52hpf.rds")

#QC 
FeaturePlot(SCTransformed_52hpf, features = c("percent.mt", "percent.rb", "S.Score", "G2M.Score", "CC.difference"))& theme_prism() +
  theme(axis.title = element_blank())+
  theme(axis.text = element_text(color = "black", size = 10))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 11))

############################################################################
##Clustering
DefaultAssay(SCTransformed_52hpf)<-"SCT"
SCTransformed_52hpf <- FindNeighbors(object = SCTransformed_52hpf, dims = 1:40)
SCTransformed_52hpf <- FindClusters(object = SCTransformed_52hpf,
                                    resolution = c(0.6, 0.8, 1.0, 1.4))
Idents(SCTransformed_52hpf)<- "SCT_snn_res.0.8"
DimPlot(SCTransformed_52hpf)

##Cluster annotation
#SCT_snn_res.0.8
markers <- FindAllMarkers(SCTransformed_52hpf, only.pos = TRUE,
                          logfc.threshold = 0.25)

# Add ensembl gene descriptions to markers table. 
library(tibble)
ah <- AnnotationHub()
# Access the Ensembl database for organism
ahDb <- query(ah, pattern = c("Danio rerio", "EnsDb"), ignore.case = TRUE)
# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)
# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

markers <- markers %>% 
  rownames_to_column(var="genes") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description","gene_biotype")]),
            by = c("gene" = "gene_name"))
View(markers)
write.csv(markers, file = "IP/markers_52hpf.csv")
write.csv(annotations, file = "raw/Daniorerio_annotations.csv")

###Examine markers
DefaultAssay(SCTransformed_52hpf) <- "RNA"
SCTransformed_52hpf <- NormalizeData(SCTransformed_52hpf, verbose = FALSE)

#general fibroblast
FeaturePlot(SCTransformed_52hpf, 
            reduction = "umap", 
            features = c( "col1a2", "col5a1", "pdgfra", "dcn"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE) & theme_prism() +
  theme(axis.title = element_blank())
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 10))

VlnPlot(SCTransformed_52hpf, features = c("mCherry", "col1a2", "col5a1", "pdgfra", "dcn", "lum"))& theme_prism()+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", size = 10, color = "black"))+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 11))+
  theme(legend.position = "none")
  
#general sclerotome
FeaturePlot(SCTransformed_52hpf, 
            reduction = "umap", 
            features = c("postnb", "twist1b", "tgfbi", "matn4"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = FALSE) & theme_prism() +
  theme(axis.title = element_blank())+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 10))

Sclerotome_top30 <- c("nkx3-1", "mCherry","pax9", "nkx3.2", "pax1a", "twist1b", "vcanb", "foxc1b", "tnn", "tgfbi", "col1a2", "matn4", "pmp22a", "fmoda", "wnt11r", "col1a1b", "mfap2", "tpm4a", "vwde", "cd82a", "col5a1", "postnb", "ca6", "dcn", "nkx3.2","cthrc1a", "crabp1b", "col1a1a", "hapln1a", "postnb", "scxa", "prelp", "hgd", "hmcn4")
DoHeatmap(SCTransformed_52hpf_clean, features = Sclerotome_top30)

VlnPlot(SCTransformed_52hpf, features = Sclerotome_miller)& theme_prism()+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", size = 10, color = "black"))+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 11))+
  theme(legend.position = "none")

VlnPlot(SCTransformed_52hpf, features = c("nkx3-1", "pax9", "nkx3.2", "pax1a", "twist1b", "vcanb", "foxc1b"))& theme_prism()+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", size = 10, color = "black"))+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 11))+
  theme(legend.position = "none")

#tenocytes
FeaturePlot(SCTransformed_52hpf, 
            reduction = "umap", 
            features = c("scxa", "tnmd", "prelp", "cilp", "tcf15"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE) & theme_prism() +
  theme(axis.title = element_blank())+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 10))

VlnPlot(SCTransformed_52hpf, features = c("scxa", "tnmd", "prelp", "cilp", "tcf15"))& theme_prism()+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", size = 10, color = "black"))+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 11))+
  theme(legend.position = "none")

#Notochord Associated
FeaturePlot(SCTransformed_52hpf, 
            reduction = "umap", 
            features = c("ptch1", "ptch2", "nkx3.2","pax1a"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE) & theme_prism() +
  theme(axis.title = element_blank())+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 10))

VlnPlot(SCTransformed_52hpf, features = c("ptch1", "ptch2", "nkx3.2","pax1a"))& theme_prism()+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", size = 10, color = "black"))+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 11))+
  theme(legend.position = "none")

#fin
FeaturePlot(SCTransformed_52hpf, 
            reduction = "umap", 
            features = c("fbln1", "hmcn2", "hmcn1", "and1", "hgd"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE) & theme_prism() +
  theme(axis.title = element_blank())+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 10))

VlnPlot(SCTransformed_52hpf, features = c("fbln1", "hmcn2", "hmcn1", "and2", "hgd"))& theme_prism()+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", size = 10, color = "black"))+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 11))+
  theme(legend.position = "none")

#muscle progenitor
FeaturePlot(SCTransformed_52hpf, 
            reduction = "umap", 
            features = c("fras1", "bmp4", "bmp2b", "bmp6"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE) & theme_prism() +
  theme(axis.title = element_blank())+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 10))

VlnPlot(SCTransformed_52hpf, features = c("fras1", "bmp4", "bmp2b", "bmp6"))& theme_prism()+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", size = 10, color = "black"))+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 10))+
  theme(legend.position = "none")

#endothelial cells
FeaturePlot(SCTransformed_52hpf, 
            reduction = "umap", 
            features = c("kdrl", "fli1a", "etv2","flt4", "ephb4a", "ephb2a"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE) & theme_prism() +
  theme(axis.title = element_blank())+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 10))

VlnPlot(SCTransformed_52hpf, features = c("kdrl", "fli1a", "etv2","flt4", "ephb4a", "ephb2a"))& theme_prism()+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", size = 10, color = "black"))+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 10))+
  theme(legend.position = "none")


SCTransformed_52hpf <- RenameIdents(object = SCTransformed_52hpf, 
                                    "0" = "Scl6",
                                    "1" = "Stressed",
                                    "2" = "Scl3",
                                    "3" = "Cycling",
                                    "4" = "Scl4",
                                    "5" = "Scl1",
                                    "6" = "Scl2",
                                    "7" = "Myotome",
                                    "8" = "Stressed 2",
                                    "9" = "Scl5",
                                    "10" = "Periderm",
                                    "11" = "Neuron",
                                    "12" = "Muscle",
                                    "13" = "Neural Crest",
                                    "14" = "Macrophage",
                                    "15" = "Endothelium")
DimPlot(SCTransformed_52hpf, label = TRUE, pt.size = 1, repel = TRUE, label.size = 4) + theme_prism()+
  theme(axis.title = element_text(face = "bold", size = 11))+
  theme(axis.text = element_text(color = "black", size = 10))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 11))+
  theme(legend.text = element_text(face = "bold", size = 10))+
  xlab("UMAP 1") +
  ylab("UMAP 2")


# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(SCTransformed_52hpf, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
View(n_cells)

# Fin Mesenchymal Cells/scl6 = 410 
# Stressed = 303
# Fibroblast 1/Scl3 = 258
# Cycling Fibroblasts = 256,
# Fibroblast 2/Scl4 = 255
# Notochord Associated Fibroblasts/Scl1 = 250
# Fibroblast 3//Scl2 = 229
# Myotome = 207
# Stressed 2 = 140
# Tenocytes/scl5 = 98
# Periderm = 98
# Neuron = 79
# Muscle = 47
# Neural Crest = 35
# Macrophage = 23
# Endothelium = 17

SCTransformed_52hpf@meta.data$identity <- Idents(SCTransformed_52hpf)
saveRDS(SCTransformed_52hpf, file = "IP/SCTransformed_52hpf.rds")

##Visualize top markers
top10 <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(SCTransformed_52hpf_clean, features=top10$gene)

#################################################################
#Other general analysis
library(readxl)
core_matrisome_genes <- read_excel("raw/core matrisome genes.xls")
matrisome_assoc_genes <- read_excel("raw/matrisome assoc genes.xls")
matrisome <- core_matrisome_genes$Gene
matrisome_assoc <- matrisome_assoc_genes$Gene
expressed <- rownames(SCTransformed_52hpf@assays$RNA@counts)
expressed_matrisome <- matrisome[matrisome %in% expressed]
expressed_matrisome_assoc <- matrisome_assoc[matrisome_assoc%in% expressed]
SCTransformed_52hpf[["percent.matrisome"]] <- PercentageFeatureSet(object = SCTransformed_52hpf, features = expressed_matrisome)
SCTransformed_52hpf[["percent.matrisome_assoc"]] <- PercentageFeatureSet(object = SCTransformed_52hpf, features = expressed_matrisome_assoc)

##expressed matrisome = 192/331
##expressed matrisome assoc = 217/364

SCTransformed_clean <- subset(SCTransformed_52hpf, idents = c("Stressed", "Stressed 2"), invert = TRUE)
DefaultAssay(SCTransformed_52hpf_clean)<- "RNA"
FeaturePlot(SCTransformed_52hpf_clean, features = "percent.matrisome", min.cutoff = 0, max.cutoff = 25)& theme_prism() +
  theme(axis.title = element_blank())+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 11))

VlnPlot(SCTransformed_52hpf_clean, features = "percent.matrisome")+ theme_prism()+
  theme(axis.title.x = element_blank(), axis.title.y.left = element_text(face = "bold", size = 10, color = "black"))+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(axis.text.x = element_text(angle = 90))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 11))+
  theme(legend.position = "none")+
  ylab("Proportion Matrisome Genes (%)")
  

FeaturePlot(SCTransformed_52hpf_clean, features = "percent.matrisome_assoc", min.cutoff = 0, max.cutoff = 25)& theme_prism() +
  theme(axis.title = element_blank())+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 11))

VlnPlot(SCTransformed_52hpf_clean, features = "percent.matrisome_assoc") + theme_prism()+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", size = 10, color = "black"))+
  theme(axis.text = element_text(color = "black", size = 8))+
  theme(axis.text.x = element_text(angle = 90))+
  theme(plot.title = element_blank())+
  theme(legend.position = "none")+
  ylab("Proportion of Matrisome-Associated Genes (%)")


##########################################################
##Subsetting

Fibroblasts_52hpf <- subset(SCTransformed_52hpf, idents = c("Scl1", "Scl2", "Scl3", "Scl4", "Scl5", "Scl6", "Cycling"))

##Subsetted total of 1756 cells

Fibroblasts_52hpf <- SCTransform(Fibroblasts_52hpf, vars.to.regress = "percent.mt")
Fibroblasts_52hpf <- RunPCA(Fibroblasts_52hpf)
ElbowPlot(Fibroblasts_52hpf, ndims = 50)
Fibroblasts_52hpf <- RunUMAP(Fibroblasts_52hpf, dims= 1:40)

DimPlot(Fibroblasts_52hpf, reduction = "umap", label = FALSE, pt.size = 1, cols = c("darkseagreen4", "coral", "violet", "darkslateblue", "grey", "deepskyblue", "deeppink2"))+ theme_prism()+
  theme(axis.title = element_text(face = "bold", size = 11))+
  theme(axis.text = element_text(color = "black", size = 10))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 11))+
  theme(legend.text = element_text(face = "bold", size = 10))+
  xlab("UMAP 1") +
  ylab("UMAP 2")
saveRDS(Fibroblasts_52hpf, file = "results/Fibroblasts_52hpf.rds")

fib_markers <- FindAllMarkers(Fibroblasts_52hpf, only.pos = TRUE,
                          logfc.threshold = 0.25)
fib_topmarkers <- FindAllMarkers(Fibroblasts_52hpf, only.pos = TRUE,
                              logfc.threshold = 0.25, min.diff.pct = 0.3)

Fibroblasts_52hpf <- RenameIdents(object = Fibroblasts_52hpf, 
                                    "Scl1" = "Notochord Associated Fibroblast",
                                    "Scl2" = "Dorsal Fibroblast",
                                    "Scl3" = "Interstitial/Perivascular Fibroblast",
                                    "Scl4" = "Stromal Reticular Cell",
                                    "Scl5" = "Tenocyte",
                                    "Scl6" = "Fin Mesenchymal Cell")

levels(Fibroblasts_52hpf) <- c("Interstitial/Perivascular Fibroblast", "Stromal Reticular Cell", "Dorsal Fibroblast", "Fin Mesenchymal Cell", "Cycling", "Notochord Associated Fibroblast", "Tenocyte")
DimPlot(Fibroblasts_52hpf, reduction = "umap", label = FALSE, pt.size = 1, cols = c("darkseagreen4", "coral", "violet", "darkslateblue", "grey", "deepskyblue", "deeppink2"))+ theme_prism()+
  theme(axis.title = element_text(face = "bold", size = 11))+
  theme(axis.text = element_text(color = "black", size = 10))+
  theme(plot.title = element_text(face = "bold", color = "black", size = 11))+
  theme(legend.text = element_text(face = "bold", size = 10))+
  xlab("UMAP 1") +
  ylab("UMAP 2")
saveRDS(Fibroblasts_52hpf, file = "results/Fibroblasts_52hpf.rds")

top_fin <- c ("hgd", "fbln1", "hpdb")
top_tenocyte <- c("prp", "cilp", "prelp")
top_notochord <- c("nkx3.2", "crispld1b", "ctgfb")
top_src <- c("clec19a", "sfrp1a", "ccl25b")
top_dorsal <- c("wnt11r", "zic1", "ca6")
top_ifpf <- c("hapln1a", "cdh11", "pcdh18b")
top_cycling <- c("mki67", "pcna", "top2a")

top <- c(top_fin, top_tenocyte, top_src, top_notochord, top_dorsal, top_ifpf, top_cycling)
levels(Fibroblasts_52hpf) <- c("Cycling", "Interstitial/Perivascular Fibroblast", "Dorsal Fibroblast", "Notochord Associated Fibroblast", "Stromal Reticular Cell", "Tenocyte", "Fin Mesenchymal Cell")


DoHeatmap(Fibroblasts_52hpf, features = fib_markers$gene, label = FALSE)+ theme_prism()+
  theme(axis.text = element_blank())+
  theme(legend.title = element_text(face = "bold", size = 8))+
  theme(legend.text = element_text(face="bold", size =8))+
  theme(line = element_blank())

DotPlot(Fibroblasts_52hpf, features = top) + theme_prism()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_text(color = "black", size = 10))+
  theme(axis.text.x = element_text(angle = 90))+
  theme(plot.title = element_blank())+
  theme(legend.text = element_text(face = "bold", size = 10))+
  theme(legend.title = element_text(face = "bold", size = 10))



##########################################################
#RNA Velocity
library(SeuratDisk)
DefaultAssay(Fibroblasts_52hpf) <- "RNA"
Idents(Fibroblasts_52hpf)<- "identity"
##did not try to run this one. 
#SaveH5Seurat(fibroblasts_52hpf, filename = "results/fibroblasts.h5Seurat")
#Convert("results/fibroblasts.h5Seurat", dest = "h5ad")
##did not work well
#ibroblasts_52hpf <- as.loom(fibroblasts_52hpf, filename = "results/Fibroblasts.loom")

#####making input anndata file for use with scvelo ----worked well
# save metadata table:
setwd("./IP")
seurat_obj <- Fibroblasts_52hpf
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file='metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0(out_data_dir, 'counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_obj@reductions$pca@cell.embeddings, file='pca.csv', quote= F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)

setwd("..")
#############################################################
####Gene Ontology analysis - FishEnrichR

library(enrichR)
setEnrichrSite("FishEnrichr")

setwd("./results")
fibroblastsubtypes <- c("Tenocyte", "Cycling", "Notochord Associated Fibroblast", "Stromal Reticular Cell", "Dorsal Fibroblast", "Fin Mesenchymal Cell")
for (identity in fibroblastsubtypes) {
  go <- DEenrichRPlot(
    Fibroblasts_52hpf,
    ident.1 = identity,
    balanced = TRUE,
    logfc.threshold = 0.25,
    assay = "RNA",
    max.genes = 300,
    test.use = "wilcox",
    p.val.cutoff = 0.05,
    cols = NULL,
    enrich.database = "GO_Biological_Process_2018",
    num.pathway = 25,
    return.gene.list = TRUE)
  
  goplot <- DEenrichRPlot(
    Fibroblasts_52hpf,
    ident.1 = identity,
    balanced = TRUE,
    logfc.threshold = 0.25,
    assay = "RNA",
    max.genes = 500,
    test.use = "wilcox",
    p.val.cutoff = 0.05,
    cols = NULL,
    enrich.database = "GO_Biological_Process_2018",
    num.pathway = 25,
    return.gene.list = FALSE)
  
  print(goplot)
  
  write.csv(go, paste0("GO", identity, ".csv"))
  
}

setwd("..")