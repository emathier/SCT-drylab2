# install.packages("BiocManager")
# install.packages("remotes")
# install.packages("devtools")
# 
# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
# remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
# remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
# remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
# remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
# 
# install.packages("Signac")
# install.packages("Matrix")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("purrr")
# install.packages("doParallel")
# 
# BiocManager::install("AnnotationHub")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# BiocManager::install("TFBSTools")
# BiocManager::install("JASPAR2020")
# BiocManager::install("chromVAR")
# 
# devtools::install_github('quadbio/Pando')
# devtools::install_github("immunogenomics/presto")

library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg38)
library(presto)
library(ggplot2)
library(purrr)
library(TFBSTools)
library(JASPAR2020)
library(Pando)
library(doParallel)

setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/ETH/MSc Biotechnology/2025_SS/Single-Cell_Technologies/Drylab_Task2")

#-------------- Loading Data and create Seurat object -------------- 
#Load sequencing data
counts <- Read10X_h5("midbrain_organoid/filtered_feature_bc_matrix.h5")

#Retrieve gene annotation
ah <- AnnotationHub()
ensdbs <- query(ah, c("EnsDb.Hsapiens"))

ensdb_id <- ensdbs$ah_id[grep(paste0(" 98 EnsDb"), ensdbs$title)]
ensdb <- ensdbs[[ensdb_id]]

seqlevelsStyle(ensdb) <- "UCSC"
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotations) <- "hg38"

#Create seurat object
seurat <- CreateSeuratObject(counts = counts$`Gene Expression`,
                                 assay = "RNA",
                                 project = "midbrain_organoid")

seurat[['ATAC']] <- CreateChromatinAssay(counts = counts$`Peaks`,
                                         annotation = annotations,
                                         fragments = "midbrain_organoid/atac_fragments.tsv.gz",
                                         sep = c(":", "-"),
                                         genome = 'hg38')


#union the peak lists of different samples called by Cell Ranger ARC to 
#get the unified peak list and filter peaks
peaks <- Signac::reduce(unlist(as(c(seurat@assays$ATAC@ranges),
                          "GRangesList")))
peakwidths <- width(peaks)
peaks <- peaks[peakwidths < 10000 & peakwidths > 20]

counts_atac_merged <- FeatureMatrix(seurat@assays$ATAC@fragments,
                                    features = peaks,
                                    cells = colnames(seurat))
seurat[['ATAC']] <- CreateChromatinAssay(counts_atac_merged,
                                         fragments = seurat@assays$ATAC@fragments,
                                         annotation = seurat@assays$ATAC@annotation,
                                         sep = c(":","-"),
                                         genome = "hg38")

#Subset peaks to focus on standard chromosomes
standard_chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
idx_standard_chroms <- which(as.character(seqnames(granges(seurat[['ATAC']]))) %in% standard_chroms)
seurat[["ATAC"]] <- subset(seurat[["ATAC"]],
                           features = rownames(seurat[["ATAC"]])[idx_standard_chroms])
seqlevels(seurat[['ATAC']]@ranges) <- intersect(seqlevels(granges(seurat[['ATAC']])),
                                                unique(seqnames(granges(seurat[['ATAC']]))))

#-------------- Quality control -------------- 
seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")
seurat <- NucleosomeSignal(seurat, assay = "ATAC")
seurat <- TSSEnrichment(seurat, assay = "ATAC")

VlnPlot(seurat,
        features = c("nFeature_RNA",
                     "percent.mt",
                     "nFeature_ATAC",
                     "TSS.enrichment",
                     "nucleosome_signal"),
        ncol = 5,
        pt.size = 0)

##ggsave("Plots/1_Quality_Control.pdf", units = "mm",  dpi = "print", width = 300, height = 210)

seurat <- subset(seurat,
                 subset = nFeature_RNA > 1000 &
                   nFeature_RNA < 8500 &
                   percent.mt < 20 &
                   nFeature_ATAC > 1000 &
                   nFeature_ATAC < 30000 &
                   TSS.enrichment > 1 &
                   nucleosome_signal < 1
)

#-------------- Analysis of the RNA assay -------------- 
DefaultAssay(seurat) <- "RNA"

#data normalization, highly variable genes identification, data scaling, 
#principal component analysis (PCA), and then generate the UMAP embedding
seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  CellCycleScoring(s.features = cc.genes.updated.2019$s.genes,
                   g2m.features = cc.genes.updated.2019$g2m.genes) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:30, reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

p1 <- DimPlot(seurat, group.by = "orig.ident", reduction = "umap_rna") & NoAxes()
p2 <- FeaturePlot(seurat,
                  c("MAP2","OTX2","FOXA2","TH"),
                  reduction = "umap_rna") & NoAxes() & NoLegend()
#p1 | p2
p1

##ggsave("Plots/2_RNA_UMAP.pdf", units = "mm",  dpi = "print", width = 200, height = 150)

#--> don't need to do any data integration as we only have one replicate

#clustering and cluster marker identification
seurat <- FindNeighbors(seurat,
                        reduction = "umap_rna",
                        dims = 1:ncol(Embeddings(seurat,"umap_rna"))) %>%
  FindClusters(resolution = 0.2)

DE_cl_rna <- presto::wilcoxauc(seurat, "RNA_snn_res.0.2")
top_markers <- DE_cl_rna %>%
  dplyr::filter(logFC > log(1.2) &
           auc > 0.7 &
           padj < 0.01 &
           pct_in - pct_out > 30 &
           pct_out < 30) %>%
  group_by(group) %>%
  top_n(1, wt = auc)

p1 <- DimPlot(seurat,
              group.by="RNA_snn_res.0.2",
              reduction="umap_rna", label=T) & NoAxes() & NoLegend()
p2 <- FeaturePlot(seurat,
                  features = unique(top_markers$feature),
                  reduction="umap_rna",
                  order = T,
                  ncol=3) & NoAxes() & NoLegend()
(p1 | p2) + patchwork::plot_layout(widths = c(2,3))
p1

##ggsave("Plots/3_RNA_Clustering.pdf", units = "mm",  dpi = "print", width = 200, height = 150)


#-------------- Analysis of the ATAC assay -------------- 
#Feature selection
DefaultAssay(seurat) <- "ATAC"
seurat <- FindTopFeatures(seurat, min.cutoff = 50)

#Normalization
seurat <- RunTFIDF(seurat, method = 1)

#Linear dimensionality reduction
seurat <- RunSVD(seurat, n = 50)

#Non-linear dimensionality reduction with UMAP for visualization
p1 <- ElbowPlot(seurat, ndims = 30, reduction="lsi")
p2 <- DepthCor(seurat, n = 30)
p1 | p2

##ggsave("Plots/4_ElbowPlot_ATAC.pdf", units = "mm",  dpi = "print", width = 300, height = 150)


# First SVD component, which explains the most variance of the data, is usually 
# highly correlated with total counts of the ATAC fragments per cell, and 
# therefore represents more technical variance than biological variance. 
# In this case, we would exclude the first SVD component and start from the 
# second one.

seurat <- RunUMAP(seurat,
                  reduction = "lsi",
                  dims = 2:30,
                  reduction.name = "umap_atac",
                  reduction.key = "UMAPATAC_")
p1 <- DimPlot(seurat,
              group.by = "orig.ident",
              reduction = "umap_atac") & NoAxes()
p2 <- FeaturePlot(seurat,
                  c("MAP2","OTX2","FOXA2","TH"),
                  reduction = "umap_atac") & NoAxes() & NoLegend()
p1 | p2

p1

##ggsave("Plots/5_ATAC_UMAP.pdf", units = "mm",  dpi = "print", width = 200, height = 150)


#predict gene expression from ATAC seq data alone
gene_act <- GeneActivity(seurat)
seurat[['RNA_inferred']] <- CreateAssayObject(gene_act) %>% NormalizeData()

DefaultAssay(seurat) <- "RNA_inferred"
beach_colscheme <- colorRampPalette(c("#cdcdcd","#edf8b1","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#0c2c84"))
p3 <- FeaturePlot(seurat,
                  c("MAP2","OTX2","FOXA2","TH"),
                  reduction = "umap_atac",
                  cols = beach_colscheme(30)) & NoAxes() & NoLegend()
p1 | p2 | p3

#--> don't need to do any data integration as we only have one replicate


#-------------- Bi-modal integrative analysis of the RNA-ATAC scMultiome data -------------- 
#Weighted nearest neighbor analysis
seurat <- FindMultiModalNeighbors(seurat,
                                  reduction.list = list("pca", "lsi"),
                                  dims.list = list(1:ncol(Embeddings(seurat,"pca")),
                                                   1:ncol(Embeddings(seurat,"lsi"))),
                                  modality.weight.name = c("RNA.weight","ATAC.weight"),
                                  verbose = TRUE)

#pca = rna lin dim red
#lsi = atac lin dim red

seurat <- RunUMAP(seurat, nn.name = "weighted.nn", assay = "RNA")
seurat <- FindClusters(seurat, graph.name = "wsnn", resolution = 0.2)

#p1 <- UMAPPlot(seurat, group.by = "orig.ident") & NoAxes()
p2 <- UMAPPlot(seurat, group.by = "wsnn_res.0.2", label=T) & NoAxes()
p3 <- FeaturePlot(seurat,
                  c("MAP2","OTX2","FOXA2","TH"),
                  reduction = "umap") & NoAxes() & NoLegend()
p2
#ggsave("Plots/6_Int_RNA_ATAC_Clustering.pdf", units = "mm",  dpi = "print", width = 200, height = 150)

#Cell type gene marker identification based on RNA data
DefaultAssay(seurat) <- "RNA"
DE_ct <- wilcoxauc(seurat, "wsnn_res.0.2", seurat_assay = "RNA")
top_markers_ct <- DE_ct %>%
  dplyr::filter(abs(logFC) > log(1.2) &
                  padj < 0.01 &
                  auc > 0.65 &
                  pct_in - pct_out > 30 &
                  pct_out < 20) %>%
  group_by(group) %>%
  top_n(10, wt = auc)

#View(top_markers_ct)

#Cell type peak marker identification based on ATAC data
DefaultAssay(seurat) <- "ATAC"
DA_ct <- wilcoxauc(seurat, "wsnn_res.0.2", seurat_assay = "ATAC")
top_peaks_ct <- DA_ct %>%
  dplyr::filter(abs(logFC) > log(1.1) &
                  padj < 0.01 &
                  auc > 0.55) %>%
  group_by(group) %>%
  top_n(100, wt = auc)

marker_peak_ct <- top_peaks_ct %>% top_n(10, wt=auc)

#View(top_peaks_ct)


#cell type annotation

# neural_marker_genes = list(
#   'Neural Progenitor Cell' = c("SOX2","NES", "PAX6", 'SFRP1', 'PTN', 'FABP7', 'VIM', 'SLC1A3', 'CST3'),
#   'Astrocytes' = c('GFAP','AQP4','S100B', 'NDRG2', 'ID2', 'SPARCL1', 'APOE', 'ANXA2', 'S100A10'),
#   'Oligodendrocytes' = c('OLIG1','OLIG2', 'PDGFRA', 'IRX1', 'IRX2', 'CD9', 'PMP2'),
#   'Glutamatergic Neurons' = c('SLC17A7', 'GRIN1', 'GRIN2B', 'NEUROD2', 'NEUROD6', 'SLA', 'BHLHE22', 'SCG2', 'GAP43', 'SPINT2'),
#   'GABAergic Neurons' = c('SLC32A1','GAD1','GAD2', 'SLC12A5', 'DLX2', 'DLX5', 'DLX6-AS1', 'MEIS2', 'RUNX1T1', 'CALB2', 'CELF4', 'GOLGA8A', 'GRIA2'),
#   'G2M' = c('MKI67'),
#   'Proliferating Cells' = c('TOP2A', 'CENPF', 'MKI67', 'NUSAP1'))

neural_marker_genes = list(
  'NPC' = c('NES', 'SOX2'),
  'Neuroblast' = c('NHLH1'),
  'Neuron' = c('DCX', 'MAP2', 'MAPT'),
  #'Telencephalon' = c('FOXG1'),
  #'dTelen (cortical)' = c('EMX1', 'EMX2'),
  #'Cortical intermediate progenitor' = c('EOMES'),
  #'dTelen (cortical) glutameric neuron' = c('NEUROD6', 'SLC17A7'),
  #'Deeper layer cortical neuron' = c('BCL11B'),
  #'Upper layer cortical neuron' = c('SATB2'),
  #'Cajal-Retzius' = c('RELN'),
  #'Ganglionic eminence (GE)' = c('DLX2', 'DLX5'),
  #'Lateral ganglionic eminece (LGE) inhibitory neuron' = c('ISL1'),
  #'Medial ganglionic eminence (MGE) inhibitory neuron' = c('NKX2-1'),
  #'Diencephalon' = c('RSPO3', 'TCF7L2', 'LHX5', 'LHX9'),
  #'Midbrain' = c('OTX2', 'LMX1A', 'EN1'),
  #'Purkinje (cerebellar) progenitor' = c('CYP26A1'),
  'Purkinje' = c('TFAP2A', 'CA8'),
  #'Hindbrain (medulla/pons) and spinal cord' = c('HOXB2', 'HOXB5'),
  'Glutamatergic neuron' = c('SLC17A6'),
  'GABAergic neuron' = c('SLC32A1', 'GAD1', 'GAD2'),
  'Dopaminergic neuron' = c('TH'),
  #'Cholinergic neuron' = c('CHAT', 'ACHE'),
  #'Choroid plexus' = c('TTR'),
  'Astrocyte' = c('GFAP', 'AQP4', 'S100B'),
  #'Oligodendrocyte precursor' = c('OLIG1'),
  #'Oligodendrocyte' = c('MBP', 'SOX10'),
  'Neural crest derivative' = c('SOX10'),
  #'Microglia' = c('AIF1'),
  #'Endothelial cell' = c('CLDN5'),
  'Mesenchymal cell' = c('DCN'),
  'G2M' = c('MKI67'))

neuron_rna_exp <- DE_ct %>% dplyr::filter(feature %in% list_c(neural_marker_genes))

for(cat in names(neural_marker_genes)){
  neuron_rna_exp$neuron_marker[neuron_rna_exp$feature %in% neural_marker_genes[[cat]]] <- cat
}

ggplot(neuron_rna_exp, aes(x = feature, y = group, col=avgExpr, size = pct_in)) +
  geom_point() +
  facet_wrap(~neuron_marker, strip.position = "bottom", scales = "free_x") +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#ggsave("Plots/7_DotPlot_CellType_Markers.pdf", units = "mm",  dpi = "print", width = 300, height = 200)

celltypes <- c("Purkinje",
              "GABAergic neuron",
              "NPC",
              "Dopaminergic neuron",
              "NPC",
              "Mesenchymal cell",
              "Neuron",
              "GABAergic neuron",
              "NPC",
              "G2M",
              "NPC")

DE_ct$group <- celltypes[as.integer(DE_ct$group) + 1]

new_ident <- setNames(celltypes,levels(seurat$wsnn_res.0.2))

seurat$celltype <- new_ident[seurat$wsnn_res.0.2] %>% setNames(colnames(seurat))

DimPlot(seurat, reduction = "umap", group.by = "celltype") & NoAxes()
#ggsave("Plots/8_Int_Cluster_Annotation_1.pdf", units = "mm",  dpi = "print", width = 200, height = 150)

top_peaks_ct$celltype <- celltypes[as.integer(top_peaks_ct$group)+1]

#linking of genes to nearby peaks for the top 10 rna markers for each cluster
#(due to large number of linear regressions more would be too time-consuming)
seurat <- Signac::RegionStats(seurat,
                      genome = BSgenome.Hsapiens.UCSC.hg38)

seurat <- LinkPeaks(seurat,
                    peak.assay = "ATAC",
                    expression.assay = "RNA",
                    genes.use = top_markers_ct$feature)

#check chromatin accessibility around the marker TH for dopaminergic neurons
# and CA8 for Pukinje cells in different cell types
DefaultAssay(seurat) <- "ATAC"
p1 <- CoveragePlot(seurat,
                   region = "TH",
                   features = "TH",
                   extend.upstream = 1000,
                   extend.downstream = 1000)
p2 <- CoveragePlot(seurat,
                   region = "CA8",
                   features = "CA8",
                   extend.upstream = 1000,
                   extend.downstream = 1000)
patchwork::wrap_plots(p1, p2, ncol = 1)
##ggsave("Plots/9_CoverageMap_TH_for_Dop_C8_for_Purkinje.pdf", units = "mm",  dpi = "print", width = 250, height = 250)

p1 <- CoveragePlot(seurat,
                   region = "chr11-2171500-2172500",
                   extend.upstream = 1000,
                   extend.downstream = 1000)
p2 <- CoveragePlot(seurat,
                   region = "chr8-60280000-60281000",
                   extend.upstream = 1000,
                   extend.downstream = 1000)
patchwork::wrap_plots(p1, p2, ncol = 1)


#-------------- TF Binding motif enrichment analysis --------------

#Get and add motifs
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
df_pfm <- data.frame(t(sapply(pfm, function(x)
  c(id=x@ID, name=x@name, symbol=ifelse(!is.null(x@tags$symbol),x@tags$symbol,NA)))))

#Identify TF binding motifs in peak list
open_peaks <- AccessiblePeaks(seurat)
peaks_matched <- MatchRegionStats(meta.feature = seurat[['ATAC']]@meta.features[open_peaks, ],
                                  query.feature = seurat[['ATAC']]@meta.features[marker_peak_ct$feature, ],
                                  n = 50000)
seurat <- Signac::AddMotifs(seurat, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)

motif_enrichment_dopa <- FindMotifs(seurat,
                                    features = top_peaks_ct$feature[top_peaks_ct$celltype == "Dopaminergic neuron"],
                                    background = peaks_matched) %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol), df_pfm$id)[motif]) %>%
  mutate(padj = p.adjust(pvalue, method="BH"))

enriched_motif_dopa <- motif_enrichment_dopa %>%
  dplyr::filter(padj < 0.01 & fold.enrichment > 3) %>%
  top_n(-4, wt = padj)


motif_enrichment_purkinje <- FindMotifs(seurat,
                                     features = top_peaks_ct$feature[top_peaks_ct$celltype == "Purkinje"],
                                     background = peaks_matched) %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol), df_pfm$id)[motif]) %>%
  mutate(padj = p.adjust(pvalue, method="BH"))

enriched_motif_purkinje <- motif_enrichment_purkinje %>%
  dplyr::filter(padj < 0.01 & fold.enrichment > 3) %>%
  top_n(-4, wt = padj)

motif_enrichment_mesench <- FindMotifs(seurat,
                                        features = top_peaks_ct$feature[top_peaks_ct$celltype == "Mesenchymal cell"],
                                        background = peaks_matched) %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol), df_pfm$id)[motif]) %>%
  mutate(padj = p.adjust(pvalue, method="BH"))

enriched_motif_mesen <- motif_enrichment_mesench %>%
  dplyr::filter(padj < 0.01 & fold.enrichment > 3) %>%
  top_n(-4, wt = padj)

p1 <- MotifPlot(seurat, motifs = enriched_motif_dopa$motif[1:4], ncol=4)
p2 <- MotifPlot(seurat, motifs = enriched_motif_purkinje$motif[1:4], ncol=4)
p3 <- MotifPlot(seurat, motifs = enriched_motif_mesen$motif[1:1], ncol=4)
p1 / p2 / p3

##ggsave("Plots/10_TopPeak_MotifEnrich_Dopa_Purk_Mesen.pdf", units = "mm",  dpi = "print", width = 250, height = 200)

save.image(file="DataAnalysis_Output.RData") 
