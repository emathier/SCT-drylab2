.libPaths("/cluster/home/lokaiser/R/x86_64-pc-linux-gnu-library/4.4")

library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(TFBSTools)
library(JASPAR2020)
library(doParallel)
library(presto)

setwd("~")

load(file="DataAnalysis_Output.RData")

#-------------- ChromVAR: another motif enrichment analysis --------------
seurat <- Signac::RunChromVAR(seurat, genome = BSgenome.Hsapiens.UCSC.hg38)

DefaultAssay(seurat) <- "chromvar"
DA_motifs_ct <- wilcoxauc(seurat, group_by = "celltype", seurat_assay = "chromvar") %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol),
                           df_pfm$id)[feature])

enriched_motifs_ct <- DA_motifs_ct %>%
  filter(padj < 0.01 & auc > 0.7) %>%
  group_by(group)
top_motifs_ct <- top_n(enriched_motifs_ct, 3, wt=auc)

#plot 1 location

tfs <- read.table("midbrain_organoid/Homo_sapiens_TF.txt", sep="\t", header=T)

tf_motifs_ct <- enriched_motifs_ct %>%
  filter(symbol %in% tfs$Symbol)

marker_tfs_ct <- DE_ct %>%
  dplyr::filter(feature %in% tfs$Symbol &
           abs(logFC) > log(1.2) &
           padj < 0.01 &
           auc > 0.65 &
           pct_in - pct_out > 20) %>%
  inner_join(tf_motifs_ct,
             by = c("feature" = "symbol"),
             suffix = c("_tf","_motif")) %>%
  dplyr::filter(group_tf == group_motif)

top_tfs_ct <- group_by(marker_tfs_ct, group_tf) %>%
  top_n(3, wt = auc_motif)

save.image(file="Post_ChromVar_Output.RData") 

bluered_colscheme <- colorRampPalette(rev(c("#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4")))
FeaturePlot(seurat,
            features = top_motifs_ct$feature,
            cols = bluered_colscheme(30),
            reduction = "umap",
            ncol = 3) & NoAxes() & NoLegend()
ggsave("11_TopMofifs_ChromVAR.pdf", units = "mm",  dpi = "print", width = 500, height = 600)

beach_colscheme <- colorRampPalette(c("#cdcdcd","#edf8b1","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#0c2c84"))
DefaultAssay(seurat) <- "RNA"
p1 <- FeaturePlot(seurat,
                  top_tfs_ct$feature,
                  reduction = "umap",
                  order=T,
                  cols=beach_colscheme(30),
                  ncol=6) & NoAxes() & NoLegend()
DefaultAssay(seurat) <- "chromvar"
p2 <- FeaturePlot(seurat,
                  top_tfs_ct$feature_motif,
                  reduction = "umap",
                  order=T,
                  cols=bluered_colscheme(30),
                  ncol=6) & NoAxes() & NoLegend()
p1 / p2
ggsave("12_Enrich_TF_and_Celltype_ChromVAR.pdf", units = "mm",  dpi = "print", width = 600, height = 400)
