.libPath("/cluster/home/lokaiser/miniconda3/lib/R/library")

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

setwd("~")

load(file="Post_ChromVar_Output.RData")

#-------------- Gene regulatory network reconstruction --------------
grn <- initiate_grn(seurat,
                    regions=phastConsElements20Mammals.UCSC.hg38,
                    rna_assay = "RNA", peak_assay = "ATAC")

grn <- find_motifs(grn,
                   pfm = Pando::motifs,
                   motif_tfs = Pando::motif2tf,
                   genome = BSgenome.Hsapiens.UCSC.hg38)

registerDoParallel(20)
grn <- infer_grn(grn,
                 peak_to_gene_method='Signac',
                 parallel = T,
                 tf_cor = 0.05,
                 method="glm",
                 family="gaussian",
                 scale=F,
                 verbose=T)

coef_grn <- coef(grn_object) %>%
  filter(padj < 0.01)

coef_grn

grn_object <- find_modules(grn_object,
                           p_thresh = 0.01)
regulons <- NetworkModules(grn_object)

positive_regulons <- regulons@features[['genes_pos']]
positive_regulons <- positive_regulons[lengths(positive_regulons) > 10]
negative_regulons <- regulons@features[['genes_neg']]
negative_regulons <- negative_regulons[lengths(negative_regulons) > 10]

DefaultAssay(seurat) <- "RNA"
mod_act_pos <- AddModuleScore(seurat,
                              features = positive_regulons,
                              name = "regulon_")@meta.data
mod_act_pos <- mod_act_pos[,grep("^regulon_", colnames(mod_act_pos))] %>%
  setNames(paste0(names(positive_regulons),"(+)"))
mod_act_neg <- AddModuleScore(seurat,
                              features = negative_regulons,
                              name = "regulon_")@meta.data
mod_act_neg <- mod_act_neg[,grep("^regulon_", colnames(mod_act_neg))] %>%
  setNames(paste0(names(negative_regulons),"(-)"))

seurat[['regulon']] <- CreateAssayObject(data = t(cbind(mod_act_pos, mod_act_neg)))

DefaultAssay(seurat) <- "RNA"
p1 <- FeaturePlot(seurat,
                  top_tfs_ct$feature,
                  reduction = "umap",
                  cols = beach_colscheme(30),
                  order = T,
                  ncol = 6) & NoAxes() & NoLegend()
DefaultAssay(seurat) <- "regulon"
p2 <- FeaturePlot(seurat,
                  features = c(intersect(paste0(top_tfs_ct$feature,"(+)"), rownames(seurat)),
                               intersect(paste0(top_tfs_ct$feature,"(-)"), rownames(seurat))),
                  reduction = "umap",
                  cols = bluered_colscheme(30),
                  order = T,
                  ncol = 6) & NoAxes() & NoLegend()
(p1 / p2) + patchwork::plot_layout(height = c(1,2))

ggsave("Plots/12_GRE_DE_and_TF_regulons.pdf", units = "mm",  dpi = "print", width = 300, height = 300)