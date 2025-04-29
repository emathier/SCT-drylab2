library(Seurat)
library(Signac)
library(presto)
library(Matrix)
library(dplyr)
library(ggplot2)
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(JASPAR2020)
library(Pando)
library(doParallel)

setwd('~/Projects/SingleCellTechnologies/drylab2')

### Mono-modal data analysis of the scMultiome data
## Step 1. Load the data
# load joint file (one sparse matrix per assay)
counts = Read10X_h5("filtered_feature_bc_matrix.h5")

# retrieve gene annotation
ah <- AnnotationHub()
ensdbs <- query(ah, c("EnsDb.Hsapiens"))
ensdb_id <- ensdbs$ah_id[grep(paste0(" 98 EnsDb"), ensdbs$title)]
ensdb <- ensdbs[[ensdb_id]]
seqlevelsStyle(ensdb) <- "UCSC"
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotations) <- "hg38"

# create Seurat object
seurat = CreateSeuratObject(counts = counts$`Gene Expression`, assay = "RNA")
seurat[['ATAC']] = CreateChromatinAssay(counts = counts$`Peaks`, annotation = annotations,
                                        fragments = "atac_fragments.tsv.gz", sep = c(":", "-"), genome = 'hg38')
# subset peaks to standard chromosomes 
standard_chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
idx_standard_chroms <- which(as.character(seqnames(granges(seurat[['ATAC']]))) %in% standard_chroms)
seurat[["ATAC"]] <- subset(seurat[["ATAC"]], features = rownames(seurat[["ATAC"]])[idx_standard_chroms]) # excluded ~50 features
seqlevels(seurat[['ATAC']]@ranges) <- intersect(seqlevels(granges(seurat[['ATAC']])),
                                                unique(seqnames(granges(seurat[['ATAC']])))) 


## Step 2. Quality control
seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")
seurat <- NucleosomeSignal(seurat, assay = "ATAC")
seurat <- TSSEnrichment(seurat, assay = "ATAC")
VlnPlot(seurat, features = c("nFeature_RNA", "percent.mt", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 5, pt.size = 0)
ggsave("quality_control.png", width = 12, height = 6)
# set manual cutoffs (cf. Data analysis methods in preprint & tutorial)
seurat <- subset(seurat,
                 subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 30 &
                   nFeature_ATAC > 1000 & nFeature_ATAC < 30000 & TSS.enrichment > 1 & nucleosome_signal < 2) # before 8656, now 7942 cells


## Step 3. Analysis on the RNA assay
DefaultAssay(seurat) <- "RNA"
seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  CellCycleScoring(s.features = cc.genes.updated.2019$s.genes,
                   g2m.features = cc.genes.updated.2019$g2m.genes) %>% # not described in paper, but potentially useful (no visible difference)
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20, reduction.name = "umap_rna", reduction.key = "UMAPRNA_") %>% 
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# gene expression cf. Figure 1g
FeaturePlot(seurat, c("SOX2", "NES", "NHLH1", "TFDP2", "STMN2", 
                      "MAPT", "GFAP", "AQP4", "SOX10", "DCN",
                      "MKI67", "FOXA2", "OTX2", "LMX1A", "GATA3",
                      "FZD10", "FGF8", "GBX2", "ZIC1", "LHX9",
                      "TLX3", "SLC17A6", "GAD1", "CA8", "TH",
                      "SLC6A5", "SLC6A3"), reduction = "umap_rna", ncol = 4) & NoAxes() & NoLegend() 
# ggsave("selected_features_UMAPRNA.png", width = 12, height = 18)

# plot and annotate clusters
p1 = DimPlot(seurat, reduction = "umap_rna", group.by = "RNA_snn_res.0.5") & NoAxes() 
# TODO: improve annotation!
annotation_RNA = setNames(c("Purkinje cell", "Dopaminergic neuron", "Glycinergic neuron", "Astrocyte/ Mesenchymal cells" , "?",
                            "NPC/ NE", "NPC/ NE", "Purkinje cell", "?", "NPC/ NE",
                            "?", "GABAergic neuron", "Astrocyte/ Mesenchymal cells", "NPC/ NE", "NPC/ NE"),
                          levels(seurat$RNA_snn_res.0.5))
seurat$annotation_RNA = annotation_RNA[seurat$RNA_snn_res.0.5] %>% setNames(colnames(seurat))
p2 = DimPlot(seurat, reduction = "umap_rna", group.by = "annotation_RNA") & NoAxes() 
p1 | p2
# ggsave("annotation_UMAPRNA.png", width = 14, height = 6)


## Step 4. Analysis on the ATAC assay
DefaultAssay(seurat) <- "ATAC"
seurat <- FindTopFeatures(seurat, min.cutoff = 50) # select all the peaks with fragment detected in at least 50 cells
seurat <- RunTFIDF(seurat, method = 1) # normalization
seurat <- RunSVD(seurat, n = 50) # dimension reduction

p1 <- ElbowPlot(seurat, ndims = 50, reduction = "lsi") # elbow ~ 20
p2 <- DepthCor(seurat, n = 50)
p1 | p2
# ggsave("evaluate_SVDATAC.png", width = 14, height = 6)

seurat <- RunUMAP(seurat, reduction = "lsi", dims = 2:20, reduction.name = "umap_atac", reduction.key = "UMAPATAC_")
p1 = DimPlot(seurat, reduction = "umap_atac", group.by = "RNA_snn_res.0.5") & NoAxes() # cells from from RNA clusters also cluster here!
p2 = DimPlot(seurat, reduction = "umap_atac", group.by = "annotation_RNA") & NoAxes() 
p1 | p2
# ggsave("annotationRNA_UMAPATAC.png", width = 14, height = 6)



### Bi-modal integrative analysis of the RNA-ATAC scMultiome data
## Step 1. Weighted nearest neighbor analysis
seurat <- FindMultiModalNeighbors(seurat, 
                                  reduction.list = list("pca", "lsi"),
                                  dims.list = list(1:ncol(Embeddings(seurat,"pca")), 1:ncol(Embeddings(seurat,"lsi"))),
                                  modality.weight.name = c("RNA.weight","ATAC.weight"),
                                  verbose = TRUE)
seurat <- RunUMAP(seurat, nn.name = "weighted.nn", assay = "RNA")
seurat <- FindClusters(seurat, graph.name = "wsnn", resolution = 0.5)
p1 <- UMAPPlot(seurat, group.by = "wsnn_res.0.2") & NoAxes()
p2 <- UMAPPlot(seurat, group.by = "RNA_snn_res.0.5") & NoAxes()
p3 <- UMAPPlot(seurat, group.by = "annotation_RNA") & NoAxes()
p1 | p2 | p3
# ggsave("annotationRNA_UMAPWNN.png", width = 20, height = 6)

## Step 2. Cell type gene/peak marker identification and visualization of the chromatin accessibility profiles
DefaultAssay(seurat) <- "RNA"
DE_ct <- wilcoxauc(seurat, "annotation_RNA", seurat_assay = "RNA")
top_markers_ct <- DE_ct %>%
  dplyr::filter(abs(logFC) > log(1.2) &
           padj < 0.01 &
           auc > 0.65 &
           pct_in - pct_out > 30 &
           pct_out < 20) %>%
  group_by(group) %>%
  top_n(10, wt = auc)
top_markers_ct

DefaultAssay(seurat) <- "ATAC"
DA_ct <- wilcoxauc(seurat, "annotation_RNA", seurat_assay = "ATAC")
top_peaks_ct <- DA_ct %>%
  dplyr::filter(abs(logFC) > log(1.1) &
           padj < 0.01 &
           auc > 0.55) %>%
  group_by(group) %>%
  top_n(100, wt = auc)
top_peaks_ct %>% top_n(5, wt=auc)

# link a gene with its nearby peaks, which are its potential cis-regulatory elements
seurat <- RegionStats(seurat, genome = BSgenome.Hsapiens.UCSC.hg38)
seurat <- LinkPeaks(seurat, peak.assay = "ATAC", expression.assay = "RNA", genes.use = top_markers_ct$feature)
# example: GATA3
CoveragePlot(seurat, region = "GATA3", features = "GATA3", group.by = "annotation_RNA", extend.upstream = 1000, extend.downstream = 1000)
# ggsave("profile_GATA3.png", width = 8, height = 6)


## Step 3. TF binding motif enrichment analysis
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
df_pfm <- data.frame(t(sapply(pfm, function(x)
  c(id=x@ID, name=x@name, symbol=ifelse(!is.null(x@tags$symbol),x@tags$symbol,NA)))))

seurat <- AddMotifs(seurat, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)

# identify enriched TF binding motifs in peak list, e.g. cell-type-specific accessible peaks
open_peaks <- AccessiblePeaks(seurat)
peaks_matched <- MatchRegionStats(meta.feature = seurat[['ATAC']]@meta.features[open_peaks, ],
                                  query.feature = seurat[['ATAC']]@meta.features[top_peaks_ct$feature, ],
                                  n = 50000)

# example: GABAergic neurons
motif_enrichment_GABA <- FindMotifs(seurat,
                                    features = top_peaks_ct$feature[top_peaks_ct$group == "GABAergic neuron"],
                                    background = peaks_matched) %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol), df_pfm$id)[motif]) %>%
  mutate(padj = p.adjust(pvalue, method = "BH"))
enriched_motif_GABA <- motif_enrichment_GABA %>%
  dplyr::filter(padj < 0.01 & fold.enrichment > 3) %>%
  top_n(-4, wt = padj)
p1 <- MotifPlot(seurat, motifs = enriched_motif_GABA$motif[1:4], ncol = 4)
# compare with expression of corresponding TFs
DefaultAssay(seurat) <- "RNA"
p2 <- FeaturePlot(seurat, c("GATA1", "TAL1", "CRX","DMBX1","PITX2"), reduction = "umap", order = T, ncol = 5) & NoAxes() & NoLegend()
p1 / p2
# ggsave("TFenrichmenet_GABAergic.png", width = 12, height = 6)


## Step 4. ChromVAR: motif enrichment analysis without peak list
DefaultAssay(seurat) <- "ATAC"
seurat <- RunChromVAR(seurat, genome = BSgenome.Hsapiens.UCSC.hg38) # run long

# differential activity analysis between cell types to identify motifs that show cell type specific enrichment
DefaultAssay(seurat) <- "chromvar"
DA_motifs_ct <- wilcoxauc(seurat, group_by = "annotation_RNA", seurat_assay = "chromvar") %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol),
                           df_pfm$id)[feature])
enriched_motifs_ct <- DA_motifs_ct %>%
  dplyr::filter(padj < 0.01 & auc > 0.7) %>%
  group_by(group)
top_motifs_ct <- top_n(enriched_motifs_ct, 3, wt = auc)

# list of TFs from the animalTFDB database (downloaded from Tutorial)
tfs <- read.table("Homo_sapiens_TF.txt", sep = "\t", header = T)
tf_motifs_ct <- enriched_motifs_ct %>%
  dplyr::filter(symbol %in% tfs$Symbol)
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

# visualize the enrichment profiles directly in feauture plots
beach_colscheme <- colorRampPalette(c("#cdcdcd","#edf8b1","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#0c2c84"))
DefaultAssay(seurat) <- "RNA"
p1 <- FeaturePlot(seurat, top_tfs_ct$feature,
                  reduction = "umap", order = T, cols = beach_colscheme(30), ncol = 5) & NoAxes() & NoLegend()
bluered_colscheme <- colorRampPalette(rev(c("#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4")))
DefaultAssay(seurat) <- "chromvar"
p2 <- FeaturePlot(seurat, top_tfs_ct$feature_motif,
                  reduction = "umap", order = T, cols = bluered_colscheme(30), ncol = 5) & NoAxes() & NoLegend()
p1 / p2
# ggsave("chromVAR_topMotifs.png", width = 20, height = 24)



### Gene regulatory network reconstruction
# initialize a Pando-compatible Seurat object, insert sequence conservation information
grn <- initiate_grn(seurat,
                    regions = phastConsElements20Mammals.UCSC.hg38,
                    rna_assay = "RNA", peak_assay = "ATAC")
# incorporate the TF binding motif information
grn <- find_motifs(grn,
                   pfm = Pando::motifs,
                   motif_tfs = Pando::motif2tf,
                   genome = BSgenome.Hsapiens.UCSC.hg38)
# parallelize 
registerDoParallel(8) # otherwise my Mac crashes!
grn <- infer_grn(grn, peak_to_gene_method = 'Signac', parallel = T,
                 tf_cor = 0.05, method = "glm", family = "gaussian", scale = F, verbose = T) # run veeery long >3h

# extract the TF-peak-target trios with significant p-values
coef_grn <- coef(grn) %>%
  filter(padj < 0.01) # 14781 TFs

# generate the list of regulons, which are genes that are co-regulated positively or negatively by the same TF
grn <- find_modules(grn, p_thresh = 0.01) # found 558 TF modules
regulons <- NetworkModules(grn)
positive_regulons <- regulons@features[['genes_pos']]
positive_regulons <- positive_regulons[lengths(positive_regulons) > 10]
negative_regulons <- regulons@features[['genes_neg']]
negative_regulons <- negative_regulons[lengths(negative_regulons) > 10]

# calculate the regulon activity score for each regulon in each cell
DefaultAssay(seurat) <- "RNA"
mod_act_pos <- AddModuleScore(seurat, features = positive_regulons, name = "regulon_")@meta.data
mod_act_pos <- mod_act_pos[,grep("^regulon_", colnames(mod_act_pos))] %>%
  setNames(paste0(names(positive_regulons),"(+)"))
mod_act_neg <- AddModuleScore(seurat, features = negative_regulons, name = "regulon_")@meta.data
mod_act_neg <- mod_act_neg[,grep("^regulon_", colnames(mod_act_neg))] %>%
  setNames(paste0(names(negative_regulons),"(-)"))
seurat[['regulon']] <- CreateAssayObject(data = t(cbind(mod_act_pos, mod_act_neg)))

# do the feature plots for the same TFs we found combining DE and chromVAR analysis, and their positive and negative regulons
p1 <- FeaturePlot(seurat, top_tfs_ct$feature,
                  reduction = "umap", cols = beach_colscheme(30), order = T, ncol = 5) & NoAxes() & NoLegend()
DefaultAssay(seurat) <- "regulon"
p2 <- FeaturePlot(seurat,
                  features = c(intersect(paste0(top_tfs_ct$feature,"(+)"), rownames(seurat)),
                               intersect(paste0(top_tfs_ct$feature,"(-)"), rownames(seurat))),
                  reduction = "umap", cols = bluered_colscheme(30), order = T, ncol = 5) & NoAxes() & NoLegend()
p1 / p2
# ggsave("Pando_regulons.png", width = 20, height = 24)

# save files for reproducibility
# saveRDS(seurat, "seurat.rds")
# saveRDS(grn, "grn.rds")

# TODO (optional, not covered in the tutorial): 
# graph-based analysis,
# centrality analysis to identify TFs which may trigger huge regulatory cascades,
# visualizing the TF cross-regulatory network
