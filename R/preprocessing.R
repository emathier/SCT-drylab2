library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(hdf5r)
library(AnnotationHub)
library(ensembldb)
library(biovizBase)

# Read file
counts <- Read10X_h5("Data/raw/filtered_feature_bc_matrix.h5")


# Get annotations
ah <- AnnotationHub()
ensdbs <- query(ah, c("EnsDb.Hsapiens"))

ensdb_id <- ensdbs$ah_id[grep(paste0(" 98 EnsDb"), ensdbs$title)]
ensdb <- ensdbs[[ensdb_id]]
seqlevelsStyle(ensdb) <- "UCSC"
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotations) <- "hg38"


# Create seurat object
seurat <- CreateSeuratObject(counts = counts$`Gene Expression`,
                                 assay = "RNA",
                                 project = "BVO")
seurat[['ATAC']] <- CreateChromatinAssay(counts = counts$`Peaks`,
                                             annotation = annotations,
                                             fragments = "~/projects/SCT-drylab2/Data/raw/atac_fragments.tsv.gz",
                                             sep = c(":", "-"),
                                             genome = 'hg38')

# Save
saveRDS(seurat, "Data/seurat_obj_preprocessed.rds")

