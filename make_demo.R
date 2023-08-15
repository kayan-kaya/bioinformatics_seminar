library(scRNAseq)
library(Seurat)

# Load demo dataset from scRNAseq(Bioconductor)
demo_data <- BaronPancreasData('human')

demo_counts <- counts(demo_data)

demo_metadata <- as.data.frame(colData(demo_data))

demo_srt <- CreateSeuratObject(
  counts = demo_counts, 
  meta.data = demo_metadata
)


srt_obj <- demo_srt

Idents(srt_obj) <- "Human"

VlnPlot(
  srt_obj,
  features = c("nFeature_RNA", "nCount_RNA"),
  ncol = 2,
  pt.size=0,
  slot = "counts"
)

FeatureScatter(
  srt_obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  slot = "counts"
)

srt_obj <- SCTransform(srt_obj, verbose = F)

gc(verbose = F,reset = T)

srt_obj <- RunPCA(srt_obj, verbose = FALSE)
gc(verbose = F,reset = T)

ElbowPlot(srt_obj, ndims = 40)

srt_obj <- FindNeighbors(srt_obj, dims = 1:13)

srt_obj <- FindClusters(
  srt_obj,
  algorithm = 1,
  resolution = seq(0.1, 1, by = 0.1)
)

gc(verbose = F,reset = T)

srt_obj <- RunUMAP(srt_obj, dims = 1:13, verbose = FALSE)

DimPlot(
  srt_obj,
  reduction = "umap",
  label = T
) +
NoLegend()

Idents(srt_obj) <- srt_obj@meta.data$SCT_snn_res.0.1

res_01_label <- as.list(0:7)

markers <- lapply(
  X = res_01_label,
  FUN = function(x){FindMarkers(srt_obj, ident.1 = x)}
)

names(markers) <- 0:7

file_names <- paste0("./res01_cluster", 0:7, "markers.csv")

file_names <- as.list(file_names)

mapply(
  function(x, y){
    write.csv(x, file = y)
  },
  markers,
  file_names
)

save(
  demo_counts, 
  file = "demo_counts.rda",
  compress = "gzip",
  compression_level = 9
)

save(
  demo_metadata, 
  file = "demo_metadata.rda",
  compress = "gzip",
  compression_level = 9
)