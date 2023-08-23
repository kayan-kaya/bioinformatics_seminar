library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

srt <- Read10X(data.dir = "./")

srt <- CreateSeuratObject(srt)

srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = "^Mt-")
srt[["percent.rb"]] <- PercentageFeatureSet(srt, pattern = "^Rp[sl]")

View(srt@meta.data)

VlnPlot(
  srt,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
  ncol = 4,
  pt.size=0,
  slot = "counts"
)

FeatureScatter(
  srt,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  slot = "counts"
)

FeatureScatter(
  srt,
  feature1 = "nCount_RNA",
  feature2 = "percent.rb",
  slot = "counts"
)

FeatureScatter(
  srt,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  slot = "counts"
)

FeatureScatter(
  srt,
  feature1 = "percent.rb",
  feature2 = "percent.mt",
  slot = "counts"
)

FeatureScatter(
  srt,
  feature1 = "percent.mt",
  feature2 = "nFeature_RNA",
  slot = "counts"
)

print(
  quantile(
    unlist(srt[["percent.mt"]]),
    seq(0.1, 0.9, by = 0.1)
  )
)

print(
  quantile(
    unlist(srt[["nCount_RNA"]]),
    seq(0.1, 0.9, by = 0.1)
  )
)

print(
  quantile(
    unlist(srt[["nFeature_RNA"]]),
    seq(0.1, 0.9, by = 0.1)
  )
)
  


srt@meta.data %>% 
  ggplot(aes(x = nCount_RNA)) + 
  geom_density(alpha = 0.3, fill="blue") + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = c(500, 1000, 1500))

srt@meta.data %>% 
  ggplot(aes(x = nFeature_RNA)) + 
  geom_density(alpha = 0.3, fill="blue") + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = c(500, 1000, 1500))

srt@meta.data %>% 
  ggplot(aes(x = percent.mt)) + 
  geom_density(alpha = 0.3, fill="blue") + 
  theme_classic() +
  geom_vline(xintercept = c(5, 10, 15, 20, 25))

srt@meta.data %>% 
  ggplot(aes(x = percent.rb)) + 
  geom_density(alpha = 0.3, fill="blue") + 
  theme_classic() +
  geom_vline(xintercept = c(5, 10, 15, 20, 25))

srt <- srt %>% 
  SCTransform(verbose = F) %>% 
  RunPCA(verbose = F) 

ElbowPlot(srt, ndims = 40)

srt <- srt %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(
    algorithm = 1,
    resolution = seq(0.1, 1, by = 0.1)
  ) %>% 
  RunUMAP(dims = 1:20, verbose = FALSE)

Idents(srt) <- srt@meta.data$SCT_snn_res.0.1

DimPlot(
  srt,
  reduction = "umap",
  label = T
)


srt@meta.data %>% 
  ggplot(aes(x = nCount_RNA, color = SCT_snn_res.0.1, fill =SCT_snn_res.0.1)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() +
  geom_vline(xintercept = c(500, 1000))+
  scale_fill_brewer(palette = "Paired") +
  scale_colour_brewer(palette = "Paired")

srt@meta.data %>% 
  ggplot(aes(x = nFeature_RNA, color = SCT_snn_res.0.1, fill =SCT_snn_res.0.1)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() +
  geom_vline(xintercept = c(500, 1000)) +
  scale_fill_brewer(palette = "Paired") +
  scale_colour_brewer(palette = "Paired")

srt@meta.data %>% 
  ggplot(aes(x = percent.mt, color = SCT_snn_res.0.1, fill =SCT_snn_res.0.1)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() +
  geom_vline(xintercept = c(5, 10, 15, 20, 25)) +
  scale_fill_brewer(palette = "Paired") +
  scale_colour_brewer(palette = "Paired")

VlnPlot(srt, features = "nCount_RNA")
VlnPlot(srt, features = "nFeature_RNA")
VlnPlot(srt, features = "percent.mt")

FeatureScatter(
  srt,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  slot = "counts"
)

FeatureScatter(
  srt,
  feature1 = "nCount_RNA",
  feature2 = "percent.rb",
  slot = "counts"
)

FeatureScatter(
  srt,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  slot = "counts"
)

FeatureScatter(
  srt,
  feature1 = "percent.rb",
  feature2 = "percent.mt",
  slot = "counts"
)

FeatureScatter(
  srt,
  feature1 = "percent.mt",
  feature2 = "nFeature_RNA",
  slot = "counts"
)

