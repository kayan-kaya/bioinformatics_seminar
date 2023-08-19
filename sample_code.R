library(Seurat)

# データをロード
load("demo_counts.rda")
load("demo_metadata.rda")
srt_obj <- CreateSeuratObject(
  counts = demo_counts, 
  meta.data = demo_metadata
)

Idents(srt_obj) <- "Human"

# 総遺伝子数と総カウント数をプロットして確認
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

# シーケンサーの精度による系統的誤差の補正
srt_obj <- SCTransform(srt_obj, verbose = F)
gc(verbose = F,reset = T)

# 主成分分析
srt_obj <- RunPCA(srt_obj, verbose = FALSE)
gc(verbose = F,reset = T)

# 主成分得点を２次元にプロットしてみよう
# 第１主成分vs第２主成分
DimPlot(
  srt_obj,
  reduction = "pca",
  dims = c(1, 2)
)

# 第４９主成分vs第５０主成分
DimPlot(
  srt_obj,
  reduction = "pca",
  dims = c(49, 50)
)

# Elbow Plot（各主成分の標準偏差を大きい順に並べた図）から下流の解析で使う次元数を選択
# 例：ndim <- 10
ElbowPlot(srt_obj, ndims = 40)


# クラスタリング
srt_obj <- FindNeighbors(srt_obj, dims = 1:ndim) # グラフの構築

srt_obj <- FindClusters(
  srt_obj,
  algorithm = 1,
  resolution = seq(0.1, 1, by = 0.1)
) # グラフ構造から細胞をグループ分け

# resolutionとはグループ分けの細かさ
# 0に近いほど大雑把に、１に近いほど細かく分けられる

# グループ分けの結果はmetadataに保存されている
View(srt_obj@meta.data)　

# resolution = 0.1 によるクラスタリング結果を適用
Idents(srt_obj) <- srt_obj@meta.data$SCT_snn_res.0.1

# UMAPで高次元座標を二次元座標に変換
srt_obj <- RunUMAP(srt_obj, dims = 1:13, verbose = FALSE)

# UMAPの座標をプロット
DimPlot(
  srt_obj,
  reduction = "umap",
  label = T
)

# クラスター0 vs その他で遺伝子の発現量を統計検定で比較
marker_0 <- FindMarkers(srt_obj, ident.1 = 0)

# 結果を確認
View(marker_0)

# どのクラスターがどの細胞種か考えてみよう！
