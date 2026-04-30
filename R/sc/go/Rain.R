library(biomaRt)
library(dplyr)
library(tidyr)
library(plotly)
library(Seurat)


download.file(
  "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
  destfile = "pbmc3k.tar.gz"
)

untar("pbmc3k.tar.gz")

pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data)

save(pbmc, pbmc.data, file = "3K.RData")
#clean environment
load("3K.RData")

go <- read.delim("go_output.tsv", header = FALSE)

go_subset <- go[, c(1, 14)]
colnames(go_subset) <- c("protein", "go_terms")

go_subset$uniprot <- sub(".*\\|(.*?)\\|.*", "\\1", go_subset$protein)

#split multiple GO terms
go_long <- go_subset %>%
  separate_rows(go_terms, sep = "\\|")

#clean up
go_long <- go_long[go_long$go_terms != "-", ]
go_long <- go_long[go_long$go_terms != "", ]

colnames(go_long)[2] <- "go_id"

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  mart = mart
)

#clean up
mapping <- mapping[mapping$uniprotswissprot != "", ]
mapping <- mapping[mapping$hgnc_symbol != "", ]

#join protein → gene → domain
go_gene <- merge(
  go_long,
  mapping,
  by.x = "uniprot",
  by.y = "uniprotswissprot"
)

#count data
expr <- GetAssayData(pbmc, layer = "counts")

#keep only genes in PBMC
pbmc_genes <- rownames(expr)
go_gene <- go_gene[go_gene$hgnc_symbol %in% pbmc_genes, ]

#build gene → domain table
gene_go <- go_gene[, c("hgnc_symbol", "go_id")]

#build domain matrix
go_terms <- unique(na.omit(gene_go$go_id))

go_matrix <- matrix(
  0,
  nrow = length(go_terms),
  ncol = ncol(expr)
)

rownames(go_matrix) <- go_terms
colnames(go_matrix) <- colnames(expr)

#aggregate expression
for (g in go_terms) {
  genes_in_g <- gene_go$hgnc_symbol[gene_go$go_id == g]
  
  genes_in_g <- intersect(genes_in_g, rownames(expr))
  
  if (length(genes_in_g) > 0) {
    go_matrix[g, ] <- colSums(expr[genes_in_g, , drop = FALSE])
  }
}

#domains × cells
dim(go_matrix)

#remove truly empty cells
cs <- colSums(go_matrix)
#define valid cells
valid_cells <- names(cs[cs > 0])

#subset Seurat object
pbmc0 <- pbmc
pbmc0 <- subset(pbmc0, cells = valid_cells)
#force alignment
pbmc0 <- pbmc0[, valid_cells]

#overlay Seurat clustering
pbmc0 <- SCTransform(pbmc0, verbose = FALSE)
pbmc0 <- RunPCA(pbmc0, assay = "SCT")
pbmc0 <- FindNeighbors(pbmc0, dims = 1:30)
pbmc0 <- FindClusters(pbmc0)
pbmc0 <- RunUMAP(pbmc0, dims = 1:30)

table(pbmc0$seurat_clusters)

marker_list <- list(
  
  "CD4 T (naive)" = c("CD3D", "IL7R", "CCR7"),
  
  "CD8 T" = c("CD3D", "CD8A"),
  
  "B cell" = c("MS4A1", "CD79A"),
  
  "NK cell" = c("NKG7", "GNLY"),
  
  "CD14 Monocyte" = c("LYZ", "CD14", "S100A8"),
  
  "CD16 Monocyte" = c("FCGR3A", "MS4A7"),
  
  "Dendritic cell" = c("FCER1A", "CST3"),
  
  "Platelet" = c("PPBP", "PF4"),
  
  "Megakaryocyte" = c("PPBP", "ITGA2B"),
  
  "Plasma cell" = c("MZB1", "SDC1"),
  
  "CD4 T (memory)" = c("CD3D", "IL7R", "S100A4")
  
)

#ANCHOR FLAG
USE_ANCHOR <- FALSE

if (USE_ANCHOR) {
  #assign labels
  cluster_labels <- c(
    "0" = "CD4 T (naive)",
    "1" = "CD8 T",
    "2" = "B cell",
    "3" = "NK cell",
    "4" = "CD14 Monocyte",
    "5" = "CD16 Monocyte",
    "6" = "Dendritic cell",
    "7" = "Platelet",
    "8" = "Megakaryocyte",
    "9" = "Plasma cell",
    "10" = "CD4 T (memory)"
  )
  
  #build metadata
  meta_df <- data.frame(
    cell = colnames(pbmc0),
    cluster = as.factor(pbmc0$seurat_clusters),
    celltype = cluster_labels[as.character(pbmc0$seurat_clusters)]
  )
  
} else {
  #fallback
  meta_df <- data.frame(
    cell = colnames(pbmc0),
    cluster = as.factor(pbmc0$seurat_clusters),
    celltype = as.factor(pbmc0$seurat_clusters)
  )
}

#PCA coordinates
coords <- as.data.frame(Embeddings(pbmc0, "pca")[,1:2])
#force alignment
coords <- coords[colnames(pbmc0), ]
#check every cell
stopifnot(all(colnames(pbmc0) == rownames(coords)))

#UMAP coordinates
coords_umap <- as.data.frame(Embeddings(pbmc0, "umap")[,1:2])
#force alignment
coords_umap <- coords_umap[colnames(pbmc0), ]
#check every cell
stopifnot(all(colnames(pbmc0) == rownames(coords_umap)))

#UMAP FLAG
USE_UMAP <- FALSE
coords <- if (USE_UMAP) coords_umap else coords
axis_label <- if (USE_UMAP) "UMAP" else "PC"

#force order + identity
coords <- coords[valid_cells, ]
coords$cell <- rownames(coords)
colnames(coords)[1:2] <- c("x", "y")

#check exact cell order + identity
stopifnot(all(valid_cells == coords$cell))
stopifnot(all(valid_cells == meta_df$cell))

#build long format
df <- as.data.frame(go_matrix)
df$go_id <- rownames(df)

long <- pivot_longer(
  df,
  -go_id,
  names_to = "cell",
  values_to = "tpm"
)

#remove zeros
long <- long %>% filter(tpm > 0)

#attach PCA coordinates
long <- long %>%
  inner_join(meta_df, by = "cell") %>%
  inner_join(coords, by = "cell")

#define z-axis
long$z <- as.numeric(factor(long$go_id))

save(long, axis_label, file = "3D.RData")
#clean environment
load("3D.RData")

long$tpm <- signif(long$tpm, 4)
long$x <- signif(long$x, 4)
long$y <- signif(long$y, 4)

#trim
long <- long %>% filter(tpm > 1)

#prune
set.seed(1)
long <- long %>% sample_frac(0.6)

p <- plot_ly(
  data = long,
  x = ~x,
  y = ~y,
  z = ~z,
  type = "scatter3d",
  mode = "markers",
  color = ~celltype,
  colors = "Spectral",
  size = ~tpm,
  sizes = c(1, 10),
  marker = list(opacity = 0.5),
  showlegend = TRUE
)

p <- p %>% layout(
  title = "3D Cell Towers (GO composition)",
  scene = list(
    xaxis = list(title = paste0(axis_label, "1")),
    yaxis = list(title = paste0(axis_label, "2")),
    zaxis = list(title = "GO term index")
  )
)

p <- p %>% partial_bundle()
p
