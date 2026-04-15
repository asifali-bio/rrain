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

#log-normalize
go_matrix_norm <- log1p(go_matrix)

#scale per cell
cs <- colSums(go_matrix_norm)
cs[cs == 0] <- 1
go_matrix_norm <- t(t(go_matrix_norm) / cs)

#remove zero-variance domains
go_var <- apply(go_matrix_norm, 1, var)
go_matrix_filt <- go_matrix_norm[go_var > 0, ]

#remove zero-variance cells
cell_var <- apply(go_matrix_filt, 2, var)
go_matrix_filt <- go_matrix_filt[, cell_var > 0]

#PCA
pca <- prcomp(t(go_matrix_filt), scale. = TRUE)

#preview PCA
plot(pca$x[,1], pca$x[,2],
     col = "blue",
     pch = 16,
     main = "PCA of GO Matrix")

#define valid cells
valid_cells <- colnames(go_matrix_filt)

#subset Seurat object
pbmc0 <- pbmc
pbmc0 <- subset(pbmc0, cells = valid_cells)
#force alignment
pbmc0 <- pbmc0[, valid_cells]

#overlay Seurat clustering
pbmc0 <- NormalizeData(pbmc0)
pbmc0 <- FindVariableFeatures(pbmc0)
pbmc0 <- ScaleData(pbmc0)
pbmc0 <- RunPCA(pbmc0)
pbmc0 <- FindNeighbors(pbmc0)
pbmc0 <- FindClusters(pbmc0)

clusters <- pbmc0$seurat_clusters

plot(pca$x[,1], pca$x[,2],
     col = clusters,
     pch = 16,
     main = "GO PCA colored by Seurat clusters")

#flag
USE_AZIMUTH <- FALSE

if (USE_AZIMUTH) {
  library(SeuratData)
  
  #load reference
  reference <- LoadData("pbmcref", type = "azimuth")
  
  #run mapping
  anchors <- FindTransferAnchors(
    reference = reference,
    query = pbmc0,
    dims = 1:30
  )
  
  pbmc0 <- MapQuery(
    anchorset = anchors,
    query = pbmc0,
    reference = reference,
    refdata = list(celltype = "celltype")
  )
  
  #build metadata
  meta_df <- data.frame(
    cell = colnames(pbmc0),
    cluster = as.factor(pbmc0$seurat_clusters),
    celltype = pbmc0$predicted.celltype
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
coords <- as.data.frame(pca$x[,1:2])
#enforce order and identity
coords <- coords[valid_cells, ]
coords$cell <- rownames(coords)
colnames(coords)[1:2] <- c("x", "y")

#safety
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

save(long, file = "3D.RData")
#clean environment
load("3D.RData")

long$tpm <- signif(long$tpm, 4)
long$x <- signif(long$x, 4)
long$y <- signif(long$y, 4)

#trim
long <- long %>% filter(tpm > 1)

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
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    zaxis = list(title = "GO term index")
  )
)

p <- p %>% partial_bundle()
p
