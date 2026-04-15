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

pfam <- read.delim("pfam_output.tsv", header = FALSE)

pfam_subset <- pfam[, c(1, 5, 6)]
colnames(pfam_subset) <- c("protein", "db", "domain")

pfam_subset$uniprot <- sub(".*\\|(.*?)\\|.*", "\\1", pfam_subset$protein)


mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  mart = mart
)

#clean up
mapping <- mapping[mapping$uniprotswissprot != "", ]
mapping <- mapping[mapping$hgnc_symbol != "", ]

#join protein → gene → domain
pfam_gene <- merge(
  pfam_subset,
  mapping,
  by.x = "uniprot",
  by.y = "uniprotswissprot"
)

#count data
expr <- GetAssayData(pbmc, layer = "counts")

#keep only genes in PBMC
pbmc_genes <- rownames(expr)
pfam_gene <- pfam_gene[pfam_gene$hgnc_symbol %in% pbmc_genes, ]

#build gene → domain table
gene_domain <- pfam_gene[, c("hgnc_symbol", "domain")]

#build domain matrix
domains <- unique(gene_domain$domain)

domain_matrix <- matrix(
  0,
  nrow = length(domains),
  ncol = ncol(expr)
)

rownames(domain_matrix) <- domains
colnames(domain_matrix) <- colnames(expr)

#aggregate expression
for (d in domains) {
  genes_in_d <- gene_domain$hgnc_symbol[gene_domain$domain == d]
  
  genes_in_d <- intersect(genes_in_d, rownames(expr))
  
  if (length(genes_in_d) > 0) {
    domain_matrix[d, ] <- colSums(expr[genes_in_d, , drop = FALSE])
  }
}

#domains × cells
dim(domain_matrix)

#log-normalize
domain_matrix_norm <- log1p(domain_matrix)

#scale per cell
cs <- colSums(domain_matrix_norm)
cs[cs == 0] <- 1
domain_matrix_norm <- t(t(domain_matrix_norm) / cs)

#remove zero-variance domains
domain_var <- apply(domain_matrix_norm, 1, var)
domain_matrix_filt <- domain_matrix_norm[domain_var > 0, ]

#remove zero-variance cells
cell_var <- apply(domain_matrix_filt, 2, var)
domain_matrix_filt <- domain_matrix_filt[, cell_var > 0]

#PCA
pca <- prcomp(t(domain_matrix_filt), scale. = TRUE)

#preview PCA
plot(pca$x[,1], pca$x[,2],
     col = "blue",
     pch = 16,
     main = "PCA of Domain Matrix")

#define valid cells
valid_cells <- colnames(domain_matrix_filt)

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
     main = "Domain PCA colored by Seurat clusters")

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
df <- as.data.frame(domain_matrix)
df$domain <- rownames(df)

long <- pivot_longer(
  df,
  -domain,
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
long$z <- as.numeric(factor(long$domain))

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
  title = "3D Cell Towers (Domain composition)",
  scene = list(
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    zaxis = list(title = "Domain index")
  )
)

p <- p %>% partial_bundle()
p
