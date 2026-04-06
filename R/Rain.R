library(tidyr)
library(dplyr)
library(plotly)


load("pfam.RData")

new3 = new2
rownames(new3) = new3[,1]
new3 = new3[,-c(1,2)]
#set NA to 0
new3[is.na(new3)] = 0

new4 = new2
rownames(new4) = new4[,1]
new4 = new4[,-c(1,2)]

X = new3
X = t(X)

#remove zero-variance Pfam
X = X[, apply(X, 2, var) != 0]



#UMAP
library(uwot)

umap_coords = umap(X, n_neighbors = 3, min_dist = 0.2)
coords = as.data.frame(umap_coords)
colnames(coords) = c("UMAP1", "UMAP2")
coords$species = rownames(X)

#preview UMAP
ggplot(coords, aes(x = UMAP1, y = UMAP2, label = species)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.5) +
  theme_classic() +
  labs(title = "UMAP of species by Pfam composition")



#PCA
pca = prcomp(X, scale. = TRUE)
coords = as.data.frame(pca$x[,1:2])
coords$species = rownames(coords)

#preview PCA
ggplot(coords, aes(x = PC1, y = PC2, label = species)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.5) +
  theme_classic() +
  labs(title = "PCA of species by Pfam composition")

#3D
colnames(coords)[1:2] = c("x", "y")

df = as.data.frame(new4)
df$pfam = rownames(df)

long = pivot_longer(df, -pfam,
                    names_to = "species",
                    values_to = "tpm")

#remove NA so towers have gaps
long = long %>% filter(!is.na(tpm))

#attach PCA coordinates
long = merge(long, coords, by = "species")

long$z = as.numeric(factor(long$pfam, levels = rownames(new4)))

p = plot_ly(
  data = long,
  x = ~x,
  y = ~y,
  z = ~z,
  type = "scatter3d",
  mode = "markers",
  color = ~species,
  colors = "Spectral",
  size = ~tpm,
  sizes = c(1, 12),
  marker = list(opacity = 0.8)
)

p = p %>% layout(
  title = "3D Species Towers (Pfam composition)",
  scene = list(
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    zaxis = list(title = "Pfam index")
  )
)

p
