<p align="center">
  <img src="assets/rain_cloud_moji.svg" width="220">
</p>

# Reverse Rain (RRain)

**RRain** is a visualization framework for representing **single-cell RNA-seq data in functional space** using protein domains.

It extends the PLANT framework by projecting cells into a low-dimensional embedding and stacking their **domain-level features vertically**, creating “towers” that resemble **upside-down rain** emerging from a cloud of cells.

Instead of genes defining the primary structure, RRain aggregates expression into **protein domains (e.g., Pfam)**—allowing visualization of **functional composition per cell**.

---

## 🌐 Demo
https://asifali-bio.github.io/RRain/

---

## Concept

RRain separates the data into three intuitive dimensions:

- **(x, y)** → latent embedding of cells (domain space PCA/UMAP)
- **z** → feature index (protein domains)
- **size** → domain-level expression (aggregated counts / TPM)
- **color** → cell clusters (e.g., Seurat clustering)

Each cell becomes a **vertical tower of domains**, positioned according to similarity in functional space.

---

## Key idea

Traditional scRNA-seq pipelines operate in **gene space**, using highly variable genes for dimensionality reduction and clustering.

RRain instead constructs a parallel representation in **domain space**:

```
cells × genes → cells × domains
```

by mapping:

```
gene → protein → Pfam domain
```


This enables:

- Functional aggregation across genes
- Reduced redundancy (shared domains)
- Sharper boundaries between cell states
- Visualization of transitions via domain composition

---

## What it shows

RRain allows simultaneous visualization of:

- Global similarity between cells (embedding)
- Functional structure of each cell (tower composition)
- Domain presence/absence (gaps in towers)
- Domain magnitude (bubble size)
- Cluster identity (color via Seurat)

---

## Workflow

1. **Input data**

   - scRNA-seq count matrix (cells × genes)

2. **Domain annotation**

   - Run InterProScan on a reference proteome
   - Extract Pfam domains
   - Map UniProt → gene symbols (Ensembl / biomaRt)

3. **Aggregation**

   - Map genes → domains
   - Sum expression per domain per cell

```
cells × genes → cells × domains
```

4. **Embedding**

- PCA / UMAP on domain matrix

5. **Visualization**

- Convert to long format
- Plot 3D towers:

  - x, y = embedding
  - z = domain index
  - size = domain expression
  - color = cluster identity

---

## Relationship to PLANT

RRain builds directly on PLANT:

- **PLANT** → bulk / cross-sample domain composition
- **RRain** → single-cell domain composition

RRain extends the idea from:

> “one tower per sample”

to:

> “one tower per cell”

enabling exploration of **cellular heterogeneity and transitions**.

---

## Example use cases

- Single-cell RNA-seq (primary use case)
- Functional profiling of cell states
- Domain-level clustering and boundary detection
- Comparative analysis of cell populations

---

## 🔧 Prototype stage

- Current implementation uses point-based towers (Plotly)

Future work:
- continuous / volumetric tower rendering
- optimized large-scale visualization
- domain hierarchy-aware stacking

---

*Banner emoji designed by OpenMoji.*