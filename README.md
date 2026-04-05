<p align="center">
  <img src="assets/rain.svg" width="220">
</p>

# Reverse Rain (RRain)

**RRain** is a visualization framework for representing high-dimensional data as **vertical structures in latent space**.

It projects samples into a low-dimensional plane (e.g., PCA/UMAP), then stacks their features vertically—creating “towers” that resemble **upside-down rain** emerging from a cloud of points. Instead of points falling downward, features appear to **rise upward** from each sample’s position—like rain in reverse.

## 🌐 Demo
https://asifali-bio.github.io/RRain/

---

## Concept

RRain separates data into three intuitive dimensions:

* **(x, y)** → *latent embedding of samples*
* **z** → *feature index (e.g., Pfam domains, genes, variables)*
* **color / size** → *feature magnitude (e.g., TPM, counts)*

Each sample becomes a **vertical column (tower)** of its features, positioned according to similarity in latent space.

---

## What it shows

RRain allows you to simultaneously visualize:

* Global similarity between samples (via spatial proximity)
* Internal structure of each sample (via vertical patterns)
* Feature presence/absence (via sparsity or gaps)
* Feature magnitude (via size or color)

---

## Example use cases

* Protein domain composition across species
* Single-cell RNA-seq profiles
* Microbiome abundance data
* Any high-dimensional feature matrix

---

## Basic workflow

1. Start with a feature matrix:

   * rows = features
   * columns = samples

2. Preprocess:

   * handle missing values (e.g., NA → 0 for distance)
   * optionally filter low-variance features

3. Compute embedding:

   * PCA, UMAP, etc.

4. Convert to long format:

   * one row per (sample, feature)

5. Plot in 3D:

   * x, y = embedding
   * z = feature index
   * color/size = feature value

---

## 🔧 Prototype stage

* Eventually, continuous / volumetric tower rendering

---

*Banner emoji designed by OpenMoji.*