# Virtual Cell Challenge 2025 Codebase and Pipelines

This repository serves as a comprehensive collection of all codebases and analysis pipelines developed as part of the Virtual Cell Challenge 2025.

---
***

# Manifold Visualization Pipeline

This pipeline visualizes the underlying structure (manifold) of gene expression patterns and perturbation effects.

## Overview Card

| **Component** | **Description** |
|---------------|-----------------|
| **Purpose** | Visualize manifold structure of gene expression patterns and perturbation effects |
| **Input Format** | Gene expression matrix (cells × genes) in h5ad format |
| **Input File** | `adata_Training.h5ad` (AnnData format for single-cell data) |
| **Output** | Interactive 3D visualization of gene and perturbation manifolds |

## Pipeline Components

| **Step** | **Method** | **Purpose** | **Output** |
|----------|------------|-------------|------------|
| **1** | Principal Component Analysis (PCA) | Reduce high-dimensional gene expression data | Gene Embeddings (compact representation) |
| **2** | UMAP (Uniform Manifold Approximation and Projection) | Further reduce embeddings to 3D for visualization | 3D coordinates for each gene |
| **3** | Perturbation Embedding (PCA → UMAP → 3D visualization) | Create representative points for each perturbation | Perturbation manifold coordinates |

## Visualization Features

### Gene Manifold
- **Points**: Each point represents a gene
- **Proximity**: Genes with similar expression patterns across cells are closer together
- **Distance**: Reflects similarity in expression patterns
- **Interactivity**: Hover information showing gene names
- **Interpretation**: Far apart points have very different expression patterns

### Perturbation Manifold
- **Creation Process**:
  - Find corresponding gene's PCA embedding
  - If perturbation isn't in gene list → use average of all gene embeddings
  - Average across all samples with that perturbation
- **Result**: Representative point for each perturbation's effect
- **Proximity**: Perturbations with similar transcriptomic effects are closer
- **Analysis**: Reveals functional relationships and pathway groupings

## Biological Insights

| **Insight Type** | **Description** |
|------------------|-----------------|
| **Functional Relationships** | Structure reveals functional relationships between genes |
| **Pathway Identification** | Groups of genes affecting similar pathways can be identified based on manifold position |
| **Perturbation Effects** | Similar perturbation effects cluster together on the manifold |

## Resources

### Figure
![Manifold Visualization](https://github.com/nbahador/single_cell_genetic_perturbations/blob/main/Manifold_Visualization_Pipeline/Manifold_Visualization_Fig.jpg)

### Interactive Visualization
[Interactive Manifold Visualizations](https://github.com/nbahador/single_cell_genetic_perturbations/blob/main/Manifold_Visualization_Pipeline/manifold_visualizations.html)

## Technical Details

| **Aspect** | **Details** |
|------------|-------------|
| **Data Format** | AnnData (h5ad) - standard for single-cell gene expression |
| **Dimensionality Reduction** | PCA → UMAP → 3D visualization |
| **Embedding Strategy** | Gene-based embeddings with perturbation averaging |
| **Visualization Type** | Interactive 3D scatter plot with hover information |

---
***
