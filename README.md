# Virtual Cell Challenge 2025 Codebase and Pipelines

This repository serves as a collection of all codebases and analysis pipelines developed as part of the Virtual Cell Challenge 2025.

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

# Differential Manifold Pipeline: Mapping of Control vs. Perturbed Single-Cell States

This analysis pipeline processes single-cell gene expression data to visualize and compare the effects of genetic perturbations.

## Cell Classification

| **Group** | **Definition** | **Purpose** |
|-----------|----------------|-------------|
| **Control cells** | Cells with 'non-targeting' perturbations | Baseline/reference group |
| **Perturbed cells** | Cells with specific genetic perturbations | Treatment group for comparison |

## Dimensionality Reduction Pipeline

To visualize high-dimensional gene expression data, we reduce it to 2D/3D using a two-step approach:

### A. Principal Component Analysis (PCA)

| **Parameter** | **Value** | **Purpose** |
|---------------|-----------|-------------|
| **Components** | Top 50 (n_pca=50) | Capture most variance in the data |
| **Scope** | All cells together | Ensure comparable coordinates |
| **Function** | Identifies main axes of variation | Provides initial dimensionality reduction |

### B. UMAP (Uniform Manifold Approximation and Projection)

| **Parameter** | **Value** | **Purpose** |
|---------------|-----------|-------------|
| **Input dimensions** | 50 PCA components | From previous PCA step |
| **Output dimensions** | 3D (n_umap=3) | Final visualization space |
| **n_neighbors** | 30 | Balances local vs global structure |
| **min_dist** | 0.1 | Controls how tightly points are packed |
| **Scope** | All cells together | Maintains consistency across groups |

## Statistical Comparison Analysis

### Density Estimation
- **Method**: Gaussian Kernel Density Estimation (KDE)
- **Purpose**: Estimate probability density of control and perturbed cells
- **Output**: Smooth maps showing concentration of each group

### Density Differences
- **Calculation**: Subtract control density from perturbed density at each point
- **Significance threshold**: Top 5% of absolute values marked as significant
- **Visualization**:
  - **Red regions**: Higher density in perturbed vs control
  - **Blue regions**: Higher density in control vs perturbed
  - **Black dots**: Statistically significant differences

## Interpretation

| **Analysis Component** | **What It Shows** |
|------------------------|-------------------|
| **3D UMAP space** | Where perturbed cells concentrate differently from controls |
| **2D heatmap** | Summary view of distribution differences |
| **Contour lines** | Statistically significant regions (p < 0.05) |

## Resources

### Figure
![2D Density Difference Map](https://github.com/nbahador/single_cell_genetic_perturbations/blob/main/Differential_Manifold_Mapping_of_Control_vs._Perturbed_Single-Cell_States/Differential_Manifold.jpg)

### Interactive Visualization
[Interactive Manifold Comparison Report](https://github.com/nbahador/single_cell_genetic_perturbations/blob/main/Differential_Manifold_Mapping_of_Control_vs._Perturbed_Single-Cell_States/manifold_comparison_report.html)

---
***

# Perturbation Effects on Gene Expression

> This pipeline visualizes and analyzes the effects of genetic perturbations on single-cell gene expression profiles.

| Feature | Description |
|---------|-------------|
| **UMAP Positioning** | Reveals similarity in gene expression changes - cells with similar perturbation responses cluster together |
| **Effect Magnitude** | Quantifies how strongly each perturbation modifies the overall expression landscape |
| **Gene-Level Changes** | Identifies which specific genes are most significantly affected by each perturbation |

## Applications

<div align="center">

| **Mechanistic Insights** | **Hit Prioritization** | **Gene Function Prediction** |
|:------------------------:|:----------------------:|:----------------------------:|
| Discover perturbations with similar molecular effects through clustering analysis | Rank perturbations by their overall impact magnitude for compound screening | Identify genes that consistently co-regulate across different perturbations |

</div>

## Analysis Workflow

### 1. Perturbation Effect Calculation
- **Differential Expression Analysis**: Computes expression deltas between perturbed and control (non-targeting) cells
- **Effect Magnitude Quantification**: Measures perturbation strength using L2 norm of delta expression vectors  
- **Gene Prioritization**: Identifies top up- and down-regulated genes for each perturbation condition

### 2. 3D Manifold Visualization
An interactive 3D UMAP projection revealing perturbation landscapes:

| Element | Representation |
|---------|----------------|
| **Points** | Individual cells in expression space |
| **Position** | UMAP coordinates based on expression change similarity |
| **Color** | Perturbation effect magnitude (viridis colormap) |
| **Size** | Effect strength (larger points = stronger perturbations) |
| **Interactivity** | Hover details: cell ID, perturbation type, effect size, top changed genes |

**Features**: Dropdown filtering by specific perturbations for focused analysis

### 3. Differential Expression Heatmap
Overview of gene expression changes across all perturbations:
- **Rows**: Individual perturbations
- **Columns**: Top differentially expressed genes
- **Color Scale**: Mean delta expression (red = upregulation, blue = downregulation)

### 4. Effect Magnitude Profiling
Statistical summary showing mean ± standard deviation of perturbation effect sizes across conditions.

## Visualizations
![Perturbation Effects on Gene Expression Changes](https://github.com/nbahador/single_cell_genetic_perturbations/blob/main/perturbation_effects_on_gene_expression/perturbation_effect_on_gene_expression_changes.jpg)

*Visualization showing the impact of genetic perturbations on cellular gene expression profiles*

### Interactive Visualizations

| Visualization Type | Preview | Link |
|-------------------|---------|------|
| **3D Manifold Explorer** | [**Launch Interactive 3D Manifold**](https://github.com/nbahador/single_cell_genetic_perturbations/blob/main/perturbation_effects_on_gene_expression/enhanced_delta_manifold.html) | Explore perturbation effects in 3D expression space with hover details and filtering |
| **Expression Heatmap** | [**Open Interactive Heatmap**](https://github.com/nbahador/single_cell_genetic_perturbations/blob/main/perturbation_effects_on_gene_expression/top_genes_heatmap.html) | Interactive heatmap of top differentially expressed genes across all perturbations |

---
***

# Transcriptomic Impact of Gene Perturbations Analysis Pipeline

## Core Goal

The primary aim of this analysis pipeline is to quantify the effect of individual gene perturbations on the transcriptome, focusing on identifying genes that become significantly upregulated or downregulated in response to each perturbation.

## Purpose of the Pipeline

This code is designed to summarize and visualize how different gene perturbations impact global gene expression. It helps answer key questions such as:

- How many genes are significantly altered by each perturbation?
- Are the changes primarily upregulation or downregulation of gene expression?
- What is the distribution of effect sizes—are the changes mostly modest or do some perturbations cause extreme shifts?
- How does a perturbation's strength (maximum effect size) relate to its breadth (number of significantly affected genes)?

Through these analyses, the pipeline characterizes the diversity of transcriptional responses triggered by different perturbations, reflecting their mechanisms of action and potency.

## Key Analytical Outputs

### Number of Significantly Affected Genes per Perturbation

Quantifies the total number of genes whose expression significantly changes due to each gene perturbation. Helps distinguish between strong, broad regulators and weak or targeted perturbations.

### Counts of Upregulated vs. Downregulated Genes

Shows the directionality of the effect: Does a perturbation mainly activate or suppress gene expression? This can reflect the biological role of the perturbed gene (e.g., activator vs. repressor).

### Effect Size Distributions

Examines the magnitude of transcriptional changes. Distributions are often bell-shaped or skewed, showing that:

- Most genes change modestly
- A few genes may experience extreme upregulation or downregulation, suggesting strong regulatory control

This helps identify master regulators or genes with broad vs. narrow influence.

### Strength vs. Coverage Plot

Plots the maximum effect size (strength) against the number of significantly affected genes (coverage). Bigger dots represent perturbations supported by more cell data, indicating higher confidence.

This plot reveals categories of perturbations:

- **Broad but modest** (many genes, small changes)
- **Focused but potent** (few genes, large changes)
- **Globally strong** (many genes, large changes)

## Interpretation and Insights

- Stronger or broader perturbations generally lead to more significant gene expression changes
- Perturbations can be grouped based on their transcriptional impact patterns
- Some perturbations show balanced effects, while others exhibit asymmetry, mainly activating or suppressing downstream genes
- By comparing gene expression in perturbed vs. control cells, the pipeline helps uncover gene regulatory roles and potential pathways involved
- The diversity in responses reflects the underlying biology and mechanisms of action of each gene perturbation

## Visualizations

### Static Figure
![Transcriptomic Impact Analysis](https://github.com/nbahador/single_cell_genetic_perturbations/blob/main/Transcriptomic_Impact_of_Gene_Perturbations/img.png)

### Interactive Visualization
[Interactive Perturbation Summary](https://github.com/nbahador/single_cell_genetic_perturbations/blob/main/Transcriptomic_Impact_of_Gene_Perturbations/perturbation_summary.html)

---
***
