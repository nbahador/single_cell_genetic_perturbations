import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from umap import UMAP
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import seaborn as sns
import matplotlib.pyplot as plt
import os
from datetime import datetime
import time
import warnings
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import networkx as nx
warnings.filterwarnings('ignore')

# Configuration
DATA_PATH = "adata_Training.h5ad"
OUTPUT_DIR = os.path.normpath("interpretable_perturbation_analysis")  # Normalize path for OS
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Timer decorator
def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print(f'{method.__name__} executed in {te-ts:.2f} seconds')
        return result
    return timed

# ==================== DATA LOADING AND PREPROCESSING ====================
@timeit
def load_and_preprocess_data():
    print("Loading and preprocessing data...")
    adata = sc.read_h5ad(DATA_PATH)
    
    # Convert to dense matrix if sparse
    if hasattr(adata.X, 'toarray'):
        expression_matrix = adata.X.toarray()
    else:
        expression_matrix = adata.X
    
    return adata, expression_matrix, adata.var_names

adata, expression_matrix, gene_names = load_and_preprocess_data()

# ==================== BIOLOGICAL PERTURBATION ANALYSIS ====================
@timeit
def analyze_biological_perturbation_effects(adata, expression_matrix, gene_names):
    """Analyze perturbation effects with biological interpretation focus"""
    print("Analyzing biological perturbation effects...")
    
    # Calculate control statistics
    control_mask = adata.obs['target_gene'] == 'non-targeting'
    control_expression = expression_matrix[control_mask]
    control_mean = control_expression.mean(axis=0)
    control_std = control_expression.std(axis=0)
    
    # Avoid division by zero
    control_std = np.where(control_std == 0, 1, control_std)
    
    biological_effects = {}
    
    for pert in adata.obs['target_gene'].unique():
        if pert == 'non-targeting':
            continue
            
        print(f"  Analyzing {pert}...")
        pert_mask = adata.obs['target_gene'] == pert
        pert_expression = expression_matrix[pert_mask]
        
        if pert_expression.shape[0] < 3:  # Skip if too few cells
            continue
        
        # Calculate fold changes and z-scores
        pert_mean = pert_expression.mean(axis=0)
        fold_change = np.log2((pert_mean + 1) / (control_mean + 1))  # Log2 fold change
        z_score = (pert_mean - control_mean) / control_std
        
        # Statistical significance testing
        p_values = []
        for i in range(len(gene_names)):
            if control_expression.shape[0] > 5 and pert_expression.shape[0] > 5:
                _, p_val = stats.ttest_ind(control_expression[:, i], pert_expression[:, i])
                p_values.append(p_val)
            else:
                p_values.append(1.0)
        
        p_values = np.array(p_values)
        
        # Identify significantly changed genes (|fold_change| > 0.5 and p < 0.05)
        significant_mask = (np.abs(fold_change) > 0.5) & (p_values < 0.05)
        
        # Categorize genes
        upregulated = (fold_change > 0.5) & (p_values < 0.05)
        downregulated = (fold_change < -0.5) & (p_values < 0.05)
        
        # Get top changed genes
        top_up_idx = np.argsort(fold_change)[-20:][::-1]
        top_down_idx = np.argsort(fold_change)[:20]
        
        biological_effects[pert] = {
            'fold_change': fold_change,
            'z_score': z_score,
            'p_values': p_values,
            'significant_genes': gene_names[significant_mask],
            'upregulated_genes': gene_names[upregulated],
            'downregulated_genes': gene_names[downregulated],
            'top_upregulated': [(gene_names[i], fold_change[i], p_values[i]) for i in top_up_idx if fold_change[i] > 0],
            'top_downregulated': [(gene_names[i], fold_change[i], p_values[i]) for i in top_down_idx if fold_change[i] < 0],
            'n_cells': pert_expression.shape[0],
            'n_significant': np.sum(significant_mask),
            'n_upregulated': np.sum(upregulated),
            'n_downregulated': np.sum(downregulated),
            'mean_expression': pert_mean,
            'expression_data': pert_expression
        }
    
    return biological_effects, control_mean, control_expression

# ==================== PERTURBATION-CENTRIC VISUALIZATION ====================
def create_perturbation_effect_summary(biological_effects):
    """Create comprehensive perturbation effect summary"""
    print("Creating perturbation effect summary...")
    
    perturbations = list(biological_effects.keys())
    
    # Prepare data
    n_significant = [biological_effects[p]['n_significant'] for p in perturbations]
    n_upregulated = [biological_effects[p]['n_upregulated'] for p in perturbations]
    n_downregulated = [biological_effects[p]['n_downregulated'] for p in perturbations]
    n_cells = [biological_effects[p]['n_cells'] for p in perturbations]
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[
            'Significantly Changed Genes per Perturbation',
            'Up vs Down-regulated Genes',
            'Effect Size Distribution',
            'Perturbation Strength vs Coverage'
        ],
        specs=[[{"type": "bar"}, {"type": "bar"}],
               [{"type": "histogram"}, {"type": "scatter"}]]
    )
    
    # 1. Bar chart of significant genes
    fig.add_trace(
        go.Bar(
            x=perturbations,
            y=n_significant,
            name='Significant Genes',
            marker_color='lightblue',
            hovertemplate='<b>%{x}</b><br>Significant genes: %{y}<extra></extra>'
        ),
        row=1, col=1
    )
    
    # 2. Stacked bar chart of up/down regulation
    fig.add_trace(
        go.Bar(
            x=perturbations,
            y=n_upregulated,
            name='Upregulated',
            marker_color='red',
            opacity=0.7
        ),
        row=1, col=2
    )
    
    fig.add_trace(
        go.Bar(
            x=perturbations,
            y=[-x for x in n_downregulated],  # Negative for visual separation
            name='Downregulated',
            marker_color='blue',
            opacity=0.7
        ),
        row=1, col=2
    )
    
    # 3. Distribution of effect sizes
    all_fold_changes = []
    all_pert_labels = []
    for pert, data in biological_effects.items():
        significant_fc = data['fold_change'][np.abs(data['fold_change']) > 0.5]
        all_fold_changes.extend(significant_fc)
        all_pert_labels.extend([pert] * len(significant_fc))
    
    fig.add_trace(
        go.Histogram(
            x=all_fold_changes,
            nbinsx=50,
            name='Effect Sizes',
            marker_color='green',
            opacity=0.7,
            hovertemplate='Log2 FC: %{x:.2f}<br>Count: %{y}<extra></extra>'
        ),
        row=2, col=1
    )
    
    # 4. Scatter plot: strength vs coverage
    max_fold_changes = [np.max(np.abs(biological_effects[p]['fold_change'])) for p in perturbations]
    
    fig.add_trace(
        go.Scatter(
            x=n_significant,
            y=max_fold_changes,
            mode='markers+text',
            text=perturbations,
            textposition='top center',
            marker=dict(
                size=[np.sqrt(x/10) for x in n_cells],  # Size by cell count
                color='purple',
                opacity=0.7
            ),
            name='Perturbations',
            hovertemplate='<b>%{text}</b><br>' +
                         'Significant genes: %{x}<br>' +
                         'Max effect size: %{y:.2f}<br>' +
                         '<extra></extra>'
        ),
        row=2, col=2
    )
    
    # Update layout
    fig.update_layout(
        height=800,
        width=1400,
        title_text="Perturbation Effects: Biological Impact Summary",
        showlegend=True
    )
    
    # Update axis labels
    fig.update_xaxes(title_text="Perturbations", row=1, col=1, tickangle=45)
    fig.update_xaxes(title_text="Perturbations", row=1, col=2, tickangle=45)
    fig.update_xaxes(title_text="Log2 Fold Change", row=2, col=1)
    fig.update_xaxes(title_text="Number of Significant Genes", row=2, col=2)
    
    fig.update_yaxes(title_text="Count", row=1, col=1)
    fig.update_yaxes(title_text="Gene Count", row=1, col=2)
    fig.update_yaxes(title_text="Frequency", row=2, col=1)
    fig.update_yaxes(title_text="Maximum Effect Size", row=2, col=2)
    
    return fig

# ==================== EXECUTION ====================
print("Starting interpretable biological analysis...")

# Run biological analysis
biological_effects, control_mean, control_expression = analyze_biological_perturbation_effects(
    adata, expression_matrix, gene_names
)

print(f"Analyzed {len(biological_effects)} perturbations")

# Create perturbation effect summary
fig_pert_summary = create_perturbation_effect_summary(biological_effects)

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Save the figure
output_path = os.path.join(OUTPUT_DIR, "perturbation_biological_summary.html")
fig_pert_summary.write_html(output_path, include_plotlyjs='cdn')
print(f"✓ Perturbation biological summary saved to {output_path}")

print("\n=== ANALYSIS COMPLETE ===")
print(f"Results saved to: {os.path.abspath(OUTPUT_DIR)}")