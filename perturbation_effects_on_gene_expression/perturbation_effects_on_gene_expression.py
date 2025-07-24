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
warnings.filterwarnings('ignore')

# Configuration
DATA_PATH = "adata_Training.h5ad"
OUTPUT_DIR = "manifold_visualizations"
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

# ==================== PERTURBATION EFFECT ANALYSIS ====================
@timeit
def analyze_perturbation_effects(adata, expression_matrix, gene_names):
    """Analyze perturbation effects with magnitude calculations"""
    print("Analyzing perturbation effects...")
    
    # Calculate control mean
    control_mask = adata.obs['target_gene'] == 'non-targeting'
    control_mean = expression_matrix[control_mask].mean(axis=0)
    control_std = expression_matrix[control_mask].std(axis=0)
    
    # Calculate delta expression
    delta_expression = expression_matrix - control_mean
    
    # Calculate perturbation effect magnitudes
    perturbation_effects = {}
    perturbation_stats = {}
    
    for pert in adata.obs['target_gene'].unique():
        if pert == 'non-targeting':
            continue
            
        pert_mask = adata.obs['target_gene'] == pert
        pert_delta = delta_expression[pert_mask]
        
        # Calculate statistics
        effect_magnitude = np.sqrt(np.sum(pert_delta**2, axis=1))  # L2 norm
        mean_effect = np.mean(effect_magnitude)
        std_effect = np.std(effect_magnitude)
        
        # Top differentially expressed genes
        mean_delta_per_gene = np.mean(pert_delta, axis=0)
        top_up_genes_idx = np.argsort(mean_delta_per_gene)[-10:][::-1]
        top_down_genes_idx = np.argsort(mean_delta_per_gene)[:10]
        
        perturbation_effects[pert] = {
            'magnitude_per_cell': effect_magnitude,
            'mean_magnitude': mean_effect,
            'std_magnitude': std_effect,
            'delta_expression': pert_delta,
            'mean_delta_per_gene': mean_delta_per_gene,
            'top_up_genes': [(gene_names[i], mean_delta_per_gene[i]) for i in top_up_genes_idx],
            'top_down_genes': [(gene_names[i], mean_delta_per_gene[i]) for i in top_down_genes_idx]
        }
        
        perturbation_stats[pert] = {
            'n_cells': np.sum(pert_mask),
            'mean_magnitude': mean_effect,
            'std_magnitude': std_effect
        }
    
    return delta_expression, perturbation_effects, perturbation_stats, control_mean

# ==================== DELTA GENE EXPRESSION MANIFOLD (ENHANCED) ====================
@timeit
def create_delta_gene_expression_manifold(adata, expression_matrix, subsample=None):
    print("Creating delta gene expression manifold...")
    
    # Subsample if requested
    if subsample and subsample < expression_matrix.shape[0]:
        rng = np.random.default_rng(42)
        sample_idx = rng.choice(expression_matrix.shape[0], size=subsample, replace=False)
        expression_matrix = expression_matrix[sample_idx]
        adata = adata[sample_idx].copy()
    
    # Calculate control mean
    control_mask = adata.obs['target_gene'] == 'non-targeting'
    control_mean = expression_matrix[control_mask].mean(axis=0)
    
    # Calculate delta expression
    delta_expression = expression_matrix - control_mean
    
    # Reduce dimensionality using UMAP
    print("Running UMAP...")
    reducer = UMAP(
        n_components=3,
        n_neighbors=15,
        min_dist=0.1,
        metric='cosine',
        random_state=42,
        verbose=True
    )
    delta_3d = reducer.fit_transform(delta_expression)
    
    return delta_3d, delta_expression, adata

# ==================== VISUAL ENCODING AND PLOT PREPARATION ====================
def prepare_enhanced_plot_data(delta_3d, adata, delta_expression, gene_names, perturbation_effects):
    print("Preparing enhanced plot data with visual encodings...")
    
    perturbations = adata.obs['target_gene'].values
    cell_ids = adata.obs.index.values
    
    # Calculate perturbation effect magnitudes for each cell
    effect_magnitudes = []
    for i, pert in enumerate(perturbations):
        if pert == 'non-targeting':
            effect_magnitudes.append(0)
        else:
            # Calculate L2 norm of delta expression for this cell
            magnitude = np.sqrt(np.sum(delta_expression[i]**2))
            effect_magnitudes.append(magnitude)
    
    effect_magnitudes = np.array(effect_magnitudes)
    
    # Prepare hover text with detailed information
    top_n = 5
    abs_delta = np.abs(delta_expression)
    top_genes_idx = np.argpartition(abs_delta, -top_n, axis=1)[:, -top_n:]
    
    hover_texts = []
    for i in range(len(delta_3d)):
        top_idx = top_genes_idx[i]
        top_genes = gene_names[top_idx]
        top_values = delta_expression[i, top_idx]
        
        hover_text = (
            f"<b>Cell ID:</b> {cell_ids[i]}<br>"
            f"<b>Perturbation:</b> {perturbations[i]}<br>"
            f"<b>Effect Magnitude:</b> {effect_magnitudes[i]:.2f}<br>"
            f"<b>Top Changed Genes:</b><br>"
        )
        hover_text += "<br>".join(
            f"&nbsp;&nbsp;{gene}: {val:.2f}" 
            for gene, val in zip(top_genes, top_values)
        )
        hover_texts.append(hover_text)
    
    # Create enhanced DataFrame
    plot_df = pd.DataFrame({
        'x': delta_3d[:, 0],
        'y': delta_3d[:, 1],
        'z': delta_3d[:, 2],
        'perturbation': perturbations,
        'effect_magnitude': effect_magnitudes,
        'hover_text': hover_texts,
        'is_control': perturbations == 'non-targeting',
        'cell_id': cell_ids
    })
    
    return plot_df

# ==================== DIFFERENTIAL GENE EXPRESSION HEATMAP ====================
def create_differential_expression_heatmap(perturbation_effects, gene_names, top_n=20):
    """Create heatmap of top differentially expressed genes across perturbations"""
    print("Creating differential expression heatmap...")
    
    # Collect top genes across all perturbations
    all_top_genes = set()
    for pert_data in perturbation_effects.values():
        top_up = [gene for gene, _ in pert_data['top_up_genes'][:top_n//2]]
        top_down = [gene for gene, _ in pert_data['top_down_genes'][:top_n//2]]
        all_top_genes.update(top_up + top_down)
    
    all_top_genes = list(all_top_genes)
    
    # Create matrix for heatmap
    heatmap_data = []
    perturbation_names = []
    
    for pert, pert_data in perturbation_effects.items():
        mean_delta = pert_data['mean_delta_per_gene']
        gene_indices = [np.where(gene_names == gene)[0][0] for gene in all_top_genes if gene in gene_names]
        row_data = [mean_delta[idx] for idx in gene_indices]
        heatmap_data.append(row_data)
        perturbation_names.append(pert)
    
    # Create Plotly heatmap
    fig_heatmap = go.Figure(data=go.Heatmap(
        z=heatmap_data,
        x=[gene for gene in all_top_genes if gene in gene_names][:len(heatmap_data[0])],
        y=perturbation_names,
        colorscale='RdBu_r',
        zmid=0,
        colorbar=dict(title="Mean Delta Expression")
    ))
    
    fig_heatmap.update_layout(
        title='Top Differentially Expressed Genes Across Perturbations',
        xaxis_title='Genes',
        yaxis_title='Perturbations',
        height=600,
        width=1000
    )
    
    return fig_heatmap

# ==================== SUMMARY PLOTS ====================
def create_summary_plots(perturbation_effects, perturbation_stats):
    """Create summary plots for perturbation effects"""
    print("Creating summary plots...")
    
    # Prepare data for summary plots
    pert_names = list(perturbation_stats.keys())
    mean_magnitudes = [perturbation_stats[p]['mean_magnitude'] for p in pert_names]
    std_magnitudes = [perturbation_stats[p]['std_magnitude'] for p in pert_names]
    n_cells = [perturbation_stats[p]['n_cells'] for p in pert_names]
    
    # Create subplots
    fig_summary = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Effect Magnitude by Perturbation', 'Cell Count by Perturbation',
                       'Effect Magnitude Distribution', 'Effect Magnitude vs Cell Count'),
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}]]
    )
    
    # Bar plot of effect magnitudes
    fig_summary.add_trace(
        go.Bar(x=pert_names, y=mean_magnitudes, 
               error_y=dict(type='data', array=std_magnitudes),
               name='Effect Magnitude',
               marker_color='lightblue'),
        row=1, col=1
    )
    
    # Bar plot of cell counts
    fig_summary.add_trace(
        go.Bar(x=pert_names, y=n_cells, 
               name='Cell Count',
               marker_color='lightcoral'),
        row=1, col=2
    )
    
    # Histogram of effect magnitudes
    all_magnitudes = []
    all_pert_labels = []
    for pert, pert_data in perturbation_effects.items():
        all_magnitudes.extend(pert_data['magnitude_per_cell'])
        all_pert_labels.extend([pert] * len(pert_data['magnitude_per_cell']))
    
    fig_summary.add_trace(
        go.Histogram(x=all_magnitudes, nbinsx=50, 
                    name='Magnitude Distribution',
                    marker_color='lightgreen',
                    opacity=0.7),
        row=2, col=1
    )
    
    # Scatter plot of magnitude vs cell count
    fig_summary.add_trace(
        go.Scatter(x=n_cells, y=mean_magnitudes,
                  mode='markers+text',
                  text=pert_names,
                  textposition='top center',
                  name='Magnitude vs Count',
                  marker=dict(size=10, color='purple')),
        row=2, col=2
    )
    
    fig_summary.update_layout(height=800, width=1200, 
                             title_text="Perturbation Effect Summary",
                             showlegend=False)
    
    return fig_summary

# ==================== FACETED VIEWS ====================
def create_faceted_views(plot_df, perturbation_effects, top_perturbations=8):
    """Create faceted views by perturbation"""
    print("Creating faceted views...")
    
    # Select top perturbations by effect magnitude
    pert_magnitudes = [(pert, data['mean_magnitude']) 
                      for pert, data in perturbation_effects.items()]
    pert_magnitudes.sort(key=lambda x: x[1], reverse=True)
    top_perts = [p[0] for p in pert_magnitudes[:top_perturbations]]
    
    # Create subplot grid
    n_cols = 4
    n_rows = (len(top_perts) + n_cols - 1) // n_cols
    
    fig_faceted = make_subplots(
        rows=n_rows, cols=n_cols,
        subplot_titles=top_perts,
        specs=[[{"type": "scatter3d"} for _ in range(n_cols)] for _ in range(n_rows)],
        vertical_spacing=0.1,
        horizontal_spacing=0.05
    )
    
    for i, pert in enumerate(top_perts):
        row = i // n_cols + 1
        col = i % n_cols + 1
        
        # Get data for this perturbation
        pert_data = plot_df[plot_df['perturbation'] == pert]
        control_data = plot_df[plot_df['is_control']]
        
        # Add control points (gray)
        if not control_data.empty:
            fig_faceted.add_trace(
                go.Scatter3d(
                    x=control_data['x'], y=control_data['y'], z=control_data['z'],
                    mode='markers',
                    marker=dict(size=2, color='gray', opacity=0.3),
                    name=f'{pert}_control',
                    showlegend=False
                ),
                row=row, col=col
            )
        
        # Add perturbation points (colored by magnitude)
        if not pert_data.empty:
            fig_faceted.add_trace(
                go.Scatter3d(
                    x=pert_data['x'], y=pert_data['y'], z=pert_data['z'],
                    mode='markers',
                    marker=dict(
                        size=4,
                        color=pert_data['effect_magnitude'],
                        colorscale='Viridis',
                        opacity=0.8,
                        showscale=i == 0  # Show colorbar only for first subplot
                    ),
                    text=pert_data['hover_text'],
                    hoverinfo='text',
                    name=pert,
                    showlegend=False
                ),
                row=row, col=col
            )
    
    fig_faceted.update_layout(
        height=300 * n_rows,
        width=1400,
        title_text="Faceted View: Top Perturbations by Effect Magnitude"
    )
    
    return fig_faceted

# ==================== MAIN ENHANCED VISUALIZATION ====================
def create_enhanced_interactive_manifold(plot_df, perturbation_effects):
    """Create enhanced interactive 3D manifold with visual encoding"""
    print("Creating enhanced interactive 3D manifold visualization...")
    
    fig = go.Figure()
    
    # Add control cells
    control_df = plot_df[plot_df['is_control']]
    if not control_df.empty:
        fig.add_trace(go.Scatter3d(
            x=control_df['x'], y=control_df['y'], z=control_df['z'],
            mode='markers',
            marker=dict(
                size=4,
                color='rgba(150, 150, 150, 0.5)',
                line=dict(width=0)
            ),
            text=control_df['hover_text'],
            hoverinfo='text',
            name='Control (non-targeting)',
            hoverlabel=dict(bgcolor='white', font_size=12, font_family="Arial")
        ))
    
    # Add perturbation cells with size and color encoding
    for pert, pert_data in perturbation_effects.items():
        pert_df = plot_df[plot_df['perturbation'] == pert]
        if pert_df.empty:
            continue
            
        # Use effect magnitude for both size and color intensity
        magnitudes = pert_df['effect_magnitude'].values
        normalized_magnitudes = (magnitudes - magnitudes.min()) / (magnitudes.max() - magnitudes.min() + 1e-8)
        
        # Size encoding: 3-8 pixels based on magnitude
        sizes = 3 + 5 * normalized_magnitudes
        
        fig.add_trace(go.Scatter3d(
            x=pert_df['x'], y=pert_df['y'], z=pert_df['z'],
            mode='markers',
            marker=dict(
                size=sizes,
                color=magnitudes,
                colorscale='Viridis',
                opacity=0.8,
                line=dict(width=0.5, color='white'),
                showscale=False,
                colorbar=dict(title="Effect Magnitude") if pert == list(perturbation_effects.keys())[0] else None
            ),
            text=pert_df['hover_text'],
            hoverinfo='text',
            name=f'{pert} (n={len(pert_df)})',
            hoverlabel=dict(bgcolor='white', font_size=12, font_family="Arial")
        ))
    
    # Create dropdown for filtering
    dropdown_buttons = []
    
    # All perturbations visible
    all_visible = [True] * len(fig.data)
    dropdown_buttons.append(dict(
        label="All Perturbations",
        method="update",
        args=[{"visible": all_visible}]
    ))
    
    # Individual perturbation filters
    for i, pert in enumerate(['Control (non-targeting)'] + list(perturbation_effects.keys())):
        visible = [False] * len(fig.data)
        visible[i] = True  # Show only this perturbation
        visible[0] = True  # Always show control for reference
        
        dropdown_buttons.append(dict(
            label=f"Show {pert}",
            method="update",
            args=[{"visible": visible}]
        ))
    
    # Update layout with enhanced features
    fig.update_layout(
        title=dict(
            text='<b>Enhanced 3D Manifold: Delta Gene Expression with Visual Encoding</b><br>'
                 '<i>Point size and color represent perturbation effect magnitude</i>',
            x=0.5, y=0.95, xanchor='center', yanchor='top', font=dict(size=16)
        ),
        scene=dict(
            xaxis_title='UMAP 1 (Delta Expression Axis)',
            yaxis_title='UMAP 2 (Delta Expression Axis)',
            zaxis_title='UMAP 3 (Delta Expression Axis)',
            camera=dict(eye=dict(x=1.5, y=1.5, z=0.8)),
            xaxis=dict(backgroundcolor='rgba(240,240,240,0.5)', gridcolor='white', showbackground=True),
            yaxis=dict(backgroundcolor='rgba(240,240,240,0.5)', gridcolor='white', showbackground=True),
            zaxis=dict(backgroundcolor='rgba(240,240,240,0.5)', gridcolor='white', showbackground=True)
        ),
        margin=dict(l=0, r=0, b=0, t=120),
        height=900,
        width=1200,
        legend=dict(
            title='<b>Perturbations</b>',
            itemsizing='constant',
            yanchor="top", y=0.99,
            xanchor="left", x=0.01
        ),
        updatemenus=[
            dict(
                buttons=dropdown_buttons,
                direction="down",
                showactive=True,
                x=0.02, xanchor="left",
                y=0.85, yanchor="top",
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="rgba(0,0,0,0.2)",
                font=dict(size=12)
            )
        ]
    )
    
    # Add annotation
    fig.add_annotation(
        text=f"Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}<br>"
             f"Visual encoding: Size ∝ Effect magnitude, Color ∝ Effect magnitude",
        xref="paper", yref="paper",
        x=0.02, y=0.02,
        showarrow=False,
        font=dict(size=10)
    )
    
    return fig

# ==================== EXECUTION ====================
print("Starting enhanced analysis...")

# Analyze perturbation effects
delta_expression, perturbation_effects, perturbation_stats, control_mean = analyze_perturbation_effects(
    adata, expression_matrix, gene_names
)

# Create manifold
delta_3d, delta_expression_sub, subsampled_adata = create_delta_gene_expression_manifold(
    adata, expression_matrix, subsample=10000
)

# Prepare enhanced plot data
plot_df = prepare_enhanced_plot_data(delta_3d, subsampled_adata, delta_expression_sub, 
                                   gene_names, perturbation_effects)

print("Creating all visualizations...")

# Create main enhanced visualization
fig_main = create_enhanced_interactive_manifold(plot_df, perturbation_effects)
output_file_main = os.path.join(OUTPUT_DIR, "enhanced_delta_manifold.html")
fig_main.write_html(output_file_main, include_plotlyjs='cdn')
print(f"Enhanced manifold saved to {output_file_main}")

# Create differential expression heatmap
fig_heatmap = create_differential_expression_heatmap(perturbation_effects, gene_names)
output_file_heatmap = os.path.join(OUTPUT_DIR, "differential_expression_heatmap.html")
fig_heatmap.write_html(output_file_heatmap, include_plotlyjs='cdn')
print(f"Differential expression heatmap saved to {output_file_heatmap}")

# Create summary plots
fig_summary = create_summary_plots(perturbation_effects, perturbation_stats)
output_file_summary = os.path.join(OUTPUT_DIR, "perturbation_summary.html")
fig_summary.write_html(output_file_summary, include_plotlyjs='cdn')
print(f"Summary plots saved to {output_file_summary}")

# Create faceted views
fig_faceted = create_faceted_views(plot_df, perturbation_effects)
output_file_faceted = os.path.join(OUTPUT_DIR, "faceted_perturbation_views.html")
fig_faceted.write_html(output_file_faceted, include_plotlyjs='cdn')
print(f"Faceted views saved to {output_file_faceted}")

# Create summary statistics report
print("\n=== PERTURBATION ANALYSIS SUMMARY ===")
print(f"Total perturbations analyzed: {len(perturbation_effects)}")
print(f"Control cells: {np.sum(adata.obs['target_gene'] == 'non-targeting')}")
print(f"Total genes: {len(gene_names)}")

print("\nTop 5 perturbations by effect magnitude:")
sorted_perts = sorted(perturbation_stats.items(), key=lambda x: x[1]['mean_magnitude'], reverse=True)
for i, (pert, stats) in enumerate(sorted_perts[:5]):
    print(f"{i+1}. {pert}: {stats['mean_magnitude']:.3f} ± {stats['std_magnitude']:.3f} "
          f"(n={stats['n_cells']} cells)")

print(f"\nAll visualizations saved to {OUTPUT_DIR}/")
print("Analysis complete!")