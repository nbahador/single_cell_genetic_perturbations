import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from umap import UMAP
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import gaussian_kde
import os
from datetime import datetime
import pickle

# Configuration
DATA_PATH = "adata_Training.h5ad"
OUTPUT_DIR = "enhanced_manifold_comparison"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==================== DATA LOADING ====================
print("Loading data...")
adata = sc.read_h5ad(DATA_PATH)
expression_matrix = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
control_mask = adata.obs['target_gene'] == 'non-targeting'
perturbed_mask = ~control_mask

# ==================== MANIFOLD CONSTRUCTION ====================
print("Building manifolds...")

def build_manifolds(expression_matrix, control_mask, perturbed_mask, n_pca=50, n_umap=3):
    # Combined PCA
    pca = PCA(n_components=n_pca)
    combined_pca = pca.fit_transform(expression_matrix)
    
    # Split into control and perturbed
    control_pca = combined_pca[control_mask]
    perturbed_pca = combined_pca[perturbed_mask]
    
    # Combined UMAP
    umap = UMAP(n_components=n_umap, random_state=42, n_neighbors=30, min_dist=0.1)
    combined_umap = umap.fit_transform(combined_pca)
    
    control_umap = combined_umap[control_mask]
    perturbed_umap = combined_umap[perturbed_mask]
    
    return control_umap, perturbed_umap, pca, umap

control_umap, perturbed_umap, pca_model, umap_model = build_manifolds(
    expression_matrix, control_mask, perturbed_mask
)

# ==================== DIFFERENCE ANALYSIS ====================
print("Calculating manifold differences...")

def calculate_differences(control_umap, perturbed_umap):
    # Density differences
    kde_control = gaussian_kde(control_umap[:,:2].T)
    kde_perturbed = gaussian_kde(perturbed_umap[:,:2].T)
    
    # Create grid for visualization
    grid_size = 100
    x_min, x_max = min(control_umap[:,0].min(), perturbed_umap[:,0].min()), max(control_umap[:,0].max(), perturbed_umap[:,0].max())
    y_min, y_max = min(control_umap[:,1].min(), perturbed_umap[:,1].min()), max(control_umap[:,1].max(), perturbed_umap[:,1].max())
    
    xx, yy = np.meshgrid(np.linspace(x_min, x_max, grid_size),
                        np.linspace(y_min, y_max, grid_size))
    
    # Evaluate densities
    grid_points = np.vstack([xx.ravel(), yy.ravel()])
    z_control = kde_control(grid_points)
    z_perturbed = kde_perturbed(grid_points)
    z_diff = z_perturbed - z_control
    
    # Find regions of significant difference
    threshold = np.percentile(np.abs(z_diff), 95)
    sig_diff = np.abs(z_diff) > threshold
    
    return {
        'grid': (xx, yy),
        'control_density': z_control,
        'perturbed_density': z_perturbed,
        'density_diff': z_diff,
        'sig_diff': sig_diff.reshape(xx.shape),  # Reshape to match grid
        'threshold': threshold
    }

diff_results = calculate_differences(control_umap, perturbed_umap)

# ==================== ENHANCED VISUALIZATIONS ====================
print("Creating enhanced visualizations...")

def create_3d_difference_plot(control_umap, perturbed_umap, diff_results):
    fig = go.Figure()
    
    # Control points (smaller and more transparent)
    fig.add_trace(go.Scatter3d(
        x=control_umap[:,0], y=control_umap[:,1], z=control_umap[:,2],
        mode='markers',
        marker=dict(size=2, color='blue', opacity=0.1),
        name='Control'
    ))
    
    # Perturbed points colored by density difference
    grid_pts = np.c_[diff_results['grid'][0].ravel(), diff_results['grid'][1].ravel()]
    tree = NearestNeighbors(n_neighbors=1).fit(grid_pts)
    _, indices = tree.kneighbors(perturbed_umap[:,:2])
    perturbed_colors = diff_results['density_diff'][indices].ravel()
    
    fig.add_trace(go.Scatter3d(
        x=perturbed_umap[:,0], y=perturbed_umap[:,1], z=perturbed_umap[:,2],
        mode='markers',
        marker=dict(
            size=4,
            color=perturbed_colors,
            colorscale='RdBu',
            cmin=-diff_results['threshold'],
            cmax=diff_results['threshold'],
            opacity=0.7,
            colorbar=dict(title='Density Difference')
        ),
        name='Perturbed'
    ))
    
    # Highlight significant difference regions
    sig_mask = diff_results['sig_diff'].ravel()  # Flatten for indexing
    sig_x = diff_results['grid'][0].ravel()[sig_mask]
    sig_y = diff_results['grid'][1].ravel()[sig_mask]
    
    fig.add_trace(go.Scatter3d(
        x=sig_x, y=sig_y, z=np.zeros_like(sig_x),
        mode='markers',
        marker=dict(size=3, color='black', opacity=0.3),
        name='Significant Difference'
    ))
    
    fig.update_layout(
        title='3D Manifold with Density Differences',
        scene=dict(
            xaxis_title='UMAP1',
            yaxis_title='UMAP2',
            zaxis_title='UMAP3',
            camera=dict(eye=dict(x=1.5, y=1.5, z=0.8))
        ),
        margin=dict(l=0, r=0, b=0, t=40)
    )
    
    return fig

def create_2d_difference_map(diff_results):
    fig = go.Figure()
    
    # Density difference contour
    fig.add_trace(go.Contour(
        x=diff_results['grid'][0][0],
        y=diff_results['grid'][1][:,0],
        z=diff_results['density_diff'].reshape(diff_results['grid'][0].shape),
        colorscale='RdBu',
        contours=dict(
            coloring='heatmap',
            showlabels=True
        ),
        colorbar=dict(title='Density Difference (Perturbed - Control)'),
        name='Density Difference'
    ))
    
    # Significant difference regions
    sig_diff = diff_results['sig_diff']
    fig.add_trace(go.Contour(
        x=diff_results['grid'][0][0],
        y=diff_results['grid'][1][:,0],
        z=sig_diff,
        contours=dict(
            coloring=None,
            showlines=False
        ),
        line=dict(width=0),
        showscale=False,
        opacity=0.3,
        name='Significant Regions'
    ))
    
    fig.update_layout(
        title='2D Density Difference Map',
        xaxis_title='UMAP1',
        yaxis_title='UMAP2',
        margin=dict(l=0, r=0, b=0, t=40)
    )
    
    return fig

def create_perturbation_vectors(control_umap, perturbed_umap, perturbed_obs):
    # Calculate control centroid
    control_centroid = np.median(control_umap, axis=0)
    
    # Group perturbations by target gene
    unique_genes = perturbed_obs['target_gene'].unique()
    gene_vectors = {}
    
    for gene in unique_genes:
        if gene != 'non-targeting':
            gene_mask = perturbed_obs['target_gene'] == gene
            if gene_mask.sum() > 0:
                gene_centroid = np.median(perturbed_umap[gene_mask], axis=0)
                gene_vectors[gene] = gene_centroid - control_centroid
    
    # Create vector plot
    fig = go.Figure()
    
    # Add control centroid
    fig.add_trace(go.Scatter3d(
        x=[control_centroid[0]],
        y=[control_centroid[1]],
        z=[control_centroid[2]],
        mode='markers',
        marker=dict(size=10, color='blue', symbol='cross'),
        name='Control Centroid'
    ))
    
    # Add perturbation vectors
    for gene, vector in gene_vectors.items():
        fig.add_trace(go.Scatter3d(
            x=[control_centroid[0], control_centroid[0] + vector[0]],
            y=[control_centroid[1], control_centroid[1] + vector[1]],
            z=[control_centroid[2], control_centroid[2] + vector[2]],
            mode='lines+markers',
            line=dict(width=3, color='red'),
            marker=dict(size=5),
            name=f'{gene} Perturbation'
        ))
    
    fig.update_layout(
        title='Perturbation Vectors from Control Centroid',
        scene=dict(
            xaxis_title='UMAP1',
            yaxis_title='UMAP2',
            zaxis_title='UMAP3'
        ),
        showlegend=True,
        margin=dict(l=0, r=0, b=0, t=40)
    )
    
    return fig

# Generate all visualizations
print("Generating 3D difference plot...")
diff_3d_fig = create_3d_difference_plot(control_umap, perturbed_umap, diff_results)

print("Generating 2D difference map...")
diff_2d_fig = create_2d_difference_map(diff_results)

print("Generating perturbation vectors...")
vector_fig = create_perturbation_vectors(control_umap, perturbed_umap, adata.obs[perturbed_mask])

# ==================== HTML REPORT ====================
print("Generating HTML report...")

html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Enhanced Manifold Comparison</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .plot-container {{ margin: 30px 0; border: 1px solid #ddd; padding: 15px; }}
        .plot-title {{ font-size: 18px; font-weight: bold; color: #2c3e50; }}
        .plot-description {{ margin: 10px 0; color: #555; }}
    </style>
</head>
<body>
    <h1>Enhanced Manifold Comparison Report</h1>
    <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    
    <div class="plot-container">
        <div class="plot-title">3D Manifold with Density Differences</div>
        <div class="plot-description">
            Blue points: Control cells. Colored points: Perturbed cells colored by local density difference
            (red = higher density in perturbed, blue = higher in control). Black dots mark regions with
            statistically significant differences.
        </div>
        {diff_3d_fig.to_html(full_html=False, include_plotlyjs='cdn')}
    </div>
    
    <div class="plot-container">
        <div class="plot-title">2D Density Difference Map</div>
        <div class="plot-description">
            Heatmap showing where perturbed cell density differs from control. Red areas have higher
            perturbed density, blue areas have higher control density. Contour lines show statistically
            significant regions (p < 0.05).
        </div>
        {diff_2d_fig.to_html(full_html=False, include_plotlyjs='cdn')}
    </div>
    
    <div class="plot-container">
        <div class="plot-title">Perturbation Vector Field</div>
        <div class="plot-description">
            Vectors show the average direction and magnitude of change from the control centroid
            (blue cross) for each perturbation type.
        </div>
        {vector_fig.to_html(full_html=False, include_plotlyjs='cdn')}
    </div>
</body>
</html>
"""

with open(os.path.join(OUTPUT_DIR, 'manifold_comparison_report.html'), 'w') as f:
    f.write(html_content)

# ==================== SAVE RESULTS ====================
print("Saving results...")
np.save(os.path.join(OUTPUT_DIR, 'control_umap.npy'), control_umap)
np.save(os.path.join(OUTPUT_DIR, 'perturbed_umap.npy'), perturbed_umap)

with open(os.path.join(OUTPUT_DIR, 'pca_model.pkl'), 'wb') as f:
    pickle.dump(pca_model, f)

with open(os.path.join(OUTPUT_DIR, 'umap_model.pkl'), 'wb') as f:
    pickle.dump(umap_model, f)

print(f"Analysis complete. Results saved to {OUTPUT_DIR}")