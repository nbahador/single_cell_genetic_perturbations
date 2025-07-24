import torch
import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from umap import UMAP
from torch.utils.data import Dataset, DataLoader
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
from tqdm import tqdm
import warnings
import anndata
import pickle
import html
from datetime import datetime

warnings.filterwarnings('ignore')

# Configuration
DATA_PATH = "adata_Training.h5ad"
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
OUTPUT_DIR = "manifold_visualizations"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==================== DATA LOADING AND PREPROCESSING ====================
print("Loading and preprocessing data...")

def load_and_preprocess_data():
    adata = sc.read_h5ad(DATA_PATH)
    
    # Basic preprocessing
    if hasattr(adata.X, 'toarray'):
        expression_matrix = adata.X.toarray()
    else:
        expression_matrix = adata.X
    
    # Create control mask
    control_mask = adata.obs['target_gene'] == 'non-targeting'
    control_data = expression_matrix[control_mask]
    
    return adata, expression_matrix, control_data, adata.var_names

adata, expression_matrix, control_data, gene_names = load_and_preprocess_data()

# ==================== EMBEDDING GENERATION ====================
print("Generating embeddings...")

def generate_embeddings(expression_matrix, control_data, gene_names):
    # Create gene embeddings using PCA
    pca = PCA(n_components=50)
    gene_embeddings = pca.fit_transform(expression_matrix.T)  # Genes as features
    
    # Create perturbation embeddings by averaging gene embeddings for each perturbation
    perturbations = adata.obs['target_gene'].unique()
    pert_to_idx = {pert: idx for idx, pert in enumerate(perturbations)}
    pert_embeddings = np.zeros((len(perturbations), gene_embeddings.shape[1]))
    pert_counts = np.zeros(len(perturbations))
    
    for i, row in adata.obs.iterrows():
        pert_idx = pert_to_idx[row['target_gene']]
        if row['target_gene'] in adata.var_names:
            gene_idx = adata.var_names.get_loc(row['target_gene'])
            pert_embeddings[pert_idx] += gene_embeddings[gene_idx]
        else:
            pert_embeddings[pert_idx] += np.mean(gene_embeddings, axis=0)
        pert_counts[pert_idx] += 1
    
    # Average the embeddings
    pert_embeddings = pert_embeddings / pert_counts[:, np.newaxis]
    
    return gene_embeddings, pert_embeddings, perturbations

gene_embeddings, pert_embeddings, perturbations = generate_embeddings(expression_matrix, control_data, gene_names)

# ==================== 3D MANIFOLD VISUALIZATION ====================
print("Creating 3D manifold visualizations...")

def create_3d_manifold(embeddings, labels, title, color_map='viridis'):
    # Reduce to 3D using UMAP for non-linear manifold
    reducer = UMAP(n_components=3, random_state=42)
    emb_3d = reducer.fit_transform(embeddings)
    
    # Create interactive 3D plot
    fig = go.Figure()
    
    # Add main scatter plot
    fig.add_trace(go.Scatter3d(
        x=emb_3d[:, 0],
        y=emb_3d[:, 1],
        z=emb_3d[:, 2],
        mode='markers',
        marker=dict(
            size=4,
            color=np.arange(len(emb_3d)),  # Color by index to show spread
            colorscale=color_map,
            opacity=0.8
        ),
        text=labels,
        hoverinfo='text',
        name='Points'
    ))
    
    # Add lines to show manifold structure (simplified)
    if len(emb_3d) > 5:  # Only add connections if we have enough points
        knn = NearestNeighbors(n_neighbors=min(5, len(emb_3d)-1))
        knn.fit(emb_3d)
        knn_indices = knn.kneighbors(emb_3d, return_distance=False)
        
        for i in range(min(200, len(emb_3d))):  # Limit number of lines for performance
            for j in knn_indices[i]:
                if i != j:  # Don't connect point to itself
                    fig.add_trace(go.Scatter3d(
                        x=[emb_3d[i, 0], emb_3d[j, 0]],
                        y=[emb_3d[i, 1], emb_3d[j, 1]],
                        z=[emb_3d[i, 2], emb_3d[j, 2]],
                        mode='lines',
                        line=dict(width=1, color='rgba(150,150,150,0.2)'),
                        hoverinfo='none',
                        showlegend=False
                    ))
    
    # Layout settings
    fig.update_layout(
        title=dict(
            text=f'<b>{title}</b><br><sup>3D UMAP Manifold Visualization</sup>',
            x=0.5,
            y=0.95,
            xanchor='center',
            yanchor='top'
        ),
        scene=dict(
            xaxis_title='UMAP 1',
            yaxis_title='UMAP 2',
            zaxis_title='UMAP 3',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=0.8)
            )
        ),
        margin=dict(l=0, r=0, b=0, t=40),
        height=800,
        width=1000
    )
    
    return fig, emb_3d

# Create gene manifold visualization
gene_fig, gene_3d = create_3d_manifold(
    gene_embeddings[:2000],  # Limit to top 2000 genes for performance
    gene_names[:2000],
    'Gene Expression Manifold',
    'plasma'
)

# Create perturbation manifold visualization
pert_fig, pert_3d = create_3d_manifold(
    pert_embeddings,
    perturbations,
    'Perturbation Effect Manifold',
    'rainbow'
)

# ==================== COMBINED VISUALIZATION ====================
print("Creating combined visualization...")

def create_combined_visualization(gene_3d, gene_labels, pert_3d, pert_labels):
    # Create subplots
    fig = make_subplots(
        rows=1, cols=2,
        specs=[[{'type': 'scatter3d'}, {'type': 'scatter3d'}]],
        subplot_titles=('Gene Expression Manifold', 'Perturbation Effect Manifold')
    )
    
    # Add gene plot
    fig.add_trace(go.Scatter3d(
        x=gene_3d[:, 0],
        y=gene_3d[:, 1],
        z=gene_3d[:, 2],
        mode='markers',
        marker=dict(
            size=3,
            color=np.arange(len(gene_3d)),
            colorscale='plasma',
            opacity=0.7
        ),
        text=gene_labels,
        hoverinfo='text',
        name='Genes'
    ), row=1, col=1)
    
    # Add perturbation plot
    fig.add_trace(go.Scatter3d(
        x=pert_3d[:, 0],
        y=pert_3d[:, 1],
        z=pert_3d[:, 2],
        mode='markers',
        marker=dict(
            size=5,
            color=np.arange(len(pert_3d)),
            colorscale='rainbow',
            opacity=0.8
        ),
        text=pert_labels,
        hoverinfo='text',
        name='Perturbations'
    ), row=1, col=2)
    
    # Update layout
    fig.update_layout(
        title=dict(
            text='<b>3D Manifold Visualization of Gene and Perturbation Spaces</b>',
            x=0.5,
            y=0.95,
            xanchor='center',
            yanchor='top'
        ),
        height=600,
        width=1200,
        showlegend=False
    )
    
    # Update scene properties
    fig.update_scenes(
        aspectmode='cube',
        camera=dict(eye=dict(x=1.5, y=1.5, z=0.8)),
        row=1, col=1
    )
    fig.update_scenes(
        aspectmode='cube',
        camera=dict(eye=dict(x=1.5, y=1.5, z=0.8)),
        row=1, col=2
    )
    
    return fig

combined_fig = create_combined_visualization(
    gene_3d[:1000], gene_names[:1000],  # Further reduce for combined view
    pert_3d, perturbations
)

# ==================== HTML OUTPUT GENERATION ====================
print("Generating HTML output...")

def generate_html_output(gene_fig, pert_fig, combined_fig):
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>3D Manifold Visualizations</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 0;
                padding: 20px;
                background-color: #f5f5f5;
            }}
            .container {{
                max-width: 1200px;
                margin: 0 auto;
                background-color: white;
                padding: 20px;
                border-radius: 8px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }}
            h1 {{
                color: #2c3e50;
                text-align: center;
                margin-bottom: 30px;
            }}
            .plot-container {{
                margin-bottom: 40px;
                border: 1px solid #e0e0e0;
                border-radius: 5px;
                padding: 15px;
                background-color: white;
            }}
            .plot-title {{
                font-size: 18px;
                font-weight: bold;
                margin-bottom: 15px;
                color: #3498db;
            }}
            .description {{
                margin-bottom: 15px;
                line-height: 1.6;
                color: #555;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>3D Manifold Visualizations of Gene and Perturbation Spaces</h1>
            
            <div class="plot-container">
                <div class="plot-title">Gene Expression Manifold</div>
                <div class="description">
                    This visualization shows the 3D manifold of gene expression patterns reduced using UMAP.
                    Each point represents a gene, and the spatial relationships show similarity in expression patterns.
                </div>
                {gene_fig.to_html(full_html=False, include_plotlyjs='cdn')}
            </div>
            
            <div class="plot-container">
                <div class="plot-title">Perturbation Effect Manifold</div>
                <div class="description">
                    This visualization shows the 3D manifold of perturbation effects reduced using UMAP.
                    Each point represents a perturbation, and the spatial relationships show similarity in their effects.
                </div>
                {pert_fig.to_html(full_html=False, include_plotlyjs='cdn')}
            </div>
            
            <div class="plot-container">
                <div class="plot-title">Combined View</div>
                <div class="description">
                    Side-by-side comparison of the gene and perturbation manifolds.
                    Note the different scales and structures between the two spaces.
                </div>
                {combined_fig.to_html(full_html=False, include_plotlyjs='cdn')}
            </div>
            
            <div style="margin-top: 30px; font-size: 12px; color: #777; text-align: center;">
                Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
            </div>
        </div>
    </body>
    </html>
    """
    
    output_path = os.path.join(OUTPUT_DIR, "manifold_visualizations.html")
    with open(output_path, 'w') as f:
        f.write(html_content)
    
    return output_path

html_path = generate_html_output(gene_fig, pert_fig, combined_fig)
print(f"HTML visualization saved to: {html_path}")

# ==================== ADDITIONAL VISUALIZATIONS ====================
print("Creating additional visualizations...")

def create_curved_surface_visualization(embeddings_3d, labels, title):
    # Create a more sophisticated curved surface visualization
    fig = go.Figure()
    
    # Add the main scatter plot
    fig.add_trace(go.Scatter3d(
        x=embeddings_3d[:, 0],
        y=embeddings_3d[:, 1],
        z=embeddings_3d[:, 2],
        mode='markers',
        marker=dict(
            size=5,
            color=np.linalg.norm(embeddings_3d, axis=1),
            colorscale='viridis',
            opacity=0.8,
            colorbar=dict(title='Norm')
        ),
        text=labels,
        hoverinfo='text',
        name='Points'
    ))
    
    # Add a surface to approximate the manifold
    try:
        from scipy.interpolate import griddata
        
        # Create grid for surface
        x = embeddings_3d[:, 0]
        y = embeddings_3d[:, 1]
        z = embeddings_3d[:, 2]
        
        xi = np.linspace(x.min(), x.max(), 30)
        yi = np.linspace(y.min(), y.max(), 30)
        xi, yi = np.meshgrid(xi, yi)
        
        # Interpolate z values
        zi = griddata((x, y), z, (xi, yi), method='cubic')
        
        fig.add_trace(go.Surface(
            x=xi,
            y=yi,
            z=zi,
            colorscale='viridis',
            opacity=0.3,
            showscale=False,
            name='Manifold Surface'
        ))
    except Exception as e:
        print(f"Could not create surface: {str(e)}")
    
    # Add lines to show manifold connections
    for i in np.random.choice(len(embeddings_3d), size=min(50, len(embeddings_3d)), replace=False):
        dists = np.linalg.norm(embeddings_3d - embeddings_3d[i], axis=1)
        closest = np.argsort(dists)[1:4]  # Get 3 closest points
        
        for j in closest:
            fig.add_trace(go.Scatter3d(
                x=[embeddings_3d[i, 0], embeddings_3d[j, 0]],
                y=[embeddings_3d[i, 1], embeddings_3d[j, 1]],
                z=[embeddings_3d[i, 2], embeddings_3d[j, 2]],
                mode='lines',
                line=dict(width=1, color='rgba(255,100,100,0.4)'),
                hoverinfo='none',
                showlegend=False
            ))
    
    fig.update_layout(
        title=dict(
            text=f'<b>{title}</b><br><sup>Curved Manifold with Surface Approximation</sup>',
            x=0.5,
            y=0.95
        ),
        scene=dict(
            xaxis_title='UMAP 1',
            yaxis_title='UMAP 2',
            zaxis_title='UMAP 3',
            camera=dict(
                eye=dict(x=1.8, y=1.8, z=0.6)
            )
        ),
        margin=dict(l=0, r=0, b=0, t=60),
        height=700
    )
    
    return fig

# Create curved surface visualization for genes
gene_surface_fig = create_curved_surface_visualization(
    gene_3d[:1000],  # Use subset for better performance
    gene_names[:1000],
    'Gene Expression Curved Manifold'
)

# Create curved surface visualization for perturbations
pert_surface_fig = create_curved_surface_visualization(
    pert_3d,
    perturbations,
    'Perturbation Effect Curved Manifold'
)

# Save additional visualizations
gene_surface_path = os.path.join(OUTPUT_DIR, "gene_curved_manifold.html")
pert_surface_path = os.path.join(OUTPUT_DIR, "perturbation_curved_manifold.html")

gene_surface_fig.write_html(gene_surface_path)
pert_surface_fig.write_html(pert_surface_path)

print(f"Gene curved manifold saved to: {gene_surface_path}")
print(f"Perturbation curved manifold saved to: {pert_surface_path}")

print("\nAll visualizations generated successfully!")