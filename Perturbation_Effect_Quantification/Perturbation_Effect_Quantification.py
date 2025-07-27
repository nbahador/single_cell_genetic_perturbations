import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA, IncrementalPCA
from sklearn.manifold import TSNE
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import os
from tqdm import tqdm
import warnings
from datetime import datetime
from sklearn.cluster import MiniBatchKMeans, DBSCAN
from sklearn.preprocessing import StandardScaler
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.feature_selection import SelectKBest, f_classif
import seaborn as sns
import umap

warnings.filterwarnings('ignore')

# Configuration
DATA_PATH = "adata_Training.h5ad"
OUTPUT_DIR = "improved_cellular_state_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==================== DATA LOADING AND PREPROCESSING ====================
print("Loading and preprocessing data...")

def load_and_preprocess_data(data_path, max_cells=None):
    """Load and preprocess data with memory efficiency"""
    # Read data in chunks if needed
    adata = sc.read_h5ad(data_path)
    
    # Get expression matrix
    if hasattr(adata.X, 'toarray'):
        if max_cells is not None:
            expression_matrix = adata[:max_cells].X.toarray()
            perturbation_targets = adata.obs['target_gene'].values[:max_cells]
        else:
            expression_matrix = adata.X.toarray()
            perturbation_targets = adata.obs['target_gene'].values
    else:
        if max_cells is not None:
            expression_matrix = adata[:max_cells].X
            perturbation_targets = adata.obs['target_gene'].values[:max_cells]
        else:
            expression_matrix = adata.X
            perturbation_targets = adata.obs['target_gene'].values
    
    gene_names = adata.var_names
    
    # Basic filtering
    print(f"Original data shape: {expression_matrix.shape}")
    
    # Remove genes with zero expression in all cells
    non_zero_genes = np.array(expression_matrix.sum(axis=0) > 0).flatten()
    expression_matrix = expression_matrix[:, non_zero_genes]
    gene_names = gene_names[non_zero_genes]
    
    print(f"After filtering zero genes: {expression_matrix.shape}")
    
    return expression_matrix, gene_names, perturbation_targets

# Load data with a reasonable subset if needed
try:
    expression_matrix, gene_names, perturbation_targets = load_and_preprocess_data(DATA_PATH)
except MemoryError:
    print("Memory error with full dataset, using subset of 50,000 cells")
    expression_matrix, gene_names, perturbation_targets = load_and_preprocess_data(DATA_PATH, max_cells=50000)

# Normalize data
print("Normalizing data...")
scaler = StandardScaler()
expression_matrix = scaler.fit_transform(expression_matrix)

# ==================== MEMORY-EFFICIENT CELLULAR STATE ANALYSIS ====================
class ImprovedCellularStateAnalyzer:
    def __init__(self, data, gene_names, perturbation_targets, output_dir, max_cells_subset=30000):
        self.full_data = data
        self.gene_names = gene_names
        self.perturbation_targets = perturbation_targets
        self.output_dir = output_dir
        self.max_cells_subset = max_cells_subset
        os.makedirs(output_dir, exist_ok=True)
        
        # Initialize analysis results
        self.embedding = None
        self.clusters = None
        self.perturbation_effects = None
        self.selected_genes = None
        self.data_subset = None
        self.subset_indices = None
        self.control_center = None
        
    def analyze_cellular_states(self):
        """Perform complete cellular state analysis with memory efficiency"""
        print("\nAnalyzing cellular state dynamics (memory-efficient approach)...")
        
        # 1. Feature selection and data subsetting
        self._prepare_data_for_analysis()
        
        # 2. Create embedding of cellular states
        self._create_cellular_state_embedding()
        
        # 3. Identify cellular clusters/states
        self._identify_cellular_clusters()
        
        # 4. Analyze perturbation effects
        self._analyze_perturbation_effects()
        
        # 5. Create comprehensive visualizations
        self._create_comprehensive_visualizations()
        
    def _prepare_data_for_analysis(self):
        """Prepare data with feature selection and intelligent subsampling"""
        print("\nPreparing data for analysis...")
        
        # 1. Feature selection - select most informative genes
        print("Selecting most informative genes...")
        
        # First reduce number of genes by variance
        gene_vars = np.var(self.full_data, axis=0)
        top_var_indices = np.argsort(gene_vars)[-5000:]  # Top 5000 most variable genes
        self.full_data = self.full_data[:, top_var_indices]
        self.gene_names = self.gene_names[top_var_indices]
        
        # Then use statistical test for final selection
        selector = SelectKBest(score_func=f_classif, k=min(1000, self.full_data.shape[1]))
        
        # Create binary labels for perturbation vs control
        y_binary = (self.perturbation_targets != 'non-targeting').astype(int)
        selected_features = selector.fit_transform(self.full_data, y_binary)
        selected_gene_indices = selector.get_support(indices=True)
        self.selected_genes = self.gene_names[selected_gene_indices]
        
        print(f"Selected {len(self.selected_genes)} most informative genes")
        
        # 2. Intelligent subsampling - ensure representation from all perturbations
        print("Creating representative subset...")
        subset_indices = []
        
        # Sample from each perturbation group
        unique_perts = np.unique(self.perturbation_targets)
        cells_per_pert = min(1000, self.max_cells_subset // len(unique_perts))
        
        for pert in unique_perts:
            pert_indices = np.where(self.perturbation_targets == pert)[0]
            if len(pert_indices) > cells_per_pert:
                # Stratified sampling to maintain distribution
                if len(pert_indices) > 2 * cells_per_pert:
                    # Use k-means to find representative samples
                    kmeans = MiniBatchKMeans(n_clusters=cells_per_pert, random_state=42)
                    kmeans.fit(self.full_data[pert_indices])
                    _, sampled_indices = NearestNeighbors(n_neighbors=1).fit(kmeans.cluster_centers_).kneighbors(self.full_data[pert_indices])
                    sampled_indices = np.unique(sampled_indices)
                    sampled = pert_indices[sampled_indices]
                else:
                    # Simple random sampling
                    sampled = np.random.choice(pert_indices, cells_per_pert, replace=False)
            else:
                sampled = pert_indices
            subset_indices.extend(sampled)
        
        self.subset_indices = np.array(subset_indices)
        self.data_subset = selected_features[self.subset_indices]
        
        # Normalize the subset
        scaler = StandardScaler()
        self.data_subset = scaler.fit_transform(self.data_subset)
        
        print(f"Created subset with {len(self.subset_indices)} cells and {self.data_subset.shape[1]} genes")
        
    def _create_cellular_state_embedding(self):
        """Create embedding using PCA + UMAP with memory efficiency"""
        print("\nCreating cellular state embedding...")
        
        # First, reduce dimensions with incremental PCA (memory efficient)
        print("Performing incremental PCA dimensionality reduction...")
        n_components = min(50, self.data_subset.shape[1], self.data_subset.shape[0])
        pca = IncrementalPCA(n_components=n_components, batch_size=1000)
        pca_result = pca.fit_transform(self.data_subset)
        
        print(f"PCA explained variance ratio (top 5): {pca.explained_variance_ratio_[:5]}")
        print(f"Total explained variance: {sum(pca.explained_variance_ratio_):.2f}")
        
        # Then use UMAP for final embedding with conservative parameters
        print("Creating UMAP embedding...")
        reducer = umap.UMAP(
            n_components=2,
            n_neighbors=min(15, len(self.data_subset)//10),
            min_dist=0.1,
            metric='cosine',
            random_state=42,
            verbose=True
        )
        self.embedding = reducer.fit_transform(pca_result)
        
        # Store PCA results for additional analysis
        self.pca_embedding = pca_result
        self.pca_components = pca.components_
        self.pca_variance_ratio = pca.explained_variance_ratio_
        
    def _identify_cellular_clusters(self):
        """Identify cellular clusters using density-based clustering"""
        print("\nIdentifying cellular clusters...")
        
        # Use HDBSCAN for robust clustering (handles varying densities)
        try:
            import hdbscan
            clusterer = hdbscan.HDBSCAN(
                min_cluster_size=50,
                min_samples=10,
                cluster_selection_method='eom'
            )
            self.clusters = clusterer.fit_predict(self.embedding)
        except ImportError:
            print("HDBSCAN not available, using DBSCAN instead")
            clusterer = DBSCAN(eps=1.5, min_samples=20)
            self.clusters = clusterer.fit_predict(self.embedding)
        
        # Handle noise points
        n_clusters = len(set(self.clusters)) - (1 if -1 in self.clusters else 0)
        n_noise = list(self.clusters).count(-1)
        
        print(f"Found {n_clusters} clusters with {n_noise} noise points")
        
        # Calculate cluster centers and properties
        self.cluster_centers = {}
        self.cluster_properties = {}
        
        for cluster_id in set(self.clusters):
            if cluster_id == -1:  # Skip noise
                continue
                
            mask = self.clusters == cluster_id
            cluster_points = self.embedding[mask]
            cluster_data = self.data_subset[mask]
            
            center = np.median(cluster_points, axis=0)  # More robust than mean
            size = np.sum(mask)
            
            # Calculate cluster-specific gene expression signatures
            median_expression = np.median(cluster_data, axis=0)
            
            self.cluster_centers[cluster_id] = center
            self.cluster_properties[cluster_id] = {
                'size': size,
                'median_expression': median_expression,
                'perturbations': self.perturbation_targets[self.subset_indices][mask],
                'points': cluster_points
            }
    
    def _analyze_perturbation_effects(self):
        """Analyze effects of perturbations on cellular states"""
        print("\nAnalyzing perturbation effects...")
        
        # Get control center (median is more robust)
        control_mask = self.perturbation_targets[self.subset_indices] == 'non-targeting'
        if np.sum(control_mask) == 0:
            print("Warning: No control cells found in subset")
            return
            
        self.control_center = np.median(self.embedding[control_mask], axis=0)
        
        # Calculate perturbation effects
        self.perturbation_effects = {}
        perturbation_stats = {}
        
        unique_perts = np.unique(self.perturbation_targets[self.subset_indices])
        
        for pert in tqdm(unique_perts, desc="Analyzing perturbations"):
            if pert == 'non-targeting':
                continue
                
            pert_mask = self.perturbation_targets[self.subset_indices] == pert
            if np.sum(pert_mask) < 5:  # Skip perturbations with too few cells
                continue
                
            pert_points = self.embedding[pert_mask]
            pert_center = np.median(pert_points, axis=0)
            
            # Calculate displacement vector
            displacement = pert_center - self.control_center
            displacement_magnitude = np.linalg.norm(displacement)
            
            # Calculate dispersion (how spread out the perturbed cells are)
            distances_to_center = np.linalg.norm(pert_points - pert_center, axis=1)
            dispersion = np.median(distances_to_center)
            
            # Calculate cluster distribution
            pert_clusters = self.clusters[pert_mask]
            cluster_dist = {}
            for cluster_id in set(pert_clusters):
                if cluster_id != -1:
                    cluster_dist[cluster_id] = np.sum(pert_clusters == cluster_id) / len(pert_clusters)
            
            # Calculate perturbation signature (differential expression)
            pert_data = self.data_subset[pert_mask]
            control_data = self.data_subset[control_mask][:len(pert_data)]  # Match sizes
            diff_expression = np.median(pert_data, axis=0) - np.median(control_data, axis=0)
            
            self.perturbation_effects[pert] = {
                'displacement_vector': displacement,
                'displacement_magnitude': displacement_magnitude,
                'dispersion': dispersion,
                'cluster_distribution': cluster_dist,
                'n_cells': np.sum(pert_mask),
                'differential_expression': diff_expression
            }
        
        print(f"Analyzed {len(self.perturbation_effects)} perturbations")
    
    def _create_comprehensive_visualizations(self):
        """Create comprehensive and interpretable visualizations"""
        print("\nCreating comprehensive visualizations...")
        
        # 1. Main cellular state landscape
        self._visualize_cellular_landscape()
        
        # 2. Perturbation effect analysis
        self._visualize_perturbation_effects()
        
        # 3. Cluster analysis
        self._visualize_cluster_analysis()
        
        # 4. Gene expression analysis
        self._visualize_gene_expression()
        
        # 5. Create summary dashboard
        self._create_summary_dashboard()
        
        print("\nAll visualizations saved to:", self.output_dir)
    
    def _visualize_cellular_landscape(self):
        """Create main cellular landscape visualization with clear annotations"""
        print("Creating cellular landscape visualization...")
        
        # Create figure with multiple views
        fig = make_subplots(
            rows=2, cols=2,
            specs=[[{"type": "scatter"}, {"type": "scatter"}],
                   [{"type": "scatter"}, {"type": "bar"}]],
            subplot_titles=(
                "Cellular State Landscape by Cluster",
                "Perturbation Effects Overview",
                "Cluster Characteristics",
                "Top Perturbations by Effect Size"
            ),
            horizontal_spacing=0.1,
            vertical_spacing=0.15
        )
        
        # 1. Cells colored by cluster (with density)
        cluster_ids = sorted([cid for cid in self.cluster_properties.keys() if cid != -1])
        cluster_colors = px.colors.qualitative.Plotly
        
        for cluster_id in cluster_ids:
            props = self.cluster_properties[cluster_id]
            points = props['points']
            
            # Add density contour
            x, y = points[:, 0], points[:, 1]
            kde = gaussian_kde([x, y])
            x_grid = np.linspace(min(x), max(x), 50)
            y_grid = np.linspace(min(y), max(y), 50)
            X, Y = np.meshgrid(x_grid, y_grid)
            Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
            
            fig.add_trace(
                go.Contour(
                    x=x_grid,
                    y=y_grid,
                    z=Z,
                    colorscale=[(0, 'rgba(0,0,0,0)'), (1, cluster_colors[cluster_id % len(cluster_colors)])],
                    showscale=False,
                    opacity=0.3,
                    name=f'Cluster {cluster_id} Density',
                    hoverinfo='skip'
                ),
                row=1, col=1
            )
            
            # Add points with gene expression information
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode='markers',
                    marker=dict(
                        size=4,
                        color=cluster_colors[cluster_id % len(cluster_colors)],
                        opacity=0.7
                    ),
                    name=f'Cluster {cluster_id}',
                    text=[self._get_cell_hover_text(i, cluster_id) for i in range(len(x))],
                    hoverinfo='text'
                ),
                row=1, col=1
            )
        
        # Add cluster centers
        for cluster_id, center in self.cluster_centers.items():
            fig.add_trace(
                go.Scatter(
                    x=[center[0]],
                    y=[center[1]],
                    mode='markers',
                    marker=dict(
                        size=12,
                        color='black',
                        symbol='diamond',
                        line=dict(width=2, color='white')
                    ),
                    name=f'Cluster {cluster_id} Center',
                    hoverinfo='text',
                    text=[f"Cluster {cluster_id} Center<br>Median Expression: {self._get_top_genes_text(cluster_id)}"],
                    showlegend=False
                ),
                row=1, col=1
            )
        
        # 2. Perturbation effects overview
        if self.perturbation_effects:
            # Show control center
            fig.add_trace(
                go.Scatter(
                    x=[self.control_center[0]],
                    y=[self.control_center[1]],
                    mode='markers',
                    marker=dict(
                        size=15,
                        color='red',
                        symbol='star'
                    ),
                    name='Control Center',
                    hoverinfo='text',
                    text=['Control Center (non-targeting perturbations)'],
                    showlegend=False
                ),
                row=1, col=2
            )
            
            # Show top 10 perturbations by effect size
            top_perts = sorted(self.perturbation_effects.items(), 
                              key=lambda x: x[1]['displacement_magnitude'], 
                              reverse=True)[:10]
            
            for pert, effects in top_perts:
                end_point = self.control_center + effects['displacement_vector']
                
                # Add arrow
                fig.add_trace(
                    go.Scatter(
                        x=[self.control_center[0], end_point[0]],
                        y=[self.control_center[1], end_point[1]],
                        mode='lines',
                        line=dict(width=2, color='orange'),
                        name=f'{pert[:15]}...',
                        hoverinfo='text',
                        text=[f"{pert}<br>Effect size: {effects['displacement_magnitude']:.2f}<br>Top affected genes: {self._get_top_affected_genes_text(pert)}"],
                        showlegend=False
                    ),
                    row=1, col=2
                )
                
                # Add end point
                fig.add_trace(
                    go.Scatter(
                        x=[end_point[0]],
                        y=[end_point[1]],
                        mode='markers',
                        marker=dict(
                            size=8,
                            color='blue',
                            symbol='circle'
                        ),
                        name=f'{pert[:15]}...',
                        hoverinfo='text',
                        text=[f"{pert}<br>Effect size: {effects['displacement_magnitude']:.2f}<br>Top affected genes: {self._get_top_affected_genes_text(pert)}"],
                        showlegend=False
                    ),
                    row=1, col=2
                )
        
        # 3. Cluster characteristics
        cluster_stats = []
        for cluster_id, props in self.cluster_properties.items():
            perts = props['perturbations']
            unique_perts, counts = np.unique(perts, return_counts=True)
            dominant_pert = unique_perts[np.argmax(counts)]
            dominant_frac = counts.max() / len(perts)
            
            cluster_stats.append({
                'Cluster': cluster_id,
                'Size': props['size'],
                'Dominant Perturbation': dominant_pert,
                'Dominant Fraction': dominant_frac,
                'Unique Perturbations': len(unique_perts)
            })
        
        cluster_df = pd.DataFrame(cluster_stats)
        
        # Add size vs diversity
        fig.add_trace(
            go.Scatter(
                x=cluster_df['Size'],
                y=cluster_df['Unique Perturbations'],
                mode='markers+text',
                marker=dict(
                    size=cluster_df['Dominant Fraction'] * 20 + 10,
                    color=cluster_df['Cluster'],
                    colorscale='Viridis',
                    showscale=True,
                    colorbar=dict(title='Cluster ID')
                ),
                text=[f"Cluster {cid}" for cid in cluster_df['Cluster']],
                textposition='top center',
                name='Clusters',
                hovertext=[
                    f"Cluster {row['Cluster']}<br>"
                    f"Size: {row['Size']}<br>"
                    f"Dominant Pert: {row['Dominant Perturbation']}<br>"
                    f"Dominant Fraction: {row['Dominant Fraction']:.2f}<br>"
                    f"Unique Perturbations: {row['Unique Perturbations']}<br>"
                    f"Top Expressed Genes: {self._get_top_genes_text(row['Cluster'])}"
                    for _, row in cluster_df.iterrows()
                ],
                hoverinfo='text'
            ),
            row=2, col=1
        )
        
        # 4. Top perturbations by effect size
        if self.perturbation_effects:
            top_perts = sorted(self.perturbation_effects.items(), 
                             key=lambda x: x[1]['displacement_magnitude'], 
                             reverse=True)[:10]
            
            pert_names = [p[0][:15] + '...' if len(p[0]) > 15 else p[0] for p in top_perts]
            effect_sizes = [p[1]['displacement_magnitude'] for p in top_perts]
            
            fig.add_trace(
                go.Bar(
                    x=pert_names,
                    y=effect_sizes,
                    marker_color=px.colors.sequential.Plasma,
                    name='Effect Size',
                    hovertext=[
                        f"{p[0]}<br>"
                        f"Effect Size: {p[1]['displacement_magnitude']:.2f}<br>"
                        f"Cells: {p[1]['n_cells']}<br>"
                        f"Dispersion: {p[1]['dispersion']:.2f}<br>"
                        f"Top Affected Genes: {self._get_top_affected_genes_text(p[0])}"
                        for p in top_perts
                    ],
                    hoverinfo='text'
                ),
                row=2, col=2
            )
        
        # Update layout
        fig.update_layout(
            title_text="Comprehensive Cellular State Analysis",
            height=900,
            showlegend=True,
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            margin=dict(l=50, r=50, b=50, t=100, pad=4)
        )
        
        # Update axes
        fig.update_xaxes(title_text="UMAP 1", row=1, col=1)
        fig.update_yaxes(title_text="UMAP 2", row=1, col=1)
        fig.update_xaxes(title_text="UMAP 1", row=1, col=2)
        fig.update_yaxes(title_text="UMAP 2", row=1, col=2)
        fig.update_xaxes(title_text="Cluster Size", row=2, col=1)
        fig.update_yaxes(title_text="Unique Perturbations", row=2, col=1)
        fig.update_xaxes(title_text="Perturbation", row=2, col=2)
        fig.update_yaxes(title_text="Effect Size", row=2, col=2)
        
        filename = os.path.join(self.output_dir, "cellular_landscape.html")
        fig.write_html(filename)
        print(f"Saved cellular landscape to: {filename}")
    
    def _get_cell_hover_text(self, idx, cluster_id):
        """Generate hover text for a cell with gene expression information"""
        cell_idx = self.subset_indices[idx]
        pert = self.perturbation_targets[cell_idx]
        top_genes_idx = np.argsort(np.abs(self.data_subset[idx]))[-5:][::-1]  # Top 5 most expressed genes
        top_genes = self.selected_genes[top_genes_idx]
        top_values = self.data_subset[idx][top_genes_idx]
        
        text = f"Cell {idx}<br>Cluster: {cluster_id}<br>Perturbation: {pert}<br>Top Expressed Genes:"
        for gene, val in zip(top_genes, top_values):
            text += f"<br>{gene}: {val:.2f}"
        
        return text
    
    def _get_top_genes_text(self, cluster_id):
        """Generate text for top expressed genes in a cluster"""
        if cluster_id not in self.cluster_properties:
            return "N/A"
        
        median_expr = self.cluster_properties[cluster_id]['median_expression']
        top_genes_idx = np.argsort(median_expr)[-5:][::-1]  # Top 5 most expressed genes
        top_genes = self.selected_genes[top_genes_idx]
        top_values = median_expr[top_genes_idx]
        
        text = ""
        for gene, val in zip(top_genes, top_values):
            text += f"<br>{gene}: {val:.2f}"
        
        return text
    
    def _get_top_affected_genes_text(self, pert):
        """Generate text for top affected genes by a perturbation"""
        if pert not in self.perturbation_effects:
            return "N/A"
        
        diff_expr = self.perturbation_effects[pert]['differential_expression']
        top_genes_idx = np.argsort(np.abs(diff_expr))[-5:][::-1]  # Top 5 most affected genes
        top_genes = self.selected_genes[top_genes_idx]
        top_values = diff_expr[top_genes_idx]
        
        text = ""
        for gene, val in zip(top_genes, top_values):
            text += f"<br>{gene}: {val:.2f}"
        
        return text
    
    def _visualize_perturbation_effects(self):
        """Create detailed perturbation effect visualizations"""
        print("Creating perturbation effects visualization...")
        
        if not self.perturbation_effects:
            print("No perturbation effects to visualize")
            return
        
        # Prepare data for visualization
        pert_data = []
        for pert, effects in self.perturbation_effects.items():
            pert_data.append({
                'Perturbation': pert,
                'Effect Size': effects['displacement_magnitude'],
                'Dispersion': effects['dispersion'],
                'Cells': effects['n_cells'],
                'Main Cluster': max(effects['cluster_distribution'].items(), key=lambda x: x[1])[0] if effects['cluster_distribution'] else -1,
                'Main Cluster Fraction': max(effects['cluster_distribution'].values()) if effects['cluster_distribution'] else 0
            })
        
        pert_df = pd.DataFrame(pert_data)
        
        # Create figure
        fig = make_subplots(
            rows=2, cols=2,
            specs=[[{"type": "scatter"}, {"type": "scatter"}],
                   [{"type": "heatmap"}, {"type": "bar"}]],
            subplot_titles=(
                "Effect Size vs Dispersion",
                "Effect Size vs Cluster Specificity",
                "Cluster Distribution Heatmap (Top 20 Perturbations)",
                "Top Differential Genes for Selected Perturbation"
            )
        )
        
        # 1. Effect Size vs Dispersion
        fig.add_trace(
            go.Scatter(
                x=pert_df['Effect Size'],
                y=pert_df['Dispersion'],
                mode='markers',
                marker=dict(
                    size=np.log10(pert_df['Cells']) * 5,
                    color=pert_df['Main Cluster'],
                    colorscale='Viridis',
                    showscale=True,
                    colorbar=dict(title='Main Cluster')
                ),
                text=pert_df['Perturbation'],
                hoverinfo='text',
                name='Perturbations'
            ),
            row=1, col=1
        )
        
        # 2. Effect Size vs Cluster Specificity
        fig.add_trace(
            go.Scatter(
                x=pert_df['Effect Size'],
                y=pert_df['Main Cluster Fraction'],
                mode='markers',
                marker=dict(
                    size=np.log10(pert_df['Cells']) * 5,
                    color=pert_df['Dispersion'],
                    colorscale='Plasma',
                    showscale=True,
                    colorbar=dict(title='Dispersion')
                ),
                text=pert_df['Perturbation'],
                hoverinfo='text',
                name='Perturbations'
            ),
            row=1, col=2
        )
        
        # 3. Cluster distribution heatmap (top 20 perturbations)
        top_perts = pert_df.nlargest(20, 'Effect Size')
        heatmap_data = []
        heatmap_text = []
        
        for pert in top_perts['Perturbation']:
            effects = self.perturbation_effects[pert]
            row = []
            row_text = []
            for cid in sorted(self.cluster_properties.keys()):
                val = effects['cluster_distribution'].get(cid, 0)
                row.append(val)
                row_text.append(f"{pert}<br>Cluster {cid}: {val:.2f}")
            heatmap_data.append(row)
            heatmap_text.append(row_text)
        
        fig.add_trace(
            go.Heatmap(
                z=heatmap_data,
                x=[f'Cluster {cid}' for cid in sorted(self.cluster_properties.keys())],
                y=[pert[:15] + '...' for pert in top_perts['Perturbation']],
                colorscale='Blues',
                hoverinfo='text',
                text=heatmap_text,
                showscale=True
            ),
            row=2, col=1
        )
        
        # 4. Top differential genes for selected perturbation
        # Select the perturbation with largest effect size
        selected_pert = top_perts.iloc[0]['Perturbation']
        diff_genes = self.perturbation_effects[selected_pert]['differential_expression']
        top_gene_indices = np.argsort(np.abs(diff_genes))[-20:]  # Top 20 most differential genes
        top_genes = self.selected_genes[top_gene_indices]
        
        fig.add_trace(
            go.Bar(
                x=diff_genes[top_gene_indices],
                y=top_genes,
                orientation='h',
                marker_color=np.where(diff_genes[top_gene_indices] > 0, 'red', 'blue'),
                name='Differential Expression',
                hoverinfo='text',
                text=[f"{gene}<br>logFC: {val:.2f}" for gene, val in zip(top_genes, diff_genes[top_gene_indices])]
            ),
            row=2, col=2
        )
        
        # Update layout
        fig.update_layout(
            title_text=f"Detailed Perturbation Effects Analysis<br><sup>Selected Perturbation: {selected_pert}</sup>",
            height=900,
            showlegend=False
        )
        
        # Update axes
        fig.update_xaxes(title_text="Effect Size", row=1, col=1)
        fig.update_yaxes(title_text="Dispersion", row=1, col=1)
        fig.update_xaxes(title_text="Effect Size", row=1, col=2)
        fig.update_yaxes(title_text="Main Cluster Fraction", row=1, col=2)
        fig.update_xaxes(title_text="Cluster", row=2, col=1)
        fig.update_yaxes(title_text="Perturbation", row=2, col=1)
        fig.update_xaxes(title_text="Differential Expression", row=2, col=2)
        fig.update_yaxes(title_text="Gene", row=2, col=2)
        
        filename = os.path.join(self.output_dir, "perturbation_effects.html")
        fig.write_html(filename)
        print(f"Saved perturbation effects to: {filename}")
    
    def _visualize_cluster_analysis(self):
        """Create detailed cluster analysis visualizations"""
        print("Creating cluster analysis visualization...")
        
        # Prepare cluster data
        cluster_stats = []
        for cluster_id, props in self.cluster_properties.items():
            perts = props['perturbations']
            unique_perts, counts = np.unique(perts, return_counts=True)
            dominant_pert = unique_perts[np.argmax(counts)]
            dominant_frac = counts.max() / len(perts)
            
            cluster_stats.append({
                'Cluster': cluster_id,
                'Size': props['size'],
                'Dominant Perturbation': dominant_pert,
                'Dominant Fraction': dominant_frac,
                'Unique Perturbations': len(unique_perts)
            })
        
        cluster_df = pd.DataFrame(cluster_stats)
        
        # Prepare perturbation distribution data
        pert_dist_data = []
        for pert, effects in self.perturbation_effects.items():
            for cluster_id, frac in effects['cluster_distribution'].items():
                pert_dist_data.append({
                    'Perturbation': pert,
                    'Cluster': cluster_id,
                    'Fraction': frac
                })
        
        pert_dist_df = pd.DataFrame(pert_dist_data)
        
        # Create figure
        fig = make_subplots(
            rows=2, cols=2,
            specs=[[{"type": "bar"}, {"type": "pie"}],
                   [{"type": "scatter"}, {"type": "box"}]],
            subplot_titles=(
                "Cluster Sizes",
                "Perturbation Distribution Across Clusters",
                "Cluster Characteristics",
                "Perturbation Specificity by Cluster"
            )
        )
        
        # 1. Cluster sizes
        fig.add_trace(
            go.Bar(
                x=[f"Cluster {cid}" for cid in cluster_df['Cluster']],
                y=cluster_df['Size'],
                marker_color=px.colors.qualitative.Plotly[:len(cluster_df)],
                name='Cluster Size',
                hoverinfo='text',
                text=[f"Size: {size}<br>Top Genes: {self._get_top_genes_text(cid)}" for cid, size in zip(cluster_df['Cluster'], cluster_df['Size'])]
            ),
            row=1, col=1
        )
        
        # 2. Perturbation distribution pie chart (aggregate)
        fig.add_trace(
            go.Pie(
                labels=[f"Cluster {cid}" for cid in cluster_df['Cluster']],
                values=cluster_df['Size'],
                name='Cluster Distribution',
                hoverinfo='label+percent+text',
                text=[f"Top Genes: {self._get_top_genes_text(cid)}" for cid in cluster_df['Cluster']],
                textinfo='percent',
                marker_colors=px.colors.qualitative.Plotly[:len(cluster_df)]
            ),
            row=1, col=2
        )
        
        # 3. Cluster characteristics scatter
        fig.add_trace(
            go.Scatter(
                x=cluster_df['Unique Perturbations'],
                y=cluster_df['Dominant Fraction'],
                mode='markers+text',
                marker=dict(
                    size=cluster_df['Size'] / cluster_df['Size'].max() * 30 + 10,
                    color=cluster_df['Cluster'],
                    colorscale='Viridis',
                    showscale=True,
                    colorbar=dict(title='Cluster ID')
                ),
                text=[f"Cluster {cid}" for cid in cluster_df['Cluster']],
                textposition='top center',
                name='Clusters',
                hoverinfo='text',
                hovertext=[
                    f"Cluster {row['Cluster']}<br>"
                    f"Size: {row['Size']}<br>"
                    f"Dominant Pert: {row['Dominant Perturbation']}<br>"
                    f"Dominant Fraction: {row['Dominant Fraction']:.2f}<br>"
                    f"Unique Perturbations: {row['Unique Perturbations']}<br>"
                    f"Top Expressed Genes: {self._get_top_genes_text(row['Cluster'])}"
                    for _, row in cluster_df.iterrows()
                ]
            ),
            row=2, col=1
        )
        
        # 4. Perturbation specificity by cluster (box plot)
        if not pert_dist_df.empty:
            fig.add_trace(
                go.Box(
                    x=pert_dist_df['Cluster'],
                    y=pert_dist_df['Fraction'],
                    name='Perturbation Specificity',
                    boxpoints='all',
                    jitter=0.3,
                    pointpos=-1.8,
                    marker_color=px.colors.qualitative.Plotly[0],
                    hoverinfo='y',
                    hovertext=[
                        f"{row['Perturbation']}<br>"
                        f"Cluster: {row['Cluster']}<br>"
                        f"Fraction: {row['Fraction']:.2f}<br>"
                        f"Top Affected Genes: {self._get_top_affected_genes_text(row['Perturbation'])}"
                        for _, row in pert_dist_df.iterrows()
                    ]
                ),
                row=2, col=2
            )
        
        # Update layout
        fig.update_layout(
            title_text="Detailed Cluster Analysis",
            height=900,
            showlegend=False
        )
        
        # Update axes
        fig.update_xaxes(title_text="Cluster", row=1, col=1)
        fig.update_yaxes(title_text="Number of Cells", row=1, col=1)
        fig.update_xaxes(title_text="Unique Perturbations", row=2, col=1)
        fig.update_yaxes(title_text="Dominant Perturbation Fraction", row=2, col=1)
        fig.update_xaxes(title_text="Cluster", row=2, col=2)
        fig.update_yaxes(title_text="Perturbation Fraction", row=2, col=2)
        
        filename = os.path.join(self.output_dir, "cluster_analysis.html")
        fig.write_html(filename)
        print(f"Saved cluster analysis to: {filename}")
    
    def _visualize_gene_expression(self):
        """Create gene expression visualizations"""
        print("Creating gene expression visualization...")
        
        # Select top variable genes
        gene_vars = np.var(self.data_subset, axis=0)
        top_gene_indices = np.argsort(gene_vars)[-50:]  # Top 50 most variable genes
        top_genes = self.selected_genes[top_gene_indices]
        
        # Create heatmap data (sample 1000 cells)
        sample_indices = np.random.choice(len(self.data_subset), min(1000, len(self.data_subset)), replace=False)
        heatmap_data = self.data_subset[sample_indices][:, top_gene_indices]
        cell_labels = self.perturbation_targets[self.subset_indices][sample_indices]
        
        # Create figure
        fig = go.Figure()
        
        # Add heatmap with enhanced hover information
        fig.add_trace(
            go.Heatmap(
                z=heatmap_data.T,
                x=[f"Cell {i}" for i in range(len(heatmap_data))],
                y=top_genes,
                colorscale='RdBu_r',
                zmid=0,
                hoverinfo='text',
                text=[
                    [f"Gene: {gene}<br>Cell: {cell_labels[i]}<br>Expression: {val:.2f}<br>Cluster: {self.clusters[sample_indices[i]] if self.clusters is not None else 'N/A'}" 
                     for val in row] 
                    for i, row, gene in zip(range(len(heatmap_data)), heatmap_data.T, top_genes)
                ]
            )
        )
        
        # Add cluster annotation if available
        if self.clusters is not None:
            cluster_labels = self.clusters[sample_indices]
            unique_clusters = np.unique(cluster_labels)
            
            for cluster_id in unique_clusters:
                if cluster_id == -1:
                    continue
                    
                cluster_mask = cluster_labels == cluster_id
                cluster_indices = np.where(cluster_mask)[0]
                
                fig.add_trace(
                    go.Scatter(
                        x=cluster_indices,
                        y=[-1] * len(cluster_indices),
                        mode='markers',
                        marker=dict(
                            size=5,
                            color=px.colors.qualitative.Plotly[cluster_id % len(px.colors.qualitative.Plotly)]
                        ),
                        name=f'Cluster {cluster_id}',
                        hoverinfo='text',
                        text=[f"Cluster {cluster_id}<br>Top Genes: {self._get_top_genes_text(cluster_id)}"] * len(cluster_indices),
                        showlegend=True
                    )
                )
        
        # Update layout
        fig.update_layout(
            title_text="Gene Expression Heatmap (Top 50 Variable Genes)",
            xaxis_title="Cells",
            yaxis_title="Genes",
            height=800,
            margin=dict(l=100, r=50, b=100, t=100)
        )
        
        filename = os.path.join(self.output_dir, "gene_expression_heatmap.html")
        fig.write_html(filename)
        print(f"Saved gene expression heatmap to: {filename}")
    
    def _create_summary_dashboard(self):
        """Create a summary dashboard with key findings"""
        print("Creating summary dashboard...")
        
        # Prepare summary statistics
        n_clusters = len(self.cluster_properties)
        n_perturbations = len(self.perturbation_effects) if self.perturbation_effects else 0
        avg_effect_size = np.mean([eff['displacement_magnitude'] for eff in self.perturbation_effects.values()]) if self.perturbation_effects else 0
        
        # Get top clusters by size
        cluster_sizes = [(cid, props['size']) for cid, props in self.cluster_properties.items()]
        top_clusters = sorted(cluster_sizes, key=lambda x: x[1], reverse=True)[:3]
        
        # Get top perturbations by effect size
        if self.perturbation_effects:
            top_perts = sorted(self.perturbation_effects.items(), 
                             key=lambda x: x[1]['displacement_magnitude'], 
                             reverse=True)[:3]
        else:
            top_perts = []
        
        # Create dashboard
        fig = go.Figure()
        
        # Add summary text
        summary_text = [
            "<b>Cellular State Analysis Summary</b><br><br>",
            f"<b>Dataset Overview:</b>",
            f"- Analyzed {len(self.data_subset)} cells with {len(self.selected_genes)} genes",
            f"- Identified {n_clusters} distinct cellular clusters",
            f"- Analyzed {n_perturbations} perturbations with average effect size {avg_effect_size:.2f}",
            "<br>",
            "<b>Top Clusters by Size:</b>"
        ]
        
        for cid, size in top_clusters:
            summary_text.append(f"- Cluster {cid}: {size} cells<br>Top Genes: {self._get_top_genes_text(cid)}")
        
        summary_text.append("<br><b>Top Perturbations by Effect Size:</b>")
        
        for pert, effects in top_perts:
            summary_text.append(f"- {pert[:20]}...: effect size {effects['displacement_magnitude']:.2f}<br>Top Affected Genes: {self._get_top_affected_genes_text(pert)}")
        
        fig.add_trace(
            go.Scatter(
                x=[0], y=[0],
                mode='text',
                text="<br>".join(summary_text),
                textfont=dict(size=14),
                hoverinfo='none',
                showlegend=False
            )
        )
        
        # Add UMAP plot if available
        if self.embedding is not None:
            fig.add_trace(
                go.Scatter(
                    x=self.embedding[:, 0],
                    y=self.embedding[:, 1],
                    mode='markers',
                    marker=dict(
                        size=3,
                        color=self.clusters if self.clusters is not None else 'blue',
                        colorscale='Viridis',
                        opacity=0.7
                    ),
                    name='Cells',
                    hoverinfo='text',
                    text=[self._get_cell_hover_text(i, self.clusters[i] if self.clusters is not None else -1) for i in range(len(self.embedding))],
                    visible=False
                )
            )
        
        # Update layout
        fig.update_layout(
            title_text="Cellular State Analysis Dashboard",
            height=600,
            showlegend=False,
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            margin=dict(l=50, r=50, b=50, t=100, pad=4)
        )
        
        filename = os.path.join(self.output_dir, "summary_dashboard.html")
        fig.write_html(filename)
        print(f"Saved summary dashboard to: {filename}")

# ==================== MAIN EXECUTION ====================
if __name__ == "__main__":
    print("\nStarting cellular state analysis pipeline...")
    
    # Initialize analyzer
    analyzer = ImprovedCellularStateAnalyzer(
        expression_matrix,
        gene_names,
        perturbation_targets,
        OUTPUT_DIR
    )
    
    # Perform analysis
    analyzer.analyze_cellular_states()
    
    print("\nAnalysis complete! Results saved to:", OUTPUT_DIR)