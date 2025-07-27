import numpy as np
import scanpy as sc
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.decomposition import PCA
from umap import UMAP
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
from scipy.sparse import issparse
from scipy.stats import zscore
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

# Configuration
DATA_PATH = "adata_Training.h5ad"
OUTPUT_DIR = "autoencoder_trajectories"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Set device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# Set random seeds
torch.manual_seed(42)
np.random.seed(42)

# ==================== DATA LOADING ====================
def load_and_preprocess():
    print("Loading data...")
    adata = sc.read_h5ad(DATA_PATH)
    
    # Convert to dense if sparse
    X = adata.X.toarray() if issparse(adata.X) else adata.X
    
    # Basic filtering
    nonzero_genes = (X > 0).sum(axis=0) > 0
    X = X[:, nonzero_genes]
    gene_names = adata.var_names[nonzero_genes]
    
    # Log transform and standardize
    X = np.log1p(X)
    X = zscore(X, axis=0)
    
    # Subsample if needed
    if X.shape[0] > 20000:
        idx = np.random.choice(X.shape[0], 20000, replace=False)
        X = X[idx]
        adata = adata[idx].copy()
    
    return adata, X, gene_names

# ==================== IMPROVED AUTOENCODER ====================
class RobustAutoencoder(nn.Module):
    def __init__(self, input_dim, latent_dim=32):
        super().__init__()
        
        # Encoder with dropout and batch norm
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 512),
            nn.BatchNorm1d(512),
            nn.LeakyReLU(),
            nn.Dropout(0.2),
            
            nn.Linear(512, 256),
            nn.BatchNorm1d(256),
            nn.LeakyReLU(),
            nn.Dropout(0.2),
            
            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.LeakyReLU(),
            
            nn.Linear(128, latent_dim)
        )
        
        # Decoder with dropout and batch norm
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 128),
            nn.BatchNorm1d(128),
            nn.LeakyReLU(),
            nn.Dropout(0.1),
            
            nn.Linear(128, 256),
            nn.BatchNorm1d(256),
            nn.LeakyReLU(),
            nn.Dropout(0.1),
            
            nn.Linear(256, 512),
            nn.BatchNorm1d(512),
            nn.LeakyReLU(),
            
            nn.Linear(512, input_dim)
        )
        
        # Track optimization trajectory
        self.trajectory = {
            'weights': [],
            'gradients': [],
            'losses': [],
            'val_losses': [],
            'latents': [],
            'checkpoints': [],
            'epoch_data': []
        }
    
    def forward(self, x):
        z = self.encoder(x)
        x_recon = self.decoder(z)
        return z, x_recon
    
    def track_trajectory(self, epoch, X_sample, loss, val_loss, adata_sample):
        """Track optimization state at each epoch"""
        # Store weights and gradients
        weights = []
        grads = []
        for param in self.parameters():
            if param.requires_grad:
                weights.append(param.data.cpu().numpy().flatten())
                if param.grad is not None:
                    grads.append(param.grad.data.cpu().numpy().flatten())
        
        self.trajectory['weights'].append(np.concatenate(weights))
        self.trajectory['gradients'].append(np.concatenate(grads))
        self.trajectory['losses'].append(loss)
        self.trajectory['val_losses'].append(val_loss)
        
        # Store latent representations
        with torch.no_grad():
            z, _ = self.forward(X_sample)
            latent = z.cpu().numpy()
            self.trajectory['latents'].append(latent)
            
            # Get available metadata for hover info
            hover_text = []
            for i in range(len(adata_sample)):
                text_parts = []
                if 'target_gene' in adata_sample.obs:
                    text_parts.append(f"Perturbation: {adata_sample.obs['target_gene'].values[i]}")
                if 'cell_type' in adata_sample.obs and pd.notna(adata_sample.obs['cell_type'].values[i]):
                    text_parts.append(f"Cell type: {adata_sample.obs['cell_type'].values[i]}")
                if 'n_genes' in adata_sample.obs:
                    text_parts.append(f"Genes: {adata_sample.obs['n_genes'].values[i]}")
                hover_text.append("<br>".join(text_parts))
            
            # Store complete epoch data for visualization
            epoch_data = {
                'epoch': epoch,
                'latent': latent,
                'perturbations': adata_sample.obs['target_gene'].values,
                'cell_types': adata_sample.obs.get('cell_type', np.array(['']*len(adata_sample))),
                'hover_text': hover_text,
                'loss': loss,
                'val_loss': val_loss
            }
            self.trajectory['epoch_data'].append(epoch_data)
        
        self.trajectory['checkpoints'].append(epoch)

# ==================== IMPROVED TRAINING UTILITIES ====================
def train_robust(model, X, adata, batch_size=256, epochs=100, lr=1e-3):
    print("Training robust autoencoder...")
    
    # Split into train and validation
    n = X.shape[0]
    idx = np.random.permutation(n)
    train_idx, val_idx = idx[:int(0.9*n)], idx[int(0.9*n):]
    
    X_train, X_val = X[train_idx], X[val_idx]
    adata_train, adata_val = adata[train_idx].copy(), adata[val_idx].copy()
    
    # Convert to tensors
    X_train_tensor = torch.FloatTensor(X_train).to(device)
    X_val_tensor = torch.FloatTensor(X_val).to(device)
    
    train_dataset = TensorDataset(X_train_tensor)
    val_dataset = TensorDataset(X_val_tensor)
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
    
    # Optimizer with weight decay (L2 regularization)
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    criterion = nn.MSELoss()
    
    # Learning rate scheduler
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=5, factor=0.5, verbose=True)
    
    # Early stopping
    best_val_loss = np.inf
    patience = 10
    patience_counter = 0
    
    # Sample for trajectory tracking
    sample_size = min(1000, len(X_train_tensor))
    sample_idx = np.random.choice(len(X_train_tensor), sample_size, replace=False)
    X_sample = X_train_tensor[sample_idx]
    adata_sample = adata_train[sample_idx].copy()
    
    # Training loop with visualization at each epoch
    for epoch in tqdm(range(epochs)):
        model.train()
        epoch_loss = 0
        
        for batch in train_loader:
            x = batch[0]
            
            # Forward pass
            _, x_recon = model(x)
            loss = criterion(x_recon, x)
            
            # Backward pass with gradient clipping
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            
            epoch_loss += loss.item()
        
        # Validation
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for batch in val_loader:
                x = batch[0]
                _, x_recon = model(x)
                val_loss += criterion(x_recon, x).item()
        
        train_loss = epoch_loss/len(train_loader)
        val_loss = val_loss/len(val_loader)
        
        # Track trajectory and visualize
        model.track_trajectory(epoch, X_sample, train_loss, val_loss, adata_sample)
        
        # Learning rate scheduling
        scheduler.step(val_loss)
        
        # Early stopping check
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
            # Save best model
            torch.save(model.state_dict(), os.path.join(OUTPUT_DIR, "best_model.pt"))
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print(f"Early stopping at epoch {epoch}")
                break
        
        # Save interactive visualization for this epoch
        if epoch % 5 == 0 or epoch == epochs-1:
            save_interactive_visualization(model, epoch)
    
    # Load best model
    model.load_state_dict(torch.load(os.path.join(OUTPUT_DIR, "best_model.pt")))
    return model

# ==================== ENHANCED INTERACTIVE VISUALIZATION ====================
def save_interactive_visualization(model, epoch):
    """Create and save enhanced interactive Plotly visualization for this epoch"""
    if len(model.trajectory['epoch_data']) == 0:
        return
    
    # Get the data for this epoch
    epoch_data = model.trajectory['epoch_data'][-1]
    latent = epoch_data['latent']
    perturbations = epoch_data['perturbations']
    cell_types = epoch_data['cell_types']
    hover_text = epoch_data['hover_text']
    
    # Create subplots: 3D trajectory, 3D latent space, 2D projections, and loss curves
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{'type': 'scatter3d'}, {'type': 'scatter3d'}],
               [{'type': 'scatter'}, {'type': 'scatter'}]],
        subplot_titles=(
            'Parameter Space Trajectory (3D PCA)',
            'Latent Space (3D UMAP)',
            'Training & Validation Loss',
            'Perturbation Effects (2D UMAP)'
        ),
        horizontal_spacing=0.1,
        vertical_spacing=0.15
    )
    
    # 1. Parameter space trajectory (3D PCA)
    if len(model.trajectory['weights']) > 3:
        weights = np.array(model.trajectory['weights'])
        pca = PCA(n_components=3)
        weights_pca = pca.fit_transform(weights)
        
        # Create continuous color scale
        colors = plt.cm.viridis(np.linspace(0, 1, len(weights_pca)))
        
        # Main trajectory line
        fig.add_trace(go.Scatter3d(
            x=weights_pca[:, 0],
            y=weights_pca[:, 1],
            z=weights_pca[:, 2],
            mode='lines',
            line=dict(color='gray', width=4),
            name='Trajectory',
            hoverinfo='none',
            showlegend=False
        ), row=1, col=1)
        
        # Colored segments
        for i in range(1, len(weights_pca)):
            fig.add_trace(go.Scatter3d(
                x=weights_pca[i-1:i+1, 0],
                y=weights_pca[i-1:i+1, 1],
                z=weights_pca[i-1:i+1, 2],
                mode='lines',
                line=dict(color=f'rgb({colors[i][0]*255}, {colors[i][1]*255}, {colors[i][2]*255})', width=6),
                hoverinfo='none',
                showlegend=False
            ), row=1, col=1)
        
        # Add start and end markers
        fig.add_trace(go.Scatter3d(
            x=[weights_pca[0, 0]],
            y=[weights_pca[0, 1]],
            z=[weights_pca[0, 2]],
            mode='markers',
            marker=dict(size=8, color='red', symbol='diamond'),
            name='Start',
            hoverinfo='text',
            hovertext=f"Epoch 0<br>Loss: {model.trajectory['losses'][0]:.4f}"
        ), row=1, col=1)
        
        fig.add_trace(go.Scatter3d(
            x=[weights_pca[-1, 0]],
            y=[weights_pca[-1, 1]],
            z=[weights_pca[-1, 2]],
            mode='markers',
            marker=dict(size=8, color='green', symbol='diamond'),
            name='Current',
            hoverinfo='text',
            hovertext=f"Epoch {epoch}<br>Loss: {model.trajectory['losses'][-1]:.4f}"
        ), row=1, col=1)
    
    # 2. Latent space projection (3D UMAP)
    if latent.shape[0] > 0:
        reducer = UMAP(n_components=3, random_state=42)
        latent_3d = reducer.fit_transform(latent)
        
        unique_perts = np.unique(perturbations)
        for i, pert in enumerate(unique_perts):
            mask = perturbations == pert
            if np.sum(mask) > 0:
                fig.add_trace(go.Scatter3d(
                    x=latent_3d[mask, 0],
                    y=latent_3d[mask, 1],
                    z=latent_3d[mask, 2],
                    mode='markers',
                    marker=dict(
                        size=4,
                        color=px.colors.qualitative.Plotly[i % len(px.colors.qualitative.Plotly)],
                        opacity=0.7
                    ),
                    name=pert,
                    text=[hover_text[i] for i in np.where(mask)[0]],
                    hoverinfo='text',
                    showlegend=True
                ), row=1, col=2)
    
    # 3. Loss curves (training and validation)
    if len(model.trajectory['losses']) > 0:
        fig.add_trace(go.Scatter(
            x=model.trajectory['checkpoints'],
            y=model.trajectory['losses'],
            mode='lines+markers',
            name='Training Loss',
            line=dict(color='blue'),
            hoverinfo='x+y'
        ), row=2, col=1)
        
        fig.add_trace(go.Scatter(
            x=model.trajectory['checkpoints'],
            y=model.trajectory['val_losses'],
            mode='lines+markers',
            name='Validation Loss',
            line=dict(color='red'),
            hoverinfo='x+y'
        ), row=2, col=1)
    
    # 4. Perturbation effects (2D UMAP)
    if latent.shape[0] > 0:
        reducer = UMAP(n_components=2, random_state=42)
        latent_2d = reducer.fit_transform(latent)
        
        unique_perts = np.unique(perturbations)
        for i, pert in enumerate(unique_perts):
            mask = perturbations == pert
            if np.sum(mask) > 0:
                fig.add_trace(go.Scatter(
                    x=latent_2d[mask, 0],
                    y=latent_2d[mask, 1],
                    mode='markers',
                    marker=dict(
                        size=5,
                        color=px.colors.qualitative.Plotly[i % len(px.colors.qualitative.Plotly)],
                        opacity=0.7
                    ),
                    name=pert,
                    text=[hover_text[i] for i in np.where(mask)[0]],
                    hoverinfo='text',
                    showlegend=False
                ), row=2, col=2)
    
    # Update layout
    fig.update_layout(
        title=f"Optimization Dynamics - Epoch {epoch}",
        height=1000,
        width=1400,
        template='plotly_white',
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )
    
    # Update axis labels
    fig.update_scenes(
        xaxis_title='PC1',
        yaxis_title='PC2',
        zaxis_title='PC3',
        row=1, col=1
    )
    
    fig.update_scenes(
        xaxis_title='UMAP1',
        yaxis_title='UMAP2',
        zaxis_title='UMAP3',
        row=1, col=2
    )
    
    fig.update_xaxes(title_text="Epoch", row=2, col=1)
    fig.update_yaxes(title_text="Loss", row=2, col=1)
    fig.update_xaxes(title_text="UMAP1", row=2, col=2)
    fig.update_yaxes(title_text="UMAP2", row=2, col=2)
    
    # Save as interactive HTML
    fig.write_html(os.path.join(OUTPUT_DIR, f"optimization_epoch_{epoch}.html"))

# ==================== ENHANCED PERTURBATION ANALYSIS ====================
def analyze_perturbations(model, adata, X):
    """Analyze perturbation effects in latent space with enhanced interactive visualization"""
    print("Analyzing perturbation effects...")
    
    # Get latent representations
    with torch.no_grad():
        X_tensor = torch.FloatTensor(X).to(device)
        z, _ = model(X_tensor)
        z = z.cpu().numpy()
    
    # Calculate control mean
    control_mask = adata.obs['target_gene'] == 'non-targeting'
    control_mean = z[control_mask].mean(axis=0)
    
    # Calculate delta latent
    delta_z = z - control_mean
    
    # Prepare hover text
    hover_text = []
    for i in range(len(adata)):
        text_parts = []
        if 'target_gene' in adata.obs:
            text_parts.append(f"Perturbation: {adata.obs['target_gene'].values[i]}")
        if 'cell_type' in adata.obs and pd.notna(adata.obs['cell_type'].values[i]):
            text_parts.append(f"Cell type: {adata.obs['cell_type'].values[i]}")
        if 'n_genes' in adata.obs:
            text_parts.append(f"Genes: {adata.obs['n_genes'].values[i]}")
        hover_text.append("<br>".join(text_parts))
    
    # Create interactive visualization
    fig = make_subplots(
        rows=1, cols=2,
        specs=[[{'type': 'scatter3d'}, {'type': 'scatter'}]],
        subplot_titles=('3D UMAP Projection', '2D UMAP Projection')
    )
    
    # 3D UMAP
    reducer_3d = UMAP(n_components=3, random_state=42)
    z_umap_3d = reducer_3d.fit_transform(z)
    
    # 2D UMAP
    reducer_2d = UMAP(n_components=2, random_state=42)
    z_umap_2d = reducer_2d.fit_transform(z)
    
    # Color by perturbation
    perturbations = adata.obs['target_gene'].values
    unique_perts = np.unique(perturbations)
    
    for i, pert in enumerate(unique_perts):
        mask = perturbations == pert
        if np.sum(mask) > 0:
            # 3D plot
            fig.add_trace(go.Scatter3d(
                x=z_umap_3d[mask, 0],
                y=z_umap_3d[mask, 1],
                z=z_umap_3d[mask, 2],
                mode='markers',
                marker=dict(
                    size=4,
                    color=px.colors.qualitative.Plotly[i % len(px.colors.qualitative.Plotly)],
                    opacity=0.7
                ),
                name=pert,
                text=[hover_text[i] for i in np.where(mask)[0]],
                hoverinfo='text',
                showlegend=True
            ), row=1, col=1)
            
            # 2D plot
            fig.add_trace(go.Scatter(
                x=z_umap_2d[mask, 0],
                y=z_umap_2d[mask, 1],
                mode='markers',
                marker=dict(
                    size=5,
                    color=px.colors.qualitative.Plotly[i % len(px.colors.qualitative.Plotly)],
                    opacity=0.7
                ),
                name=pert,
                text=[hover_text[i] for i in np.where(mask)[0]],
                hoverinfo='text',
                showlegend=False
            ), row=1, col=2)
    
    fig.update_layout(
        title="Perturbation Effects in Latent Space",
        height=600,
        width=1400,
        template='plotly_white',
        legend=dict(orientation="h", yanchor="bottom", y=1.1)
    )
    
    # Update axis labels
    fig.update_scenes(
        xaxis_title='UMAP1',
        yaxis_title='UMAP2',
        zaxis_title='UMAP3',
        row=1, col=1
    )
    
    fig.update_xaxes(title_text="UMAP1", row=1, col=2)
    fig.update_yaxes(title_text="UMAP2", row=1, col=2)
    
    fig.write_html(os.path.join(OUTPUT_DIR, "perturbation_effects.html"))

# ==================== MAIN EXECUTION ====================
if __name__ == "__main__":
    print("Starting robust autoencoder analysis...")
    
    # Load data
    adata, X, gene_names = load_and_preprocess()
    print(f"Data shape: {X.shape}")
    
    # Initialize model
    model = RobustAutoencoder(input_dim=X.shape[1], latent_dim=32).to(device)
    
    # Train with trajectory tracking and visualization
    model = train_robust(model, X, adata, epochs=100, batch_size=256)
    
    # Final perturbation analysis
    analyze_perturbations(model, adata, X)
    
    print("Analysis complete! Interactive visualizations saved to:", OUTPUT_DIR)