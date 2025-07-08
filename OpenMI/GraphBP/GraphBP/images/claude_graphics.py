import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Draw
import warnings
warnings.filterwarnings('ignore')

# Set style for better-looking plots
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

# =====================================
# 1. PIPELINE OVERVIEW FLOWCHART
# =====================================
def create_pipeline_flowchart():
    """Create a flowchart showing the pipeline steps"""
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Define boxes and their positions
    boxes = [
        {"text": "Graph Generative\nModel", "xy": (1, 6), "color": "#3498db"},
        {"text": "Generated\nLigands", "xy": (3, 6), "color": "#e74c3c"},
        {"text": "Synthesizability\nScoring", "xy": (5, 6), "color": "#f39c12"},
        {"text": "QED\nCalculation", "xy": (7, 6), "color": "#27ae60"},
        {"text": "Tanimoto\nSimilarity", "xy": (9, 6), "color": "#9b59b6"},
        {"text": "Filtering\n& Selection", "xy": (11, 6), "color": "#e67e22"},
        {"text": "AutoDock Vina\nDocking", "xy": (6, 3), "color": "#2c3e50"},
        {"text": "Best Compounds\nfor Aurora B", "xy": (6, 1), "color": "#c0392b"}
    ]
    
    # Draw boxes
    for box in boxes:
        rect = Rectangle((box["xy"][0]-0.7, box["xy"][1]-0.4), 1.4, 0.8, 
                        facecolor=box["color"], alpha=0.7, edgecolor='black')
        ax.add_patch(rect)
        ax.text(box["xy"][0], box["xy"][1], box["text"], ha='center', va='center', 
                fontsize=10, weight='bold', color='white')
    
    # Draw arrows
    arrow_props = dict(arrowstyle='->', lw=2, color='black')
    arrows = [
        ((1.7, 6), (2.3, 6)),  # Model to Ligands
        ((3.7, 6), (4.3, 6)),  # Ligands to Synth
        ((5.7, 6), (6.3, 6)),  # Synth to QED
        ((7.7, 6), (8.3, 6)),  # QED to Tanimoto
        ((9.7, 6), (10.3, 6)), # Tanimoto to Filtering
        ((11, 5.6), (6.7, 3.4)), # Filtering to Docking
        ((6, 2.6), (6, 1.4))   # Docking to Results
    ]
    
    for start, end in arrows:
        ax.annotate('', xy=end, xytext=start, arrowprops=arrow_props)
    
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Aurora Kinase B Ligand Discovery Pipeline', fontsize=16, weight='bold', pad=20)
    
    plt.tight_layout()
    return fig

# =====================================
# 2. MOLECULAR PROPERTY DISTRIBUTIONS
# =====================================
def plot_molecular_properties(df):
    """Plot distributions of molecular properties"""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    properties = ['synthesizability', 'qed', 'tanimoto_similarity', 
                 'molecular_weight', 'logp', 'num_rotatable_bonds']
    # properties = ["SA_score", "SCScore", "NP_score", "QED", 
    #               "tanimoto", "len_smiles"]
    
    for i, prop in enumerate(properties):
        row, col = i // 3, i % 3
        ax = axes[row, col]
        
        if prop in df.columns:
            # Histogram with KDE
            ax.hist(df[prop], bins=30, alpha=0.7, density=True, color='skyblue', edgecolor='black')
            
            # Add KDE curve
            from scipy import stats
            x = np.linspace(df[prop].min(), df[prop].max(), 100)
            kde = stats.gaussian_kde(df[prop])
            ax.plot(x, kde(x), 'r-', linewidth=2, label='KDE')
            
            # Add vertical line for mean
            mean_val = df[prop].mean()
            ax.axvline(mean_val, color='red', linestyle='--', alpha=0.8, 
                      label=f'Mean: {mean_val:.3f}')
            
            ax.set_title(f'{prop.replace("_", " ").title()} Distribution', fontsize=12, weight='bold')
            ax.set_xlabel(prop.replace("_", " ").title())
            ax.set_ylabel('Density')
            ax.legend()
            ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig

# =====================================
# 3. PROPERTY CORRELATION HEATMAP
# =====================================
def create_correlation_heatmap(df):
    """Create correlation heatmap of molecular properties"""
    # Select numerical columns
    numerical_cols = df.select_dtypes(include=[np.number]).columns
    correlation_matrix = df[numerical_cols].corr()
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create heatmap
    sns.heatmap(correlation_matrix, annot=True, cmap='RdYlBu_r', center=0,
                square=True, fmt='.2f', cbar_kws={"shrink": .8})
    
    ax.set_title('Molecular Properties Correlation Matrix', fontsize=14, weight='bold', pad=20)
    plt.tight_layout()
    return fig

# =====================================
# 4. FILTERING FUNNEL VISUALIZATION
# =====================================
def create_filtering_funnel(stages_data):
    """Create a funnel chart showing compound filtering stages"""
    fig = go.Figure()
    
    # Example stages - replace with your actual data
    stages = list(stages_data.keys())
    counts = list(stages_data.values())
    
    fig.add_trace(go.Funnel(
        y=stages,
        x=counts,
        textinfo="value+percent initial",
        textposition="inside",
        opacity=0.8,
        marker=dict(
            color=["#3498db", "#e74c3c", "#f39c12", "#27ae60", "#9b59b6"],
            line=dict(width=2, color="white")
        )
    ))
    
    fig.update_layout(
        title="Compound Filtering Funnel",
        title_x=0.5,
        font=dict(size=12),
        height=500,
        width=700
    )
    
    return fig

# =====================================
# 5. DOCKING SCORE ANALYSIS
# =====================================
def plot_docking_analysis(df):
    """Create comprehensive docking score analysis"""
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    # 1. Docking score distribution
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(df['docking_score'], bins=30, alpha=0.7, color='lightcoral', edgecolor='black')
    ax1.axvline(df['docking_score'].mean(), color='red', linestyle='--', 
               label=f'Mean: {df["docking_score"].mean():.2f}')
    ax1.set_title('Docking Score Distribution')
    ax1.set_xlabel('AutoDock Vina Score (kcal/mol)')
    ax1.set_ylabel('Frequency')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Docking score vs QED
    ax2 = fig.add_subplot(gs[0, 1])
    scatter = ax2.scatter(df['qed'], df['docking_score'], 
                         c=df['synthesizability'], cmap='viridis', alpha=0.7)
    ax2.set_xlabel('QED Score')
    ax2.set_ylabel('Docking Score (kcal/mol)')
    ax2.set_title('Docking Score vs QED\n(colored by Synthesizability)')
    plt.colorbar(scatter, ax=ax2, label='Synthesizability')
    
    # 3. Docking score vs Tanimoto similarity
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.scatter(df['tanimoto_similarity'], df['docking_score'], alpha=0.7, color='orange')
    ax3.set_xlabel('Tanimoto Similarity')
    ax3.set_ylabel('Docking Score (kcal/mol)')
    ax3.set_title('Docking Score vs Tanimoto Similarity')
    
    # 4. Top compounds bar plot
    ax4 = fig.add_subplot(gs[1, :])
    top_compounds = df.nsmallest(10, 'docking_score')
    bars = ax4.bar(range(len(top_compounds)), top_compounds['docking_score'], 
                   color='lightgreen', alpha=0.8, edgecolor='black')
    ax4.set_xlabel('Compound Rank')
    ax4.set_ylabel('Docking Score (kcal/mol)')
    ax4.set_title('Top 10 Compounds by Docking Score')
    ax4.set_xticks(range(len(top_compounds)))
    ax4.set_xticklabels([f'Compound {i+1}' for i in range(len(top_compounds))], rotation=45)
    
    # Add value labels on bars
    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height - 0.1,
                f'{height:.2f}', ha='center', va='top', fontweight='bold')
    
    # 5. Property radar chart for top compound
    ax5 = fig.add_subplot(gs[2, 0], projection='polar')
    top_compound = df.loc[df['docking_score'].idxmin()]
    
    properties = ['qed', 'synthesizability', 'tanimoto_similarity']
    values = [top_compound[prop] for prop in properties]
    
    angles = np.linspace(0, 2*np.pi, len(properties), endpoint=False).tolist()
    values += values[:1]  # Complete the circle
    angles += angles[:1]
    
    ax5.plot(angles, values, 'o-', linewidth=2, color='red')
    ax5.fill(angles, values, alpha=0.25, color='red')
    ax5.set_xticks(angles[:-1])
    ax5.set_xticklabels(properties)
    ax5.set_title('Top Compound Properties')
    ax5.set_ylim(0, 1)
    
    # 6. Binding affinity ranges
    ax6 = fig.add_subplot(gs[2, 1:])
    
    # Create binding affinity categories
    def categorize_binding(score):
        if score < -9:
            return 'Very Strong'
        elif score < -7:
            return 'Strong'
        elif score < -5:
            return 'Moderate'
        else:
            return 'Weak'
    
    df['binding_category'] = df['docking_score'].apply(categorize_binding)
    binding_counts = df['binding_category'].value_counts()
    
    wedges, texts, autotexts = ax6.pie(binding_counts.values, labels=binding_counts.index, 
                                      autopct='%1.1f%%', startangle=90,
                                      colors=['#2ecc71', '#f39c12', '#e74c3c', '#95a5a6'])
    ax6.set_title('Binding Affinity Categories')
    
    plt.suptitle('Comprehensive Docking Analysis', fontsize=16, weight='bold')
    return fig

# =====================================
# 6. MULTI-OBJECTIVE OPTIMIZATION PLOT
# =====================================
def plot_pareto_front(df):
    """Create Pareto front analysis for multi-objective optimization"""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # 1. QED vs Docking Score
    ax1 = axes[0]
    scatter1 = ax1.scatter(df['qed'], df['docking_score'], 
                          c=df['synthesizability'], cmap='plasma', alpha=0.7, s=50)
    ax1.set_xlabel('QED Score')
    ax1.set_ylabel('Docking Score (kcal/mol)')
    ax1.set_title('QED vs Docking Score')
    ax1.grid(True, alpha=0.3)
    plt.colorbar(scatter1, ax=ax1, label='Synthesizability')
    
    # 2. Synthesizability vs Docking Score
    ax2 = axes[1]
    scatter2 = ax2.scatter(df['synthesizability'], df['docking_score'], 
                          c=df['qed'], cmap='coolwarm', alpha=0.7, s=50)
    ax2.set_xlabel('Synthesizability Score')
    ax2.set_ylabel('Docking Score (kcal/mol)')
    ax2.set_title('Synthesizability vs Docking Score')
    ax2.grid(True, alpha=0.3)
    plt.colorbar(scatter2, ax=ax2, label='QED')
    
    # 3. 3D scatter plot
    ax3 = axes[2]
    ax3.remove()
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    
    scatter3 = ax3.scatter(df['qed'], df['synthesizability'], df['docking_score'],
                          c=df['tanimoto_similarity'], cmap='viridis', alpha=0.7, s=50)
    ax3.set_xlabel('QED Score')
    ax3.set_ylabel('Synthesizability')
    ax3.set_zlabel('Docking Score')
    ax3.set_title('3D Multi-Objective Space')
    plt.colorbar(scatter3, ax=ax3, label='Tanimoto Similarity', shrink=0.5)
    
    plt.tight_layout()
    return fig

# =====================================
# 7. INTERACTIVE PLOTLY DASHBOARD
# =====================================
def create_interactive_dashboard(df):
    """Create an interactive dashboard with multiple linked plots"""
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Property Distributions', 'Docking Score Analysis', 
                       'Correlation Network', 'Top Compounds'),
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}]]
    )
    
    # 1. Property distributions
    for i, prop in enumerate(['qed', 'synthesizability', 'tanimoto_similarity']):
        fig.add_trace(
            go.Histogram(x=df[prop], name=prop.title(), opacity=0.7),
            row=1, col=1
        )
    
    # 2. Docking score scatter
    fig.add_trace(
        go.Scatter(x=df['qed'], y=df['docking_score'], 
                  mode='markers', name='Compounds',
                  marker=dict(size=8, color=df['synthesizability'], 
                            colorscale='Viridis', showscale=True)),
        row=1, col=2
    )
    
    # 3. Top compounds
    top_10 = df.nsmallest(10, 'docking_score')
    fig.add_trace(
        go.Bar(x=list(range(len(top_10))), y=top_10['docking_score'],
               name='Top Compounds', marker_color='lightgreen'),
        row=2, col=1
    )
    
    # 4. Property correlation
    corr_matrix = df[['qed', 'synthesizability', 'tanimoto_similarity', 'docking_score']].corr()
    fig.add_trace(
        go.Heatmap(z=corr_matrix.values, x=corr_matrix.columns, y=corr_matrix.columns,
                  colorscale='RdYlBu', zmid=0),
        row=2, col=2
    )
    
    fig.update_layout(height=800, showlegend=True, 
                     title_text="Aurora Kinase B Ligand Analysis Dashboard")
    
    return fig

# =====================================
# EXAMPLE USAGE AND DATA GENERATION
# =====================================
def generate_example_data(n_compounds=1000):
    """Generate example data for demonstration"""
    np.random.seed(42)
    
    data = {
        'compound_id': [f'COMP_{i:04d}' for i in range(n_compounds)],
        'qed': np.random.beta(2, 2, n_compounds),
        'synthesizability': np.random.beta(3, 2, n_compounds),
        'tanimoto_similarity': np.random.beta(1.5, 3, n_compounds),
        'docking_score': np.random.normal(-6, 2, n_compounds),
        'molecular_weight': np.random.normal(350, 100, n_compounds),
        'logp': np.random.normal(2.5, 1.5, n_compounds),
        'num_rotatable_bonds': np.random.poisson(5, n_compounds)
    }
    
    return pd.DataFrame(data)

# Example usage
if __name__ == "__main__":
    # Generate example data
    df = generate_example_data(1000)
    # df = synth_df
    
    # Create all visualizations
    print("Creating pipeline flowchart...")
    fig1 = create_pipeline_flowchart()
    fig1.savefig('pipeline_flowchart.png', dpi=300, bbox_inches='tight')
    
    print("Creating molecular property distributions...")
    fig2 = plot_molecular_properties(df)
    fig2.savefig('molecular_properties.png', dpi=300, bbox_inches='tight')
    
    print("Creating correlation heatmap...")
    fig3 = create_correlation_heatmap(df)
    fig3.savefig('correlation_heatmap.png', dpi=300, bbox_inches='tight')
    
    print("Creating filtering funnel...")
    stages_data = {
        'Initial Generated': 10000,
        'After Synthesizability Filter': 5000,
        'After QED Filter': 2000,
        'After Tanimoto Filter': 1000,
        'Final Docked': 500
    }
    fig4 = create_filtering_funnel(stages_data)
    fig4.write_html('filtering_funnel.html')
    
    print("Creating docking analysis...")
    fig5 = plot_docking_analysis(df)
    fig5.savefig('docking_analysis.png', dpi=300, bbox_inches='tight')
    
    print("Creating Pareto front analysis...")
    fig6 = plot_pareto_front(df)
    fig6.savefig('pareto_analysis.png', dpi=300, bbox_inches='tight')
    
    print("Creating interactive dashboard...")
    fig7 = create_interactive_dashboard(df)
    fig7.write_html('interactive_dashboard.html')
    
    print("All visualizations created successfully!")
    
    plt.show()