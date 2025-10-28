import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import warnings
import io
from datetime import datetime
import time
import requests
import logging
import re
import signal
from contextlib import contextmanager
from typing import List, Tuple, Dict
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET
from Bio import Entrez
import plotly.graph_objects as go
import plotly.express as px
import networkx as nx

warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Configuration Management
class Config:
    MAX_GENES = 200
    REQUEST_TIMEOUT = 30
    CACHE_TTL = 3600
    MAX_RETRIES = 3
    REQUESTS_PER_SECOND = 2

class AppConfig:
    def __init__(self):
        self.ncbi_email = ""
        self.debug_mode = False
        self.max_workers = 3
        self.cache_ttl = 3600

# Decorator for safe API calls
def safe_api_call(func):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except requests.exceptions.Timeout:
            logger.error(f"Timeout in {func.__name__}")
            return None
        except requests.exceptions.ConnectionError:
            logger.error(f"Connection error in {func.__name__}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error in {func.__name__}: {e}")
            return None
    return wrapper

# Set page configuration
st.set_page_config(
    page_title="Integrated Gene Analysis Pipeline",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# =============================================================================
# CONSISTENT DARK THEME CSS FOR ALL MODES
# =============================================================================

def get_dark_theme_css():
    return """
<style>
  @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap');

  :root {
    --bg: #0f1419;
    --panel: #1a1f2e;
    --surface: #252b3a;
    --text: #ffffff;
    --text-muted: #b3b8c5;
    --text-secondary: #8b92a5;
    --border: #2d3748;
    --border-strong: #4a5568;
    --accent: #00d4aa;
    --accent-hover: #00b894;
    --accent-active: #00a085;
    --secondary: #667eea;
    --secondary-hover: #5a67d8;
    --danger: #f56565;
    --warning: #ed8936;
    --success: #48bb78;
    --input-bg: #2d3748;
    --shadow: rgba(0, 0, 0, 0.3);
  }

  /* Global Styles */
  .stApp {
    background: var(--bg);
    color: var(--text);
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
  }

  /* Header */
  header[data-testid="stHeader"] {
    background: var(--bg) !important;
    border-bottom: 1px solid var(--border);
  }

  .block-container {
    padding-top: 2rem !important;
    padding-bottom: 2rem;
    max-width: 1200px;
  }

  /* Sidebar */
  section[data-testid="stSidebar"] {
    background: var(--panel) !important;
    border-right: 1px solid var(--border);
  }

  section[data-testid="stSidebar"] * {
    color: var(--text) !important;
  }

  /* Mode Selection Buttons */
  .mode-button {
    background: linear-gradient(135deg, var(--panel) 0%, var(--surface) 100%) !important;
    border: 2px solid var(--border) !important;
    border-radius: 12px !important;
    padding: 2rem 1rem !important;
    color: var(--text) !important;
    font-weight: 600 !important;
    font-size: 1.1rem !important;
    transition: all 0.3s ease !important;
    height: 120px !important;
    display: flex !important;
    align-items: center !important;
    justify-content: center !important;
    text-align: center !important;
  }

  .mode-button:hover {
    border-color: var(--accent) !important;
    transform: translateY(-2px) !important;
    box-shadow: 0 8px 25px rgba(0, 212, 170, 0.2) !important;
  }

  .mode-button:active {
    transform: translateY(0) !important;
  }

  .mode-button.active {
    border-color: var(--accent) !important;
    background: linear-gradient(135deg, var(--accent), var(--secondary)) !important;
    color: white !important;
    box-shadow: 0 4px 15px rgba(0, 212, 170, 0.3) !important;
  }

  /* Hero Section */
  .hero {
    background: linear-gradient(135deg, var(--panel) 0%, var(--surface) 100%);
    border: 1px solid var(--border);
    border-radius: 16px;
    padding: 2rem;
    margin: 1rem 0 2rem 0;
    box-shadow: 0 4px 20px var(--shadow);
    display: flex;
    align-items: center;
    gap: 2rem;
    min-height: 180px;
  }

  .hero-content {
    flex: 1;
  }

  .hero-logo {
    flex-shrink: 0;
    display: flex;
    align-items: center;
    justify-content: center;
  }

  .hero h1 {
    font-size: 2.5rem;
    font-weight: 800;
    margin: 0 0 0.5rem 0;
    background: linear-gradient(135deg, var(--accent), var(--secondary));
    -webkit-background-clip: text;
    background-clip: text;
    color: transparent;
    line-height: 1.2;
  }

  .hero p {
    color: var(--text-muted);
    font-size: 1.1rem;
    margin: 0;
    line-height: 1.5;
  }

  .logo-img {
    height: 120px !important;
    width: auto !important;
    max-width: 120px !important;
    object-fit: contain;
    border-radius: 12px;
  }

  /* Section Titles */
  .section-title {
    font-size: 1.4rem;
    font-weight: 700;
    color: var(--text);
    margin: 2rem 0 1rem 0;
    display: flex;
    align-items: center;
    gap: 0.75rem;
    padding-bottom: 0.5rem;
    border-bottom: 2px solid var(--border);
  }

  .section-title .icon {
    font-size: 1.3rem;
  }

  /* Hints */
  .hint {
    color: var(--text-secondary);
    font-size: 0.95rem;
    margin-bottom: 1.5rem;
    padding: 1rem;
    background: var(--surface);
    border-radius: 8px;
    border-left: 4px solid var(--accent);
  }

  /* Input Fields */
  .stTextInput > div > div,
  .stTextArea > div > div {
    background: var(--input-bg) !important;
    border: 2px solid var(--border) !important;
    border-radius: 12px !important;
    transition: border-color 0.2s ease;
  }

  .stTextInput > div > div:focus-within,
  .stTextArea > div > div:focus-within {
    border-color: var(--accent) !important;
    box-shadow: 0 0 0 3px rgba(0, 212, 170, 0.1) !important;
  }

  .stTextInput input,
  .stTextArea textarea {
    background: transparent !important;
    color: var(--text) !important;
    font-family: 'Inter', sans-serif;
  }

  .stTextInput input::placeholder,
  .stTextArea textarea::placeholder {
    color: var(--text-secondary) !important;
  }

  /* File Uploader */
  .stFileUploader [data-testid="stFileUploaderDropzone"] {
    background: var(--surface) !important;
    border: 2px solid var(--border) !important;
    border-radius: 12px !important;
    padding: 2rem !important;
    transition: all 0.2s ease;
  }

  .stFileUploader [data-testid="stFileUploaderDropzone"]:hover {
    border-color: var(--accent) !important;
    background: var(--input-bg) !important;
  }

  .stFileUploader [data-testid="stFileUploaderDropzone"] * {
    color: var(--text) !important;
  }

  /* Selectbox */
  .stSelectbox [data-baseweb="select"] > div {
    background: var(--input-bg) !important;
    border: 2px solid var(--border) !important;
    border-radius: 12px !important;
  }

  [data-baseweb="popover"] [role="listbox"] {
    background: var(--panel) !important;
    border: 1px solid var(--border) !important;
    border-radius: 8px !important;
  }

  /* Buttons */
  .stButton > button {
    background: linear-gradient(135deg, var(--accent), var(--secondary)) !important;
    color: white !important;
    border: none !important;
    border-radius: 12px !important;
    padding: 0.75rem 2rem !important;
    font-weight: 600 !important;
    font-size: 1rem !important;
    transition: all 0.2s ease !important;
    box-shadow: 0 4px 15px rgba(0, 212, 170, 0.3) !important;
  }

  .stButton > button:hover {
    transform: translateY(-2px) !important;
    box-shadow: 0 8px 25px rgba(0, 212, 170, 0.4) !important;
    background: linear-gradient(135deg, var(--accent-hover), var(--secondary-hover)) !important;
  }

  .stButton > button:active {
    transform: translateY(0) !important;
  }

  /* Checkboxes */
  .stCheckbox > label {
    color: var(--text) !important;
    font-weight: 500;
  }

  /* Tabs */
  .stTabs [data-baseweb="tab-list"] {
    gap: 0.5rem;
  }

  .stTabs [data-baseweb="tab"] {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 8px 8px 0 0;
    color: var(--text-muted);
    font-weight: 500;
    padding: 0.75rem 1.5rem;
    transition: all 0.2s ease;
  }

  .stTabs [aria-selected="true"] {
    background: var(--panel);
    color: var(--accent);
    border-bottom-color: var(--panel);
  }

  /* DataFrames */
  .stDataFrame {
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 12px;
    overflow: hidden;
  }

  .stDataFrame table {
    background: var(--panel) !important;
  }

  .stDataFrame th {
    background: var(--surface) !important;
    color: var(--text) !important;
    font-weight: 600;
    border-bottom: 2px solid var(--border) !important;
  }

  .stDataFrame td {
    background: var(--panel) !important;
    color: var(--text) !important;
    border-bottom: 1px solid var(--border) !important;
  }

  /* Progress Bar */
  .stProgress > div > div > div {
    background: linear-gradient(90deg, var(--accent), var(--secondary)) !important;
  }

  /* Alerts */
  .stAlert {
    background: var(--surface) !important;
    border: 1px solid var(--border) !important;
    border-radius: 8px !important;
    color: var(--text) !important;
  }

  .stSuccess {
    border-left: 4px solid var(--success) !important;
  }

  .stError {
    border-left: 4px solid var(--danger) !important;
  }

  .stWarning {
    border-left: 4px solid var(--warning) !important;
  }

  .stInfo {
    border-left: 4px solid var(--secondary) !important;
  }

  /* Sidebar Content */
  .sidebar-title {
    font-weight: 700;
    font-size: 1.2rem;
    margin-bottom: 0.5rem;
    color: var(--accent);
  }

  .sidebar-tip {
    color: var(--text-muted);
    font-size: 0.9rem;
    line-height: 1.4;
  }

  /* Download Buttons */
  .stDownloadButton > button {
    background: var(--surface) !important;
    color: var(--text) !important;
    border: 1px solid var(--border) !important;
    border-radius: 8px !important;
    font-weight: 500 !important;
    transition: all 0.2s ease !important;
  }

  .stDownloadButton > button:hover {
    background: var(--border) !important;
    border-color: var(--accent) !important;
  }

  /* Plotly Charts */
  .js-plotly-plot {
    background: var(--panel) !important;
    border-radius: 12px;
    overflow: hidden;
  }

  /* Responsive Design */
  @media (max-width: 768px) {
    .hero {
      flex-direction: column;
      text-align: center;
      gap: 1rem;
      padding: 1.5rem;
    }
    
    .hero h1 {
      font-size: 2rem;
    }
    
    .logo-img {
      height: 80px !important;
      max-width: 80px !important;
    }
    
    .block-container {
      padding-left: 1rem;
      padding-right: 1rem;
    }
    
    .section-title {
      font-size: 1.2rem;
    }

    .mode-button {
      padding: 1.5rem 0.5rem !important;
      height: 100px !important;
      font-size: 1rem !important;
    }
  }

  /* Scrollbars */
  ::-webkit-scrollbar {
    width: 8px;
    height: 8px;
  }

  ::-webkit-scrollbar-track {
    background: var(--surface);
  }

  ::-webkit-scrollbar-thumb {
    background: var(--border-strong);
    border-radius: 4px;
  }

  ::-webkit-scrollbar-thumb:hover {
    background: var(--accent);
  }
</style>
"""

# Apply dark theme to entire app
st.markdown(get_dark_theme_css(), unsafe_allow_html=True)

# =============================================================================
# SESSION STATE INITIALIZATION
# =============================================================================

def initialize_session_state():
    """Initialize all session state variables with proper defaults"""
    if 'analysis_mode' not in st.session_state:
        st.session_state.analysis_mode = None
    if 'show_pathway_analysis' not in st.session_state:
        st.session_state.show_pathway_analysis = False
    if 'pathway_genes' not in st.session_state:
        st.session_state.pathway_genes = None
    if 'degs_completed' not in st.session_state:
        st.session_state.degs_completed = False
    if 'current_pipeline_step' not in st.session_state:
        st.session_state.current_pipeline_step = 'selection'  # 'selection', 'degs', 'pathway'
    if 'degs_results' not in st.session_state:
        st.session_state.degs_results = None
    if 'analysis_results' not in st.session_state:
        st.session_state.analysis_results = {}

# =============================================================================
# DEGs ANALYSIS FUNCTIONS (First Pipeline)
# =============================================================================

def calculate_differential_expression_fast(count_matrix, sample_group1, sample_group2):
    """
    Ultra-fast version using batch processing for differential expression analysis
    """
    # Select data as numpy arrays for faster computation
    data_group1 = count_matrix[sample_group1].values
    data_group2 = count_matrix[sample_group2].values
    
    # Calculate means using numpy (much faster than pandas)
    mean_group1 = np.mean(data_group1, axis=1)
    mean_group2 = np.mean(data_group2, axis=1)
    
    # Calculate logFC with pseudocount
    pseudocount = 0.1
    logFC = np.log2((mean_group2 + pseudocount) / (mean_group1 + pseudocount))
    
    # Batch t-test calculation for better performance
    p_values = []
    batch_size = 1000  # Process 1000 genes at a time to balance memory and speed
    
    total_genes = len(count_matrix)
    
    for i in range(0, total_genes, batch_size):
        end_idx = min(i + batch_size, total_genes)
        
        batch_group1 = data_group1[i:end_idx]
        batch_group2 = data_group2[i:end_idx]
        
        for j in range(len(batch_group1)):
            g1_data = batch_group1[j]
            g2_data = batch_group2[j]
            
            # Fast data validation checks
            if (np.isnan(g1_data).any() or np.isnan(g2_data).any() or 
                np.isinf(g1_data).any() or np.isinf(g2_data).any()):
                p_values.append(1.0)
                continue
            
            # Check for zero variance cases
            if (np.std(g1_data) == 0 and np.std(g2_data) == 0 and 
                np.mean(g1_data) == np.mean(g2_data)):
                p_values.append(1.0)
                continue
            
            try:
                t_stat, p_val = ttest_ind(g1_data, g2_data, equal_var=False, nan_policy='omit')
                p_values.append(p_val if not np.isnan(p_val) else 1.0)
            except:
                p_values.append(1.0)
    
    # Create results DataFrame
    results = pd.DataFrame({
        'Gene': count_matrix.index,
        'logFC': logFC,
        'p_value': p_values,
        'mean_group1': mean_group1,
        'mean_group2': mean_group2
    })
    
    return results

def filter_and_sort_degs(results, logFC_threshold=2, p_value_threshold=0.05):
    """
    Filter DEGs based on logFC and p-value thresholds
    """
    # Filter significant genes using vectorized operations
    significant_mask = (np.abs(results['logFC']) > logFC_threshold) & (results['p_value'] < p_value_threshold)
    significant_genes = results[significant_mask].copy()
    
    # Add direction column
    significant_genes['direction'] = np.where(
        significant_genes['logFC'] > 0, 'upregulated', 'downregulated'
    )
    
    # Separate upregulated and downregulated genes
    upregulated = significant_genes[significant_genes['logFC'] > 0].sort_values('logFC', ascending=False)
    downregulated = significant_genes[significant_genes['logFC'] < 0].sort_values('logFC', ascending=True)
    
    return upregulated, downregulated

def check_data_quality(count_matrix, group1_samples, group2_samples):
    """
    Check data quality and report any issues
    """
    quality_issues = []
    
    # Check for NaN values
    nan_count = count_matrix[group1_samples + group2_samples].isna().sum().sum()
    if nan_count > 0:
        quality_issues.append(f"Found {nan_count} NaN values in the data")
    
    # Check for infinite values
    inf_count = np.isinf(count_matrix[group1_samples + group2_samples].values).sum()
    if inf_count > 0:
        quality_issues.append(f"Found {inf_count} infinite values in the data")
    
    # Check for zero variance genes using vectorized operations
    all_samples = group1_samples + group2_samples
    variances = count_matrix[all_samples].var(axis=1)
    zero_var_genes = (variances == 0).sum()
    if zero_var_genes > 0:
        quality_issues.append(f"Found {zero_var_genes} genes with zero variance across all samples")
    
    # Check for negative values
    negative_count = (count_matrix[group1_samples + group2_samples] < 0).sum().sum()
    if negative_count > 0:
        quality_issues.append(f"Found {negative_count} negative values in the data")
    
    return quality_issues

@st.cache_data
def cached_calculate_de(_count_matrix, sample_group1, sample_group2, logFC_threshold, p_value_threshold):
    """
    Cached version of differential expression analysis
    """
    results = calculate_differential_expression_fast(_count_matrix, sample_group1, sample_group2)
    upregulated, downregulated = filter_and_sort_degs(results, logFC_threshold, p_value_threshold)
    return results, upregulated, downregulated

def run_degs_analysis():
    """Run the DEGs analysis pipeline"""
    st.title("🧬 Differential Expression Analyzer")
    st.markdown("""
    This tool performs differential expression analysis between two sample groups using RNA-seq count data.
    Upload your count matrix and configure the analysis parameters below.
    
    **⚡ Now with optimized performance for faster analysis!**
    """)
    
    # Sidebar for file upload and parameters
    st.sidebar.header("Upload Data")
    
    uploaded_file = st.sidebar.file_uploader(
        "Choose a count matrix file",
        type=['csv', 'tsv', 'txt'],
        help="Upload a CSV or TSV file with genes as rows and samples as columns"
    )
    
    st.sidebar.header("Analysis Parameters")
    logFC_threshold = st.sidebar.number_input(
        "logFC Threshold",
        min_value=0.0,
        max_value=10.0,
        value=2.0,
        step=0.1,
        help="Minimum absolute log2 fold change for significance"
    )
    
    p_value_threshold = st.sidebar.number_input(
        "P-value Threshold",
        min_value=0.0,
        max_value=1.0,
        value=0.05,
        step=0.01,
        help="Maximum p-value for significance"
    )
    
    # Performance settings
    st.sidebar.header("Performance Settings")
    use_caching = st.sidebar.checkbox(
        "Enable result caching", 
        value=True,
        help="Cache results for faster re-analysis with same parameters"
    )
    
    if uploaded_file is not None:
        try:
            # Determine file format and read data
            if uploaded_file.name.endswith('.csv'):
                count_matrix = pd.read_csv(uploaded_file, index_col=0)
            else:
                count_matrix = pd.read_csv(uploaded_file, sep='\t', index_col=0)
            
            st.success(f"Data loaded successfully: {count_matrix.shape[0]} genes, {count_matrix.shape[1]} samples")
            
            # Display data preview
            with st.expander("Data Preview"):
                st.dataframe(count_matrix.head(), use_container_width=True)
                st.write(f"**Shape:** {count_matrix.shape[0]} rows × {count_matrix.shape[1]} columns")
            
            # Sample selection
            st.header("Sample Group Configuration")
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Group 1")
                group1_name = st.text_input("Group 1 name", value="Control", key="group1_name")
                group1_samples = st.multiselect(
                    f"Select samples for {group1_name}",
                    options=count_matrix.columns.tolist(),
                    key="group1_samples"
                )
                st.write(f"Selected {len(group1_samples)} samples")
            
            with col2:
                st.subheader("Group 2")
                group2_name = st.text_input("Group 2 name", value="Treatment", key="group2_name")
                group2_samples = st.multiselect(
                    f"Select samples for {group2_name}",
                    options=[col for col in count_matrix.columns if col not in group1_samples],
                    key="group2_samples"
                )
                st.write(f"Selected {len(group2_samples)} samples")
            
            # Check for sample overlap
            overlap = set(group1_samples) & set(group2_samples)
            if overlap:
                st.error(f"Sample overlap detected: {list(overlap)}. Please ensure groups have distinct samples.")
                return
            
            if not group1_samples or not group2_samples:
                st.warning("Please select samples for both groups to proceed with analysis.")
                return
            
            # Data quality check
            st.header("Data Quality Check")
            quality_issues = check_data_quality(count_matrix, group1_samples, group2_samples)
            
            if quality_issues:
                for issue in quality_issues:
                    st.warning(issue)
            else:
                st.success("No major data quality issues detected.")
            
            # Performance estimation
            total_genes = count_matrix.shape[0]
            estimated_time = max(5, total_genes // 500)  # Rough estimate: ~500 genes per second
            st.info(f"⏱️ Estimated analysis time: {estimated_time} seconds for {total_genes:,} genes")
            
            # Run analysis
            if st.button("Run Differential Expression Analysis", type="primary"):
                # Setup progress tracking
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                # Analysis steps with progress updates
                status_text.text("Initializing analysis...")
                progress_bar.progress(10)
                time.sleep(0.5)
                
                status_text.text("Calculating fold changes...")
                progress_bar.progress(30)
                time.sleep(0.5)
                
                try:
                    # Perform differential expression analysis
                    if use_caching:
                        results, upregulated, downregulated = cached_calculate_de(
                            count_matrix, group1_samples, group2_samples, logFC_threshold, p_value_threshold
                        )
                    else:
                        results = calculate_differential_expression_fast(count_matrix, group1_samples, group2_samples)
                        status_text.text("Filtering significant genes...")
                        progress_bar.progress(70)
                        upregulated, downregulated = filter_and_sort_degs(results, logFC_threshold, p_value_threshold)
                    
                    progress_bar.progress(100)
                    status_text.text("Analysis complete!")
                    time.sleep(0.5)
                    
                    # Clear progress indicators
                    progress_bar.empty()
                    status_text.empty()
                    
                    # Store results in session state
                    st.session_state.degs_results = {
                        'results': results,
                        'upregulated': upregulated,
                        'downregulated': downregulated,
                        'group1_name': group1_name,
                        'group2_name': group2_name
                    }
                    
                    # Mark DEGs as completed
                    st.session_state.degs_completed = True
                    
                    # Display results
                    st.header("Analysis Results")
                    
                    # Summary statistics
                    col1, col2, col3, col4 = st.columns(4)
                    
                    with col1:
                        st.metric("Total Genes", len(results))
                    
                    with col2:
                        st.metric("Significant DEGs", len(upregulated) + len(downregulated))
                    
                    with col3:
                        st.metric("Upregulated", len(upregulated))
                    
                    with col4:
                        st.metric("Downregulated", len(downregulated))
                    
                    # Display top genes
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.subheader(f"Top Upregulated Genes in {group2_name}")
                        if not upregulated.empty:
                            display_up = upregulated.head(10)[['Gene', 'logFC', 'p_value', 'mean_group1', 'mean_group2']].round(4)
                            st.dataframe(display_up, use_container_width=True)
                        else:
                            st.info("No upregulated genes found matching the criteria.")
                    
                    with col2:
                        st.subheader(f"Top Downregulated Genes in {group2_name}")
                        if not downregulated.empty:
                            display_down = downregulated.head(10)[['Gene', 'logFC', 'p_value', 'mean_group1', 'mean_group2']].round(4)
                            st.dataframe(display_down, use_container_width=True)
                        else:
                            st.info("No downregulated genes found matching the criteria.")
                    
                    # Volcano plot
                    st.subheader("Volcano Plot")
                    try:
                        import plotly.express as px
                        
                        # Prepare data for volcano plot
                        plot_data = results.copy()
                        plot_data['-log10(p_value)'] = -np.log10(plot_data['p_value'])
                        plot_data['Significant'] = (np.abs(plot_data['logFC']) > logFC_threshold) & (plot_data['p_value'] < p_value_threshold)
                        
                        fig = px.scatter(
                            plot_data,
                            x='logFC',
                            y='-log10(p_value)',
                            color='Significant',
                            hover_data=['Gene'],
                            title=f"Volcano Plot: {group1_name} vs {group2_name}",
                            labels={'logFC': 'log2 Fold Change', '-log10(p_value)': '-log10(p-value)'},
                            color_discrete_map={True: 'red', False: 'gray'}
                        )
                        
                        # Add threshold lines
                        fig.add_vline(x=logFC_threshold, line_dash="dash", line_color="red")
                        fig.add_vline(x=-logFC_threshold, line_dash="dash", line_color="red")
                        fig.add_hline(y=-np.log10(p_value_threshold), line_dash="dash", line_color="red")
                        
                        st.plotly_chart(fig, use_container_width=True)
                    except ImportError:
                        st.warning("Install plotly to generate volcano plots: pip install plotly")
                    
                    # Download results
                    st.header("Download Results")
                    
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    
                    # All significant DEGs
                    all_degs = pd.concat([upregulated, downregulated])
                    if not all_degs.empty:
                        csv_all = all_degs.to_csv(index=False)
                        st.download_button(
                            label="Download All Significant DEGs",
                            data=csv_all,
                            file_name=f"all_significant_DEGs_{timestamp}.csv",
                            mime="text/csv"
                        )
                    
                    # Complete results
                    csv_complete = results.to_csv(index=False)
                    st.download_button(
                        label="Download Complete Results",
                        data=csv_complete,
                        file_name=f"complete_DE_results_{timestamp}.csv",
                        mime="text/csv"
                    )
                    
                    # Results summary
                    with st.expander("Detailed Results Summary"):
                        st.dataframe(results.round(4), use_container_width=True)
                    
                except Exception as e:
                    progress_bar.empty()
                    status_text.empty()
                    st.error(f"Error during analysis: {str(e)}")
        
        except Exception as e:
            st.error(f"Error processing file: {str(e)}")
            st.info("Please ensure your file is properly formatted with genes as rows and samples as columns.")
    
    else:
        # Display instructions when no file is uploaded
        st.info("👆 Please upload a count matrix file to begin analysis.")
        
        st.markdown("""
        ### Expected File Format:
        - **File type**: CSV or TSV
        - **Rows**: Genes (with gene names/IDs in first column)
        - **Columns**: Samples
        - **First row**: Column headers (sample names)
        - **First column**: Gene identifiers
        
        ### Example format:
        ```
        Gene,Sample1,Sample2,Sample3,Sample4
        GeneA,15,20,8,25
        GeneB,100,85,120,95
        GeneC,5,3,7,4
        ```
        
        ### ⚡ Performance Features:
        - **Batch processing** for faster t-test calculations
        - **Vectorized operations** for fold change calculations
        - **Optional caching** for repeated analyses
        - **Progress tracking** during analysis
        """)

# =============================================================================
# PATHWAY ANALYSIS FUNCTIONS (Second Pipeline)
# =============================================================================

# RateLimitedSession and other classes/functions from second pipeline
class RateLimitedSession:
    """Session with built-in rate limiting and retry logic"""
    
    def __init__(self, requests_per_second=2, max_retries=3):
        self.session = requests.Session()
        self.min_interval = 1.0 / requests_per_second
        self.last_request_time = 0
        self.max_retries = max_retries
        
    @safe_api_call
    def get(self, url, **kwargs):
        """Rate-limited GET request with retry logic"""
        for attempt in range(self.max_retries):
            try:
                # Rate limiting
                elapsed = time.time() - self.last_request_time
                if elapsed < self.min_interval:
                    time.sleep(self.min_interval - elapsed)
                
                self.last_request_time = time.time()
                
                # Make request
                response = self.session.get(url, timeout=30, **kwargs)
                response.raise_for_status()
                return response
                
            except requests.exceptions.RequestException as e:
                if attempt == self.max_retries - 1:
                    logger.error(f"Request failed after {self.max_retries} attempts: {e}")
                    raise
                time.sleep(2 ** attempt)  # Exponential backoff
                
    @safe_api_call
    def post(self, url, **kwargs):
        """Rate-limited POST request with retry logic"""
        for attempt in range(self.max_retries):
            try:
                elapsed = time.time() - self.last_request_time
                if elapsed < self.min_interval:
                    time.sleep(self.min_interval - elapsed)
                
                self.last_request_time = time.time()
                response = self.session.post(url, timeout=30, **kwargs)
                response.raise_for_status()
                return response
                
            except requests.exceptions.RequestException as e:
                if attempt == self.max_retries - 1:
                    logger.error(f"Request failed after {self.max_retries} attempts: {e}")
                    raise
                time.sleep(2 ** attempt)

# Global session instances
kegg_session = RateLimitedSession(requests_per_second=1)
ot_session = RateLimitedSession(requests_per_second=3)

# Enhanced Helper Functions
def load_genes_from_any(uploaded_file) -> list[str]:
    """Enhanced gene loading with better error handling"""
    name = (uploaded_file.name or "").lower()

    def _clean_series_to_genes(series: pd.Series) -> list[str]:
        try:
            vals = (
                series.dropna().astype(str)
                .str.replace(r"[,;|\t ]+", "\n", regex=True)
                .str.split("\n").explode()
            )
            vals = vals.astype(str).str.strip().str.upper()
            vals = vals[vals != ""]
            
            seen, out = set(), []
            for v in vals:
                if v and v not in seen and len(v) > 1:  # Filter out single characters
                    seen.add(v)
                    out.append(v)
                if len(out) >= Config.MAX_GENES:
                    break
            return out
        except Exception as e:
            logger.error(f"Error cleaning gene series: {e}")
            return []

    try:
        # Try different file formats
        if name.endswith((".csv", ".csv.gz")):
            df = pd.read_csv(uploaded_file, compression="infer")
        elif name.endswith((".tsv", ".tsv.gz")):
            df = pd.read_csv(uploaded_file, sep="\t", compression="infer")
        elif name.endswith((".xlsx", ".xls")):
            df = pd.read_excel(uploaded_file)
        else:
            df = None

        if isinstance(df, pd.DataFrame) and not df.empty:
            # Look for gene symbol columns
            lower_map = {str(c).lower(): c for c in df.columns}
            target_col = None
            
            # Priority order for column names
            for key in ("gene.symbol", "symbol", "gene_symbol", "gene", "geneid"):
                if key in lower_map:
                    target_col = lower_map[key]
                    break
                    
            if target_col is None:
                target_col = df.columns[0]
                
            return _clean_series_to_genes(df[target_col])
            
    except Exception as e:
        logger.error(f"Error reading structured file: {e}")

    # Fallback to text parsing
    try:
        uploaded_file.seek(0)
        raw = uploaded_file.read()
        text = raw.decode("utf-8", errors="ignore") if isinstance(raw, (bytes, bytearray)) else str(raw)
        return _clean_series_to_genes(pd.Series([text]))
    except Exception as e:
        logger.error(f"Error reading as text: {e}")
        return []

# Enhanced Caching Functions
@st.cache_data(ttl=3600, show_spinner=False)
def kegg_get(path: str) -> str:
    """Enhanced KEGG API call with better error handling"""
    try:
        response = kegg_session.get(f"https://rest.kegg.jp{path}")
        return response.text if response else ""
    except Exception as e:
        logger.error(f"KEGG API error for {path}: {e}")
        return ""

@st.cache_data(ttl=3600, show_spinner=False)
def ncbi_esearch_gene_ids(gene_symbol: str, organism_entrez: str) -> list[str]:
    """Enhanced NCBI gene search"""
    try:
        handle = Entrez.esearch(
            db="gene",
            term=f"{gene_symbol}[Gene] AND {organism_entrez}[Organism]",
            retmode="xml",
            retmax=5  # Limit results
        )
        record = Entrez.read(handle)
        handle.close()
        return record.get("IdList", [])
    except Exception as e:
        logger.error(f"NCBI search error for {gene_symbol}: {e}")
        return []

@st.cache_data(ttl=3600, show_spinner=False)
def ncbi_esummary_description(gene_id: str) -> str:
    """Enhanced NCBI gene description fetch"""
    try:
        handle = Entrez.esummary(db="gene", id=gene_id, retmode="xml")
        raw_xml = handle.read()
        handle.close()
        root = ET.fromstring(raw_xml)
        docsum = root.find(".//DocumentSummary")
        return (docsum.findtext("Description", default="") or "").strip()
    except Exception as e:
        logger.error(f"NCBI summary error for {gene_id}: {e}")
        return "Description unavailable"

@st.cache_data(ttl=3600, show_spinner=False)
def kegg_ncbi_to_kegg_gene_id(ncbi_gene_id: str, kegg_org_prefix: str) -> str | None:
    """Enhanced NCBI to KEGG ID conversion"""
    try:
        txt = kegg_get(f"/conv/genes/ncbi-geneid:{ncbi_gene_id}")
        if not txt.strip():
            return None
            
        for line in txt.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) == 2 and parts[0].endswith(f"{ncbi_gene_id}") and parts[1].startswith(f"{kegg_org_prefix}:"):
                return parts[1].strip()
        return None
    except Exception as e:
        logger.error(f"KEGG conversion error for {ncbi_gene_id}: {e}")
        return None

@st.cache_data(ttl=3600, show_spinner=False)
def kegg_gene_pathways(kegg_gene_id: str) -> list[str]:
    """Enhanced KEGG pathway fetch"""
    try:
        txt = kegg_get(f"/link/pathway/{kegg_gene_id}")
        if not txt.strip():
            return []
            
        pids = []
        for line in txt.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) == 2 and parts[1].startswith("path:"):
                pids.append(parts[1])
        return pids
    except Exception as e:
        logger.error(f"KEGG pathways error for {kegg_gene_id}: {e}")
        return []

@st.cache_data(ttl=3600, show_spinner=False)
def kegg_pathway_name(pathway_id: str) -> str | None:
    """Enhanced KEGG pathway name fetch"""
    try:
        pid = pathway_id.replace("path:", "")
        txt = kegg_get(f"/get/{pid}")
        for line in txt.split("\n"):
            if line.startswith("NAME"):
                return line.replace("NAME", "").strip()
        return None
    except Exception as e:
        logger.error(f"KEGG pathway name error for {pathway_id}: {e}")
        return "Unknown pathway"

# Enhanced OpenTargets Functions
OT_GQL = "https://api.platform.opentargets.org/api/v4/graphql"

def validate_opentargets_response(data: dict) -> bool:
    """Validate OpenTargets API response structure"""
    required_keys = ['data']
    if not all(key in data for key in required_keys):
        return False
    
    # Check for API errors
    if 'errors' in data and data['errors']:
        logger.error(f"OpenTargets API errors: {data['errors']}")
        return False
    
    return True

@st.cache_data(ttl=3600, show_spinner=False)
def ot_query(query: str, variables: dict | None = None) -> dict:
    """Enhanced OpenTargets GraphQL query"""
    try:
        response = ot_session.post(
            OT_GQL, 
            json={"query": query, "variables": variables or {}},
            headers={"Content-Type": "application/json"}
        )
        
        if not response:
            return {}
            
        data = response.json()
        
        if not validate_opentargets_response(data):
            return {}
            
        return data
    except Exception as e:
        logger.error(f"OpenTargets query error: {e}")
        return {}

@st.cache_data(ttl=3600, show_spinner=False)
def ot_target_from_symbol(symbol: str, species: str = "Homo sapiens") -> dict | None:
    """Enhanced target lookup with better matching"""
    query = """
    query FindTarget($q: String!) {
      search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 10}) {
        hits { 
          id 
          name 
          entity 
          object {
            ... on Target {
              approvedSymbol
              biotype
            }
          }
        }
      }
    }
    """
    
    try:
        data = ot_query(query, {"q": symbol})
        hits = (((data or {}).get("data", {})).get("search", {}) or {}).get("hits", [])
        
        # Prioritize exact matches
        for h in hits:
            if h.get("entity") == "target":
                obj = h.get("object", {})
                approved_symbol = obj.get("approvedSymbol", "")
                if approved_symbol.upper() == symbol.upper():
                    return {
                        "id": h.get("id"), 
                        "approvedSymbol": approved_symbol,
                        "biotype": obj.get("biotype")
                    }
        
        # Fallback to first target hit
        for h in hits:
            if h.get("entity") == "target":
                return {"id": h.get("id"), "approvedSymbol": h.get("name")}
                
        return None
    except Exception as e:
        logger.error(f"Target lookup error for {symbol}: {e}")
        return None

@st.cache_data(ttl=3600, show_spinner=False)
def ot_diseases_for_target(ensembl_id: str, size: int = 25) -> pd.DataFrame:
    """Enhanced disease associations fetch"""
    query = """
    query Associations($id: String!, $size: Int!) {
      target(ensemblId: $id) {
        id
        associatedDiseases(page: {size: $size, index: 0}) {
          rows { 
            disease { 
              id 
              name 
              therapeuticAreas { id name }
            } 
            score 
            datatypeScores {
              id
              score
            }
          }
        }
      }
    }
    """
    
    try:
        data = ot_query(query, {"id": ensembl_id, "size": size})
        rows = (((data or {}).get("data", {})).get("target", {}) or {}).get("associatedDiseases", {}).get("rows", [])
        
        out = []
        for r in rows:
            d = r.get("disease", {})
            therapeutic_areas = d.get("therapeuticAreas", [])
            ta_names = "; ".join([ta.get("name", "") for ta in therapeutic_areas if ta.get("name")])
            
            out.append({
                "target": ensembl_id,
                "disease_id": d.get("id"),
                "disease_name": d.get("name"),
                "association_score": r.get("score"),
                "therapeutic_areas": ta_names,
            })
        return pd.DataFrame(out)
    except Exception as e:
        logger.error(f"Disease associations error for {ensembl_id}: {e}")
        return pd.DataFrame(columns=["target", "disease_id", "disease_name", "association_score", "therapeutic_areas"])

@st.cache_data(ttl=3600, show_spinner=False)
def ot_drugs_for_target(ensembl_id: str, size: int = 50) -> pd.DataFrame:
    """FIXED: Enhanced drug suggestions fetch with improved error handling and correct query"""
    query = """
    query KnownDrugs($id: String!, $size: Int!) {
      target(ensemblId: $id) {
        id
        approvedSymbol
        knownDrugs(size: $size) {
          count
          rows {
            phase
            status
            mechanismOfAction
            drug { 
              id 
              name 
              drugType
              maximumClinicalTrialPhase
            }
            disease { 
              id 
              name 
              therapeuticAreas { name }
            }
          }
        }
      }
    }
    """
    
    try:
        data = ot_query(query, {"id": ensembl_id, "size": size})
        
        # Debug: Log the response structure
        logger.info(f"OpenTargets response for {ensembl_id}: {data.keys() if data else 'No data'}")
        
        target_data = ((data or {}).get("data", {})).get("target", {})
        if not target_data:
            logger.warning(f"No target data found for {ensembl_id}")
            return pd.DataFrame()
            
        known_drugs = target_data.get("knownDrugs", {})
        rows = known_drugs.get("rows", [])
        
        logger.info(f"Found {len(rows)} drug rows for target {ensembl_id}")
        
        if not rows:
            return pd.DataFrame()

        # Process data
        out = []
        for r in rows:
            drug_obj = r.get("drug") or {}
            disease_obj = r.get("disease") or {}
            
            therapeutic_areas = disease_obj.get("therapeuticAreas", [])
            ta_names = "; ".join([ta.get("name", "") for ta in therapeutic_areas if ta.get("name")])
            
            # Handle phase conversion safely
            phase = r.get("phase")
            try:
                phase_numeric = int(phase) if phase and phase.isdigit() else 0
            except:
                phase_numeric = 0
                
            max_phase = drug_obj.get("maximumClinicalTrialPhase")
            try:
                max_phase_numeric = int(max_phase) if max_phase and str(max_phase).isdigit() else 0
            except:
                max_phase_numeric = 0
            
            out.append({
                "target": ensembl_id,
                "drug_id": drug_obj.get("id"),
                "drug_name": drug_obj.get("name"),
                "drug_type": drug_obj.get("drugType"),
                "phase": phase,
                "status": r.get("status"),
                "max_phase": max_phase,
                "moa": r.get("mechanismOfAction"),
                "disease_name": disease_obj.get("name"),
                "therapeutic_areas": ta_names,
                "phase_numeric": phase_numeric,
                "max_phase_numeric": max_phase_numeric,
            })
        
        result_df = pd.DataFrame(out)
        logger.info(f"Processed {len(result_df)} drug entries for {ensembl_id}")
        return result_df
    
    except Exception as e:
        logger.error(f"Error fetching drugs for target {ensembl_id}: {e}")
        return pd.DataFrame()

# Enhanced Core Functions
def fetch_gene_metadata_and_kegg(gene_list: list[str], organism_entrez: str,
                                 kegg_org_prefix: str, progress_bar=None):
    """Enhanced metadata fetch with better progress tracking"""
    results = []
    pathway_to_genes = defaultdict(set)
    
    total_genes = len(gene_list)
    
    for i, gene in enumerate(gene_list, start=1):
        if progress_bar:
            progress_bar.progress(i / total_genes, text=f"Processing {gene} ({i}/{total_genes})")
            
        try:
            # NCBI gene search
            ids = ncbi_esearch_gene_ids(gene, organism_entrez)
            if not ids:
                results.append({
                    "Gene": gene, 
                    "NCBI_ID": None, 
                    "Description": "No NCBI match found", 
                    "KEGG_Pathways": None,
                    "Status": "No match"
                })
                continue
                
            gene_id = ids[0]
            description = ncbi_esummary_description(gene_id)
            
            # KEGG conversion
            kegg_id = kegg_ncbi_to_kegg_gene_id(gene_id, kegg_org_prefix)
            if not kegg_id:
                results.append({
                    "Gene": gene, 
                    "NCBI_ID": gene_id, 
                    "Description": description, 
                    "KEGG_Pathways": None,
                    "Status": "No KEGG match"
                })
                continue
                
            # KEGG pathways
            pids = kegg_gene_pathways(kegg_id)
            pathway_pairs = []
            
            for pid in pids:
                name = kegg_pathway_name(pid) or "Unknown"
                pathway_pairs.append(f"{pid.replace('path:', '')} - {name}")
                pathway_to_genes[pid].add(gene)
                
            pathways_str = "; ".join(pathway_pairs) if pathway_pairs else None
            
            results.append({
                "Gene": gene, 
                "NCBI_ID": gene_id, 
                "Description": description, 
                "KEGG_Pathways": pathways_str,
                "Status": f"Found {len(pids)} pathways" if pids else "No pathways"
            })
            
            # Rate limiting
            time.sleep(0.1)
            
        except Exception as e:
            logger.error(f"Error processing {gene}: {e}")
            results.append({
                "Gene": gene, 
                "NCBI_ID": None, 
                "Description": f"Error: {str(e)[:100]}", 
                "KEGG_Pathways": None,
                "Status": "Error"
            })
    
    return pd.DataFrame(results), pathway_to_genes

def compute_enrichment_counts_only(pathway_to_genes: dict) -> pd.DataFrame:
    """Enhanced enrichment computation with better sorting - REMOVED Genes column"""
    if not pathway_to_genes:
        return pd.DataFrame(columns=["Pathway_ID", "Pathway_Name", "Count", "Gene_List"])
        
    rows = []
    for pid, genes in pathway_to_genes.items():
        pathway_name = kegg_pathway_name(pid) or "Unknown pathway"
        gene_list = sorted(list(genes))
        
        rows.append({
            "Pathway_ID": pid.replace("path:", ""),
            "Pathway_Name": pathway_name,
            "Count": len(genes),
            "Gene_List": ", ".join(gene_list[:10]) + ("..." if len(gene_list) > 10 else "")
        })
    
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["Count", "Pathway_Name"], ascending=[False, True]).reset_index(drop=True)
    
    return df

@st.cache_data(ttl=3600, show_spinner=False)
def build_gene_to_ot_target_map(genes: list[str], species: str = "Homo sapiens") -> dict:
    """Enhanced gene to target mapping with progress tracking"""
    g2t = {}
    for g in genes:
        try:
            hit = ot_target_from_symbol(g, species)
            if hit:
                g2t[g] = hit
                logger.info(f"Mapped {g} to {hit.get('id')}")
            else:
                logger.warning(f"No target found for {g}")
            time.sleep(0.05)  # Rate limiting
        except Exception as e:
            logger.error(f"Error mapping {g}: {e}")
    return g2t

def collect_disease_links(gene_to_target: dict) -> pd.DataFrame:
    """Enhanced disease collection with better error handling"""
    frames = []
    
    for g, tgt in gene_to_target.items():
        tid = tgt.get("id")
        if not tid:
            continue
            
        try:
            df = ot_diseases_for_target(tid)
            if not df.empty:
                df.insert(0, "gene", g)
                df.insert(1, "gene_symbol", tgt.get("approvedSymbol", g))
                frames.append(df)
            time.sleep(0.05)
        except Exception as e:
            logger.error(f"Error collecting diseases for {g}: {e}")
    
    if frames:
        combined = pd.concat(frames, ignore_index=True)
        # Remove duplicates and sort by score
        combined = combined.drop_duplicates(subset=["gene", "disease_id"])
        combined = combined.sort_values("association_score", ascending=False)
        return combined
    
    return pd.DataFrame(columns=["gene", "gene_symbol", "target", "disease_id", "disease_name", "association_score", "therapeutic_areas"])

def collect_drug_suggestions(gene_to_target: dict) -> pd.DataFrame:
    """FIXED: Enhanced drug collection with better filtering and error handling"""
    frames = []
    
    for g, tgt in gene_to_target.items():
        tid = tgt.get("id")
        if not tid:
            continue
            
        try:
            df = ot_drugs_for_target(tid)
            if not df.empty:
                df.insert(0, "gene", g)
                df.insert(1, "gene_symbol", tgt.get("approvedSymbol", g))
                frames.append(df)
                logger.info(f"Found {len(df)} drugs for {g}")
            else:
                logger.info(f"No drugs found for {g}")
            time.sleep(0.05)  # Rate limiting
        except Exception as e:
            logger.error(f"Error collecting drugs for {g}: {e}")
    
    if frames:
        combined = pd.concat(frames, ignore_index=True)
        # Remove duplicates
        combined = combined.drop_duplicates(subset=["gene", "drug_id", "disease_name"])
        logger.info(f"Total drugs collected: {len(combined)}")
        return combined
    
    logger.warning("No drugs collected from any target")
    return pd.DataFrame()

# Enhanced Plotly Theme
def apply_plotly_dark_theme(fig: go.Figure) -> go.Figure:
    """Apply consistent dark theme to Plotly figures"""
    fig.update_layout(
        paper_bgcolor="#1a1f2e",
        plot_bgcolor="#252b3a",
        font_color="#ffffff",
        font_family="Inter, sans-serif",
        title_font_size=16,
        title_font_color="#00d4aa",
        legend=dict(
            bgcolor="rgba(26, 31, 46, 0.8)",
            bordercolor="#4a5568",
            borderwidth=1
        ),
        margin=dict(l=50, r=50, t=50, b=50)
    )
    
    # Update axes
    fig.update_xaxes(
        gridcolor="#4a5568",
        linecolor="#4a5568",
        tickcolor="#ffffff"
    )
    fig.update_yaxes(
        gridcolor="#4a5568",
        linecolor="#4a5568",
        tickcolor="#ffffff"
    )
    
    return fig

# Input validation functions
def validate_email(email: str) -> bool:
    """Basic email validation"""
    if not email:
        return False
    pattern = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
    return bool(re.match(pattern, email))

def sanitize_gene_symbol(gene: str) -> str:
    """Sanitize gene symbols"""
    # Remove potentially dangerous characters
    sanitized = ''.join(c for c in gene.upper() if c.isalnum() or c in ['-', '_'])
    return sanitized[:50]  # Limit length

def validate_gene_list(genes: List[str]) -> Tuple[bool, str]:
    """Validate gene list before processing"""
    if not genes:
        return False, "Gene list is empty"
    
    if len(genes) > Config.MAX_GENES:
        return False, f"Too many genes (max {Config.MAX_GENES})"
    
    # Check for invalid characters
    invalid_chars = set()
    for gene in genes:
        sanitized = sanitize_gene_symbol(gene)
        if sanitized != gene:
            invalid_chars.update(set(''.join(c for c in gene if not c.isalnum() and c not in ['-', '_'])))
    
    if invalid_chars:
        return False, f"Invalid characters found: {''.join(invalid_chars)}"
    
    return True, "Valid"

def run_pathway_analysis(genes_from_input=None):
    """Run the pathway enrichment and drug discovery pipeline"""
    
    # Initialize session state
    if 'analysis_results' not in st.session_state:
        st.session_state.analysis_results = {}
    
    # Hero section
    st.markdown(f"""
    <div class="hero" style="display: flex; flex-direction: column; align-items: center; justify-content: center; text-align: center; gap: 1.5rem; padding: 2rem;">
        <div>
            <h1 style="margin-bottom: 0; font-size: 2.3rem; font-weight: 800;
                    background: linear-gradient(135deg, #00d4aa, #667eea);
                    -webkit-background-clip: text; color: transparent;
                    text-align: center;">
                Gene2Therapy
            </h1>
            <p style="margin-top: 6px; color: #b3b8c5; font-size: 1.05rem;">
                Advanced gene analysis pipeline: annotations → enrichment → disease associations → drug repurposing
            </p>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Enhanced sidebar
    with st.sidebar:
        st.markdown('<div class="sidebar-title">🧬 BioContext Analytics</div>', unsafe_allow_html=True)
        st.markdown('<div class="sidebar-tip">Comprehensive gene-to-therapy pipeline which takes DEGs as an input and find pathway enrichment, disease relation and suggest drugs</div>', unsafe_allow_html=True)
        
        st.markdown("---")
        
        st.markdown("**💡 Tips & Best Practices**")
        st.markdown("""
        - **Gene Lists**: Keep under 100 genes for optimal performance
        - **Rate Limits**: Built-in throttling prevents API timeouts  
        - **Caching**: Results cached for 1 hour to improve speed
        - **Formats**: Supports CSV, TSV, XLSX, and plain text
        - **Columns**: Looks for 'Gene.symbol', 'Symbol', or first column
        """)
        
        st.markdown("---")
        
        st.markdown("**🔗 Data Sources**")
        st.markdown("""
        - **NCBI**: Gene annotations and descriptions
        - **KEGG**: Pathway enrichment analysis  
        - **OpenTargets**: Disease associations and drug data
        """)
        
        st.markdown("---")
        
        with st.expander("📊 Analysis Pipeline"):
            st.markdown("""
            1. **Gene Mapping**: NCBI gene database lookup
            2. **Pathway Analysis**: KEGG pathway enrichment
            3. **Disease Links**: OpenTargets association scores
            4. **Drug Discovery**: Therapeutic compound suggestions
            5. **Visualization**: Interactive network and charts
            """)
    
    # If genes are provided from DEGs analysis, skip input section
    if genes_from_input:
        st.info(f"🎯 Using {len(genes_from_input)} genes from DEGs analysis (top 10 upregulated + top 10 downregulated)")
        run_pathway_analysis_with_genes(genes_from_input)
        return
    
    # Input Section
    st.markdown('<div class="section-title"><span class="icon">🔧</span>Configuration & Input</div>', unsafe_allow_html=True)
    st.markdown('<div class="hint">📋 Upload your gene list or paste gene symbols to begin the analysis pipeline. All API calls are optimized with rate limiting and caching.</div>', unsafe_allow_html=True)
    
    # Email input for NCBI
    email = st.text_input(
        "📧 NCBI Entrez Email (Required)", 
        value="", 
        help="Required by NCBI for API access. Your email helps them contact you if there are issues.",
        placeholder="your.email@example.com"
    )
    
    # Validate email
    if email and not validate_email(email):
        st.error("❌ Please enter a valid email address")
        email = ""
    
    if email:
        Entrez.email = email
    
    # Organism selection
    organisms = {
        "🧑 Homo sapiens (Human)": {"entrez": "Homo sapiens", "kegg": "hsa"},
        "🐭 Mus musculus (Mouse)": {"entrez": "Mus musculus", "kegg": "mmu"},
        "🐀 Rattus norvegicus (Rat)": {"entrez": "Rattus norvegicus", "kegg": "rno"},
    }
    
    org_label = st.selectbox("🔬 Select Organism", list(organisms.keys()), index=0)
    organism_entrez = organisms[org_label]["entrez"]
    kegg_org_prefix = organisms[org_label]["kegg"]
    
    # File upload and manual input
    col_upload, col_manual = st.columns([1, 1])
    
    with col_upload:
        uploaded = st.file_uploader(
            "📁 Upload Gene List", 
            type=["csv", "tsv", "txt", "xlsx"],
            help="Supports CSV, TSV, XLSX, and plain text files. Looks for 'Gene.symbol' or 'Symbol' columns."
        )
    
    with col_manual:
        manual_input = st.text_area(
            "✏️ Or Paste Gene Symbols",
            placeholder="TP53, BRCA1, EGFR, MYC, PIK3CA\nAKT1, KRAS, PTEN...",
            height=100,
            help="Enter gene symbols separated by commas, spaces, or new lines"
        )
    
    # Drug filtering options
    st.markdown("---")
    st.markdown("**💊 Drug Analysis Filters**")
    
    col_phase, col_type = st.columns([1, 1])
    
    with col_phase:
        opt_only_phase4 = st.checkbox(
            "✅ Show only approved drugs (Phase 4+)",
            value=False,
            help="Filter to show only drugs that have completed clinical trials"
        )
    
    with col_type:
        show_investigational = st.checkbox(
            "🧪 Include investigational compounds",
            value=False,
            help="Include drugs in earlier clinical trial phases"
        )
    
    # Process input
    genes_from_user: list[str] = []
    
    if manual_input.strip():
        genes_from_user = (
            pd.Series([manual_input])
            .str.replace(r"[,;|\t\n ]+", "\n", regex=True)
            .str.split("\n").explode().str.strip().str.upper()
        )
        genes_from_user = [g for g in genes_from_user.tolist() if g and len(g) > 1][:Config.MAX_GENES]
        
    elif uploaded is not None:
        genes_from_user = load_genes_from_any(uploaded)
    
    # Validate gene list
    if genes_from_user:
        is_valid, validation_msg = validate_gene_list(genes_from_user)
        if not is_valid:
            st.warning(f"⚠️ {validation_msg}")
            if len(genes_from_user) > Config.MAX_GENES:
                genes_from_user = genes_from_user[:Config.MAX_GENES]
                st.info(f"Using first {Config.MAX_GENES} genes")
    
    # Analysis button
    st.markdown("---")
    
    col_btn, col_info = st.columns([1, 2])
    
    with col_btn:
        run_analysis = st.button(
            "🚀 Start Analysis", 
            type="primary", 
            disabled=(not genes_from_user or not email),
            use_container_width=True
        )
    
    with col_info:
        if genes_from_user:
            st.success(f"✅ {len(genes_from_user)} genes ready for analysis")
        elif not email:
            st.warning("⚠️ Please provide your email address")
        else:
            st.info("📝 Please upload a file or paste gene symbols")
    
    if run_analysis:
        run_pathway_analysis_with_genes(genes_from_user, organism_entrez, kegg_org_prefix, 
                                       opt_only_phase4, show_investigational)

def run_pathway_analysis_with_genes(genes_from_input, organism_entrez=None, 
                                   kegg_org_prefix=None, opt_only_phase4=False, 
                                   show_investigational=False):
    """Run the actual pathway analysis with provided genes"""
    
    if organism_entrez is None:
        organism_entrez = "Homo sapiens"
    if kegg_org_prefix is None:
        kegg_org_prefix = "hsa"
    
    st.markdown("---")
    
    # Create tabs
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📇 Gene Metadata", 
        "📊 Pathway Enrichment", 
        "🧬 Disease Associations", 
        "💊 Drug Suggestions", 
        "🌐 Network Visualization"
    ])
    
    # Step 1: Gene Metadata & KEGG
    with tab1:
        st.markdown('<div class="section-title"><span class="icon">📇</span>Gene Annotations & KEGG Pathways</div>', unsafe_allow_html=True)
        
        if 'metadata' not in st.session_state.analysis_results:
            progress_container = st.container()
            with progress_container:
                progress_bar = st.progress(0.0)
                status_text = st.empty()
                
            with st.spinner("🔍 Fetching gene metadata and pathway information..."):
                df_meta, pathway_to_genes = fetch_gene_metadata_and_kegg(
                    genes_from_input, organism_entrez, kegg_org_prefix, progress_bar
                )
                
            st.session_state.analysis_results['metadata'] = df_meta
            st.session_state.analysis_results['pathways'] = pathway_to_genes
            
            progress_container.empty()
        else:
            df_meta = st.session_state.analysis_results['metadata']
            pathway_to_genes = st.session_state.analysis_results['pathways']
        
        if not df_meta.empty:
            # Summary metrics
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric("Total Genes", len(df_meta))
            with col2:
                found_genes = len(df_meta[df_meta['NCBI_ID'].notna()])
                st.metric("NCBI Matches", found_genes)
            with col3:
                kegg_genes = len(df_meta[df_meta['KEGG_Pathways'].notna()])
                st.metric("KEGG Mapped", kegg_genes)
            with col4:
                total_pathways = len(pathway_to_genes)
                st.metric("Unique Pathways", total_pathways)
            
            st.markdown("---")
            
            # Display results table
            display_df = df_meta.copy()
            display_df.insert(0, "#", range(1, len(display_df) + 1))
            
            st.dataframe(
                display_df, 
                use_container_width=True, 
                hide_index=True,
                column_config={
                    "Description": st.column_config.TextColumn(width="medium"),
                    "KEGG_Pathways": st.column_config.TextColumn(width="large"),
                    "Status": st.column_config.TextColumn(width="small")
                }
            )
            
            # Download button
            csv_data = df_meta.to_csv(index=False).encode("utf-8")
            st.download_button(
                "⬇️ Download Gene Metadata",
                data=csv_data,
                file_name=f"gene_metadata_{organism_entrez.replace(' ', '_').lower()}.csv",
                mime="text/csv"
            )
        else:
            st.error("❌ No gene metadata could be retrieved. Please check your gene symbols and try again.")
    
    # Step 2: Pathway Enrichment
    with tab2:
        st.markdown('<div class="section-title"><span class="icon">📊</span>KEGG Pathway Enrichment</div>', unsafe_allow_html=True)
        
        if 'pathways' in st.session_state.analysis_results:
            pathway_to_genes = st.session_state.analysis_results['pathways']
            
            with st.spinner("📈 Computing pathway enrichment..."):
                df_enrich = compute_enrichment_counts_only(pathway_to_genes)
            
            if not df_enrich.empty:
                # Summary metrics
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.metric("Enriched Pathways", len(df_enrich))
                with col2:
                    max_genes = df_enrich['Count'].max() if not df_enrich.empty else 0
                    st.metric("Max Genes/Pathway", max_genes)
                with col3:
                    avg_genes = df_enrich['Count'].mean() if not df_enrich.empty else 0
                    st.metric("Avg Genes/Pathway", f"{avg_genes:.1f}")
                
                st.markdown("---")
                
                # Display enrichment table - REMOVED Genes column
                display_enrich = df_enrich[['Pathway_ID', 'Pathway_Name', 'Count', 'Gene_List']].copy()
                display_enrich.insert(0, "#", range(1, len(display_enrich) + 1))
                
                st.dataframe(
                    display_enrich, 
                    use_container_width=True, 
                    hide_index=True,
                    column_config={
                        "Pathway_Name": st.column_config.TextColumn(width="large"),
                        "Gene_List": st.column_config.TextColumn(width="large"),
                        "Count": st.column_config.NumberColumn(width="small")
                    }
                )
                
                # Visualization
                if len(df_enrich) > 0:
                    st.markdown("---")
                    st.markdown("**📊 Top Pathways Visualization**")
                    
                    top_pathways = df_enrich.head(15)
                    
                    fig = px.bar(
                        top_pathways, 
                        x="Count", 
                        y="Pathway_Name", 
                        orientation="h",
                        title="Top 15 Pathways by Gene Count",
                        labels={"Count": "Number of Genes", "Pathway_Name": "KEGG Pathway"},
                        color="Count",
                        color_continuous_scale="viridis"
                    )
                    
                    fig.update_layout(height=600, yaxis={'categoryorder':'total ascending'})
                    fig = apply_plotly_dark_theme(fig)
                    
                    st.plotly_chart(fig, use_container_width=True)
                
                # Download button - also remove Genes column from download
                download_enrich = df_enrich[['Pathway_ID', 'Pathway_Name', 'Count', 'Gene_List']].copy()
                csv_data = download_enrich.to_csv(index=False).encode("utf-8")
                st.download_button(
                    "⬇️ Download Enrichment Results",
                    data=csv_data,
                    file_name="pathway_enrichment_analysis.csv",
                    mime="text/csv"
                )
            else:
                st.info("ℹ️ No pathway enrichment found for the provided genes.")
    
    # Step 3: Disease Associations
    with tab3:
        st.markdown('<div class="section-title"><span class="icon">🧬</span>Disease Associations</div>', unsafe_allow_html=True)
        
        if 'diseases' not in st.session_state.analysis_results:
            with st.spinner("🔍 Mapping genes to targets and fetching disease associations..."):
                gene_to_target = build_gene_to_ot_target_map(genes_from_input)
                df_diseases = collect_disease_links(gene_to_target)
                
            st.session_state.analysis_results['gene_to_target'] = gene_to_target
            st.session_state.analysis_results['diseases'] = df_diseases
        else:
            gene_to_target = st.session_state.analysis_results['gene_to_target']
            df_diseases = st.session_state.analysis_results['diseases']
        
        if not df_diseases.empty:
            # Summary metrics
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                mapped_genes = len(gene_to_target)
                st.metric("Mapped Genes", mapped_genes)
            with col2:
                unique_diseases = df_diseases['disease_name'].nunique()
                st.metric("Associated Diseases", unique_diseases)
            with col3:
                avg_score = df_diseases['association_score'].mean()
                st.metric("Avg Association Score", f"{avg_score:.3f}")
            with col4:
                high_conf = len(df_diseases[df_diseases['association_score'] > 0.5])
                st.metric("High Confidence (>0.5)", high_conf)
            
            st.markdown("---")
            
            # Disease summary
            disease_summary = (
                df_diseases.groupby(['disease_id', 'disease_name', 'therapeutic_areas'])
                .agg({
                    'gene': lambda x: len(set(x)),
                    'association_score': ['max', 'mean']
                })
                .round(3)
                .reset_index()
            )
            
            disease_summary.columns = ['Disease_ID', 'Disease_Name', 'Therapeutic_Areas', 'Gene_Count', 'Max_Score', 'Avg_Score']
            disease_summary = disease_summary.sort_values(['Gene_Count', 'Max_Score'], ascending=[False, False])
            
            # Display disease summary
            display_diseases = disease_summary.copy()
            display_diseases.insert(0, "#", range(1, len(display_diseases) + 1))
            
            st.dataframe(
                display_diseases, 
                use_container_width=True, 
                hide_index=True,
                column_config={
                    "Disease_Name": st.column_config.TextColumn(width="large"),
                    "Therapeutic_Areas": st.column_config.TextColumn(width="medium"),
                    "Max_Score": st.column_config.NumberColumn(format="%.3f"),
                    "Avg_Score": st.column_config.NumberColumn(format="%.3f")
                }
            )
            
            # Visualization
            if len(disease_summary) > 0:
                st.markdown("---")
                st.markdown("**🎯 Top Disease Associations**")
                
                top_diseases = disease_summary.head(20)
                
                fig = px.scatter(
                    top_diseases,
                    x="Gene_Count",
                    y="Max_Score", 
                    size="Avg_Score",
                    hover_name="Disease_Name",
                    title="Disease Associations: Gene Count vs Association Score",
                    labels={
                        "Gene_Count": "Number of Associated Genes",
                        "Max_Score": "Maximum Association Score",
                        "Avg_Score": "Average Association Score"
                    }
                )
                
                fig = apply_plotly_dark_theme(fig)
                st.plotly_chart(fig, use_container_width=True)
            
            # Download buttons
            col_dl1, col_dl2 = st.columns(2)
            
            with col_dl1:
                csv_detailed = df_diseases.to_csv(index=False).encode("utf-8")
                st.download_button(
                    "⬇️ Download Detailed Associations",
                    data=csv_detailed,
                    file_name="disease_associations_detailed.csv",
                    mime="text/csv"
                )
            
            with col_dl2:
                csv_summary = disease_summary.to_csv(index=False).encode("utf-8")
                st.download_button(
                    "⬇️ Download Disease Summary",
                    data=csv_summary,
                    file_name="disease_associations_summary.csv",
                    mime="text/csv"
                )
        else:
            st.info("ℹ️ No disease associations found. This may occur with non-human genes or if the genes are not well-characterized.")
    
    # Step 4: Drug Suggestions - FIXED with corrected Avg Phase calculation
    with tab4:
        st.markdown('<div class="section-title"><span class="icon">💊</span>Therapeutic Drug Suggestions</div>', unsafe_allow_html=True)
        
        if 'drugs' not in st.session_state.analysis_results:
            if 'gene_to_target' in st.session_state.analysis_results:
                gene_to_target = st.session_state.analysis_results['gene_to_target']
                with st.spinner("💊 Collecting drug suggestions from OpenTargets..."):
                    df_drugs = collect_drug_suggestions(gene_to_target)
                st.session_state.analysis_results['drugs'] = df_drugs
            else:
                st.error("❌ Gene-to-target mapping required. Please complete the Disease Associations step first.")
                df_drugs = pd.DataFrame()
        else:
            df_drugs = st.session_state.analysis_results['drugs']
        
        if not df_drugs.empty:
            # Apply filters
            filtered_drugs = df_drugs.copy()
            
            # Debug info
            st.markdown(f"**📊 Found {len(filtered_drugs)} total drug entries**")
            
            if opt_only_phase4:
                before_filter = len(filtered_drugs)
                filtered_drugs = filtered_drugs[
                    (filtered_drugs['phase_numeric'] >= 4) | 
                    (filtered_drugs['max_phase_numeric'] >= 4)
                ]
                st.markdown(f"**🔍 Phase 4+ filter: {before_filter} → {len(filtered_drugs)} entries**")
            
            if not show_investigational:
                before_filter = len(filtered_drugs)
                filtered_drugs = filtered_drugs[
                    filtered_drugs['status'] != 'Investigational'
                ]
                st.markdown(f"**🔍 Non-investigational filter: {before_filter} → {len(filtered_drugs)} entries**")
            
            if not filtered_drugs.empty:
                # FIXED: Correct Avg Phase calculation - only consider valid phases > 0
                valid_phases = filtered_drugs[filtered_drugs['phase'] > 0]['phase']
                avg_phase = valid_phases.mean() if len(valid_phases) > 0 else 0
                
                # Summary metrics with CORRECTED Avg Phase
                col1, col2, col3, col4 = st.columns(4)
                
                with col1:
                    unique_drugs = filtered_drugs['drug_name'].nunique()
                    st.metric("Unique Drugs", unique_drugs)
                with col2:
                    unique_targets = filtered_drugs['gene'].nunique()
                    st.metric("Targeted Genes", unique_targets)
                with col3:
                    # CORRECTED: Show actual average phase
                    st.metric("Avg Phase", f"{avg_phase:.1f}")
                with col4:
                    approved_drugs = len(filtered_drugs[filtered_drugs['phase_numeric'] >= 4])
                    st.metric("Approved Drugs", approved_drugs)
                
                st.markdown("---")
                
                # Create display dataframe without phase_numeric and max_phase_numeric columns
                display_columns = [
                    'gene', 'gene_symbol', 'drug_id', 'drug_name', 'drug_type', 
                    'phase', 'status', 'max_phase', 'moa', 'disease_name', 'therapeutic_areas'
                ]
                
                # Filter to only include columns that exist in the dataframe
                available_columns = [col for col in display_columns if col in filtered_drugs.columns]
                display_drugs = filtered_drugs[available_columns].copy()
                display_drugs.insert(0, "#", range(1, len(display_drugs) + 1))
                
                st.dataframe(
                    display_drugs, 
                    use_container_width=True,
                    column_config={
                        "#": st.column_config.NumberColumn(width="small"),
                        "gene": st.column_config.TextColumn("Gene", width="small"),
                        "gene_symbol": st.column_config.TextColumn("Gene Symbol", width="small"),
                        "drug_name": st.column_config.TextColumn("Drug Name", width="medium"),
                        "drug_type": st.column_config.TextColumn("Drug Type", width="small"),
                        "phase": st.column_config.TextColumn("Phase", width="small"),
                        "status": st.column_config.TextColumn("Status", width="small"),
                        "max_phase": st.column_config.TextColumn("Max Phase", width="small"),
                        "moa": st.column_config.TextColumn("Mechanism of Action", width="large"),
                        "disease_name": st.column_config.TextColumn("Disease", width="medium"),
                        "therapeutic_areas": st.column_config.TextColumn("Therapeutic Areas", width="medium"),
                    }
                )
                
                # Visualization
                if len(filtered_drugs) > 0:
                    st.markdown("---")
                    st.markdown("**📈 Drug Development Phase Distribution**")
                    
                    phase_counts = filtered_drugs['phase'].value_counts().reset_index()
                    phase_counts.columns = ['Phase', 'Count']
                    phase_counts = phase_counts.sort_values('Phase')
                    
                    fig = px.bar(
                        phase_counts,
                        x='Phase',
                        y='Count',
                        title='Drugs by Clinical Trial Phase',
                        labels={'Phase': 'Clinical Trial Phase', 'Count': 'Number of Drugs'},
                        color='Count',
                        color_continuous_scale='viridis'
                    )
                    
                    fig = apply_plotly_dark_theme(fig)
                    st.plotly_chart(fig, use_container_width=True)
                
                # Download button - also remove numeric columns from download
                download_columns = [col for col in display_columns if col in filtered_drugs.columns]
                download_data = filtered_drugs[download_columns].to_csv(index=False).encode("utf-8")
                st.download_button(
                    "⬇️ Download Drug Suggestions",
                    data=download_data,
                    file_name="drug_suggestions.csv",
                    mime="text/csv"
                )
            else:
                st.info("ℹ️ No drugs match the current filter criteria. Try adjusting the filters or check if your genes have known drug associations.")
        else:
            st.info("ℹ️ No drug suggestions available. This could happen if:")
            st.markdown("""
            - The mapped genes don't have known drugs in OpenTargets
            - The genes are not well-characterized for drug targeting
            - There are temporary API issues
            - Try using more common cancer genes (TP53, EGFR, BRAF, etc.) for testing
            """)

    # Step 5: Network Visualization
    with tab5:
        st.markdown('<div class="section-title"><span class="icon">🌐</span>Interactive Network Visualization</div>', unsafe_allow_html=True)
        
        # Check if we have the required data
        has_diseases = 'diseases' in st.session_state.analysis_results and not st.session_state.analysis_results['diseases'].empty
        has_drugs = 'drugs' in st.session_state.analysis_results and not st.session_state.analysis_results['drugs'].empty
        has_pathways = 'pathways' in st.session_state.analysis_results and st.session_state.analysis_results['pathways']
        
        if has_diseases or has_drugs or has_pathways:
            col_left, col_right = st.columns(2)
            
            # Gene-Disease-Drug Network
            with col_left:
                if has_diseases:
                    st.markdown("**🔗 Gene-Disease-Drug Network**")
                    
                    try:
                        df_diseases = st.session_state.analysis_results['diseases']
                        
                        # Get top diseases by gene count
                        top_diseases = (
                            df_diseases.groupby('disease_name')
                            .agg({'gene': lambda x: len(set(x))})
                            .sort_values('gene', ascending=False)
                            .head(10)
                            .index.tolist()
                        )
                        
                        # Filter data
                        disease_subset = df_diseases[df_diseases['disease_name'].isin(top_diseases)]
                        genes_in_network = sorted(set(disease_subset['gene']))
                        
                        # Add drug connections if available
                        drug_connections = []
                        if has_drugs:
                            df_drugs = st.session_state.analysis_results['drugs']
                            drug_subset = df_drugs[
                                (df_drugs['gene'].isin(genes_in_network))
                            ].head(15)  # Limit drugs
                            
                            for _, row in drug_subset.iterrows():
                                drug_connections.append((f"Gene: {row['gene']}", f"Drug: {row['drug_name']}"))
                        
                        # Create network
                        nodes = (
                            [f"Gene: {g}" for g in genes_in_network] + 
                            [f"Disease: {d}" for d in top_diseases] +
                            [f"Drug: {drug.split(': ')[1]}" for _, drug in drug_connections]
                        )
                        
                        node_index = {n: i for i, n in enumerate(nodes)}
                        
                        # Create links
                        sources, targets, values = [], [], []
                        
                        # Gene-Disease links
                        for _, row in disease_subset.iterrows():
                            gene_node = f"Gene: {row['gene']}"
                            disease_node = f"Disease: {row['disease_name']}"
                            if gene_node in node_index and disease_node in node_index:
                                sources.append(node_index[gene_node])
                                targets.append(node_index[disease_node])
                                values.append(max(1, int(row['association_score'] * 10)))
                        
                        # Gene-Drug links
                        for gene_drug, drug_node in drug_connections:
                            if gene_drug in node_index and drug_node in node_index:
                                sources.append(node_index[gene_drug])
                                targets.append(node_index[drug_node])
                                values.append(3)
                        
                        if sources:
                            # Create Sankey diagram
                            fig_sankey = go.Figure(data=[go.Sankey(
                                node=dict(
                                    pad=15,
                                    thickness=20,
                                    line=dict(color="#4a5568", width=2),
                                    label=nodes,
                                    color="#2d3748"
                                ),
                                link=dict(
                                    source=sources,
                                    target=targets,
                                    value=values,
                                    color="rgba(0, 212, 170, 0.3)"
                                )
                            )])
                            
                            fig_sankey.update_layout(
                                title_text="Gene → Disease → Drug Connections",
                                height=600
                            )
                            fig_sankey = apply_plotly_dark_theme(fig_sankey)
                            
                            st.plotly_chart(fig_sankey, use_container_width=True)
                        else:
                            st.info("Not enough connections for network visualization")
                            
                    except Exception as e:
                        st.error(f"Error creating network: {e}")
            
            # Pathway Network
            with col_right:
                if has_pathways:
                    st.markdown("**🛤️ Gene-Pathway Network**")
                    
                    try:
                        pathway_to_genes = st.session_state.analysis_results['pathways']
                        
                        # Get top pathways
                        top_pathways = sorted(
                            pathway_to_genes.items(), 
                            key=lambda x: len(x[1]), 
                            reverse=True
                        )[:8]
                        
                        # Create network graph
                        G = nx.Graph()
                        
                        for pathway_id, genes in top_pathways:
                            pathway_name = kegg_pathway_name(pathway_id) or pathway_id.replace("path:", "")
                            pathway_name = pathway_name[:30] + "..." if len(pathway_name) > 30 else pathway_name
                            
                            for gene in genes:
                                G.add_edge(f"P: {pathway_name}", f"G: {gene}")
                        
                        if G.nodes():
                            # Calculate layout
                            pos = nx.spring_layout(G, k=1, iterations=50, seed=42)
                            
                            # Prepare edge traces
                            edge_x, edge_y = [], []
                            for edge in G.edges():
                                x0, y0 = pos[edge[0]]
                                x1, y1 = pos[edge[1]]
                                edge_x.extend([x0, x1, None])
                                edge_y.extend([y0, y1, None])
                            
                            # Prepare node traces
                            node_x = [pos[node][0] for node in G.nodes()]
                            node_y = [pos[node][1] for node in G.nodes()]
                            node_text = list(G.nodes())
                            node_colors = ['#00d4aa' if n.startswith('P:') else '#667eea' for n in node_text]
                            
                            # Create figure
                            fig_network = go.Figure()
                            
                            # Add edges
                            fig_network.add_trace(go.Scatter(
                                x=edge_x, y=edge_y,
                                line=dict(width=1, color='rgba(255,255,255,0.3)'),
                                hoverinfo='none',
                                mode='lines'
                            ))
                            
                            # Add nodes
                            fig_network.add_trace(go.Scatter(
                                x=node_x, y=node_y,
                                mode='markers+text',
                                hoverinfo='text',
                                text=[n.split(': ')[1] for n in node_text],
                                textposition="middle center",
                                hovertext=node_text,
                                marker=dict(
                                    size=10,
                                    color=node_colors,
                                    line=dict(width=2, color='white')
                                )
                            ))
                            
                            fig_network.update_layout(
                                title="Gene-Pathway Interaction Network",
                                showlegend=False,
                                hovermode='closest',
                                margin=dict(b=20,l=5,r=5,t=40),
                                annotations=[
                                    dict(
                                        text="Pathways (green) connected to genes (blue)",
                                        showarrow=False,
                                        xref="paper", yref="paper",
                                        x=0.005, y=-0.002,
                                        xanchor='left', yanchor='bottom',
                                        font=dict(color='#b3b8c5', size=12)
                                    )
                                ],
                                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                height=600
                            )
                            
                            fig_network = apply_plotly_dark_theme(fig_network)
                            st.plotly_chart(fig_network, use_container_width=True)
                        else:
                            st.info("No pathway connections available for visualization")
                            
                    except Exception as e:
                        st.error(f"Error creating pathway network: {e}")
        else:
            st.info("🔄 Complete the previous analysis steps to generate network visualizations.")
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: #8b92a5; font-size: 0.9rem; padding: 1rem;'>
        <p><strong>Gene2Therapy</strong> • Enhanced dark mode interface with optimized API integration</p>
        <p>Data sources: NCBI E-utilities, KEGG REST API, OpenTargets GraphQL • Results cached for 1 hour</p>
        <p>⚠️ Always validate findings with primary literature and clinical sources</p>
    </div>
    """, unsafe_allow_html=True)

# =============================================================================
# MAIN APP INTEGRATION
# =============================================================================

def main():
    # Initialize session state
    initialize_session_state()
    
    # Main title and mode selection
    st.title("🧬 Integrated Gene Analysis Pipeline")
    
    st.markdown("""
    This integrated pipeline provides a complete workflow from differential expression analysis 
    to pathway enrichment, disease association, and drug discovery.
    """)
    
    # Analysis mode selection with improved styling - REMOVED the third option
    st.header("Select Analysis Mode")
    
    col1, col2 = st.columns(2)  # Changed from 3 columns to 2
    
    with col1:
        # Create custom button with active state styling
        button_class = "mode-button active" if st.session_state.analysis_mode == 'degs_only' else "mode-button"
        st.markdown(f'<div class="{button_class}">🧬<br><strong>DEGs Analysis Only</strong><br><small>Identify differentially expressed genes from RNA-seq data</small></div>', 
                   unsafe_allow_html=True)
        if st.button("Select DEGs Analysis Only", key="degs_btn", use_container_width=True):
            st.session_state.analysis_mode = 'degs_only'
            st.session_state.current_pipeline_step = 'degs'
            st.session_state.show_pathway_analysis = False
            st.session_state.degs_completed = False
            st.rerun()
    
    with col2:
        button_class = "mode-button active" if st.session_state.analysis_mode == 'pathway_only' else "mode-button"
        st.markdown(f'<div class="{button_class}">💊<br><strong>Pathway & Drug Analysis Only</strong><br><small>Enrichment, disease links, and drug discovery</small></div>', 
                   unsafe_allow_html=True)
        if st.button("Select Pathway & Drug Analysis Only", key="pathway_btn", use_container_width=True):
            st.session_state.analysis_mode = 'pathway_only'
            st.session_state.current_pipeline_step = 'pathway'
            st.session_state.show_pathway_analysis = True
            st.session_state.degs_completed = False
            st.rerun()
    
    # Display appropriate pipeline based on selection - REMOVED the 'both' option logic
    if st.session_state.analysis_mode == 'degs_only':
        st.markdown("---")
        st.info("🔬 Running Differential Expression Analysis Only")
        run_degs_analysis()
    
    elif st.session_state.analysis_mode == 'pathway_only':
        st.markdown("---")
        st.info("💊 Running Pathway Enrichment & Drug Discovery Only")
        run_pathway_analysis()
    
    else:
        # Show instructions when no mode is selected
        st.markdown("---")
        st.subheader("Pipeline Overview")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("### 🧬 DEGs Analysis Pipeline")
            st.markdown("""
            - **Input**: RNA-seq count matrix
            - **Process**: Differential expression analysis
            - **Output**: 
                - Up/down regulated genes
                - Volcano plots
                - Statistical results
            - **Use Case**: Identify significantly expressed genes between conditions
            """)
        
        with col2:
            st.markdown("### 💊 Pathway & Drug Analysis Pipeline")
            st.markdown("""
            - **Input**: Gene list (from DEGs or manual input)
            - **Process**: 
                - Pathway enrichment (KEGG)
                - Disease associations (OpenTargets)
                - Drug suggestions
            - **Output**: 
                - Annotated pathways
                - Disease links
                - Therapeutic compounds
            - **Use Case**: Functional interpretation and drug repurposing
            """)
        
        st.markdown("---")
        st.info("👆 Select an analysis mode above to get started!")

if __name__ == "__main__":
    main()