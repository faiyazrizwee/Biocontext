import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
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
    MAX_DISPLAY_ROWS = 500
    MAX_NETWORK_NODES = 100

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
    page_title="Gene2Therapy - Integrated Gene Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# =============================================================================
# MODERN LIGHT THEME CSS
# =============================================================================

def get_light_theme_css():
    return """
<style>
  @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=Source+Sans+Pro:wght@400;600;700&display=swap');

  :root {
    --bg: #ffffff;
    --panel: #f8fafc;
    --surface: #ffffff;
    --card: #ffffff;
    --text: #1e293b;
    --text-muted: #64748b;
    --text-secondary: #94a3b8;
    --border: #e2e8f0;
    --border-strong: #cbd5e1;
    --accent: #3b82f6;
    --accent-hover: #2563eb;
    --accent-active: #1d4ed8;
    --accent-light: #dbeafe;
    --secondary: #10b981;
    --secondary-hover: #059669;
    --secondary-light: #d1fae5;
    --danger: #ef4444;
    --danger-light: #fee2e2;
    --warning: #f59e0b;
    --warning-light: #fef3c7;
    --success: #10b981;
    --success-light: #d1fae5;
    --info: #0ea5e9;
    --info-light: #e0f2fe;
    --input-bg: #ffffff;
    --shadow: rgba(0, 0, 0, 0.05);
    --shadow-strong: rgba(0, 0, 0, 0.1);
    --gradient: linear-gradient(135deg, #3b82f6 0%, #10b981 100%);
    --gradient-light: linear-gradient(135deg, #dbeafe 0%, #d1fae5 100%);
  }

  /* Global Styles */
  .stApp {
    background: var(--bg) !important;
    color: var(--text);
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    line-height: 1.6;
  }

  /* Smooth transitions */
  * {
    transition: all 0.2s ease;
  }

  /* Header */
  header[data-testid="stHeader"] {
    background: var(--bg) !important;
    border-bottom: 1px solid var(--border);
    backdrop-filter: blur(10px);
  }

  .block-container {
    padding-top: 2rem !important;
    padding-bottom: 2rem;
    max-width: 1200px;
  }

  /* Sidebar */
  section[data-testid="stSidebar"] {
    background: linear-gradient(180deg, var(--panel) 0%, #f1f5f9 100%) !important;
    border-right: 1px solid var(--border);
    box-shadow: 2px 0 12px var(--shadow);
  }

  section[data-testid="stSidebar"] * {
    color: var(--text) !important;
  }

  /* Hero Section */
  .hero {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 20px;
    padding: 2.5rem;
    margin: 1rem 0 2.5rem 0;
    box-shadow: 0 4px 20px var(--shadow);
    position: relative;
    overflow: hidden;
  }

  .hero::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    height: 4px;
    background: var(--gradient);
  }

  .hero-content h1 {
    font-size: 2.75rem;
    font-weight: 800;
    margin: 0 0 1rem 0;
    background: var(--gradient);
    -webkit-background-clip: text;
    background-clip: text;
    color: transparent;
    line-height: 1.2;
    font-family: 'Source Sans Pro', sans-serif;
  }

  .hero-content p {
    color: var(--text-muted);
    font-size: 1.1rem;
    margin: 0;
    line-height: 1.6;
    max-width: 600px;
  }

  /* Section Titles */
  .section-title {
    font-size: 1.5rem;
    font-weight: 700;
    color: var(--text);
    margin: 2.5rem 0 1.5rem 0;
    padding-bottom: 0.75rem;
    border-bottom: 2px solid var(--border);
    position: relative;
    font-family: 'Source Sans Pro', sans-serif;
  }

  .section-title::after {
    content: '';
    position: absolute;
    bottom: -2px;
    left: 0;
    width: 80px;
    height: 2px;
    background: var(--gradient);
  }

  /* Cards */
  .card {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 16px;
    padding: 1.5rem;
    margin-bottom: 1.5rem;
    box-shadow: 0 2px 8px var(--shadow);
    transition: transform 0.2s ease, box-shadow 0.2s ease;
  }

  .card:hover {
    transform: translateY(-2px);
    box-shadow: 0 8px 25px var(--shadow-strong);
  }

  /* Mode Selection Buttons */
  .mode-button {
    background: var(--surface);
    border: 2px solid var(--border);
    border-radius: 16px;
    padding: 2rem 1.5rem;
    color: var(--text);
    font-weight: 600;
    font-size: 1.1rem;
    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
    height: 140px;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    text-align: center;
    margin-bottom: 1.5rem;
    position: relative;
    overflow: hidden;
  }

  .mode-button:hover {
    border-color: var(--accent);
    transform: translateY(-4px);
    box-shadow: 0 12px 30px rgba(59, 130, 246, 0.15);
  }

  .mode-button:active {
    transform: translateY(-1px);
  }

  .mode-button.active {
    border-color: var(--accent);
    background: linear-gradient(135deg, #f0f9ff 0%, #f0fdf4 100%);
    box-shadow: 0 8px 25px rgba(59, 130, 246, 0.2);
  }

  .mode-button.active::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    height: 3px;
    background: var(--gradient);
  }

  .mode-icon {
    font-size: 2.5rem;
    margin-bottom: 1rem;
    display: block;
  }

  /* Input Fields */
  .stTextInput > div > div,
  .stTextArea > div > div {
    background: var(--input-bg) !important;
    border: 2px solid var(--border) !important;
    border-radius: 12px !important;
    transition: all 0.2s ease;
  }

  .stTextInput > div > div:focus-within,
  .stTextArea > div > div:focus-within {
    border-color: var(--accent) !important;
    box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1) !important;
  }

  .stTextInput input,
  .stTextArea textarea {
    background: transparent !important;
    color: var(--text) !important;
    font-family: 'Inter', sans-serif;
    font-size: 0.95rem;
  }

  .stTextInput input::placeholder,
  .stTextArea textarea::placeholder {
    color: var(--text-secondary) !important;
  }

  /* File Uploader */
  .stFileUploader [data-testid="stFileUploaderDropzone"] {
    background: var(--surface) !important;
    border: 2px dashed var(--border) !important;
    border-radius: 16px !important;
    padding: 3rem 2rem !important;
    transition: all 0.3s ease;
  }

  .stFileUploader [data-testid="stFileUploaderDropzone"]:hover {
    border-color: var(--accent) !important;
    background: var(--accent-light) !important;
    border-style: solid !important;
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
    background: var(--surface) !important;
    border: 1px solid var(--border) !important;
    border-radius: 12px !important;
    box-shadow: 0 4px 20px var(--shadow) !important;
  }

  /* Buttons */
  .stButton > button {
    background: var(--gradient) !important;
    color: white !important;
    border: none !important;
    border-radius: 12px !important;
    padding: 0.875rem 2rem !important;
    font-weight: 600 !important;
    font-size: 1rem !important;
    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1) !important;
    box-shadow: 0 4px 15px rgba(59, 130, 246, 0.3) !important;
    position: relative;
    overflow: hidden;
  }

  .stButton > button:hover {
    transform: translateY(-2px) !important;
    box-shadow: 0 8px 25px rgba(59, 130, 246, 0.4) !important;
  }

  .stButton > button:active {
    transform: translateY(0) !important;
  }

  .stButton > button::after {
    content: '';
    position: absolute;
    top: 50%;
    left: 50%;
    width: 5px;
    height: 5px;
    background: rgba(255, 255, 255, 0.5);
    opacity: 0;
    border-radius: 100%;
    transform: scale(1, 1) translate(-50%);
    transform-origin: 50% 50%;
  }

  .stButton > button:focus:not(:active)::after {
    animation: ripple 1s ease-out;
  }

  @keyframes ripple {
    0% {
      transform: scale(0, 0);
      opacity: 0.5;
    }
    100% {
      transform: scale(20, 20);
      opacity: 0;
    }
  }

  /* Secondary Buttons */
  .secondary-button > button {
    background: var(--surface) !important;
    color: var(--text) !important;
    border: 2px solid var(--border) !important;
    border-radius: 12px !important;
    padding: 0.875rem 2rem !important;
    font-weight: 600 !important;
    font-size: 1rem !important;
    transition: all 0.3s ease !important;
    box-shadow: 0 2px 8px var(--shadow) !important;
  }

  .secondary-button > button:hover {
    border-color: var(--accent) !important;
    background: var(--accent-light) !important;
    transform: translateY(-2px) !important;
    box-shadow: 0 4px 15px var(--shadow-strong) !important;
  }

  /* Checkboxes */
  .stCheckbox > label {
    color: var(--text) !important;
    font-weight: 500;
    padding: 0.5rem 0;
  }

  .stCheckbox [data-baseweb="checkbox"] > div {
    border-radius: 6px !important;
  }

  /* Tabs */
  .stTabs [data-baseweb="tab-list"] {
    gap: 0.5rem;
    border-bottom: 2px solid var(--border);
  }

  .stTabs [data-baseweb="tab"] {
    background: transparent;
    border: none;
    border-radius: 8px 8px 0 0;
    color: var(--text-muted);
    font-weight: 600;
    padding: 0.75rem 1.5rem;
    transition: all 0.2s ease;
    font-size: 0.95rem;
  }

  .stTabs [aria-selected="true"] {
    background: transparent;
    color: var(--accent);
    position: relative;
  }

  .stTabs [aria-selected="true"]::after {
    content: '';
    position: absolute;
    bottom: -2px;
    left: 0;
    right: 0;
    height: 3px;
    background: var(--gradient);
    border-radius: 2px 2px 0 0;
  }

  /* DataFrames */
  .stDataFrame {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 12px;
    overflow: hidden;
    box-shadow: 0 2px 8px var(--shadow);
  }

  .stDataFrame table {
    background: var(--surface) !important;
  }

  .stDataFrame th {
    background: var(--panel) !important;
    color: var(--text) !important;
    font-weight: 700;
    border-bottom: 2px solid var(--border) !important;
    text-transform: uppercase;
    font-size: 0.85rem;
    letter-spacing: 0.05em;
  }

  .stDataFrame td {
    background: var(--surface) !important;
    color: var(--text) !important;
    border-bottom: 1px solid var(--border) !important;
    font-size: 0.9rem;
  }

  .stDataFrame tr:hover td {
    background: var(--accent-light) !important;
  }

  /* Progress Bar */
  .stProgress > div > div > div {
    background: var(--gradient) !important;
    border-radius: 4px;
  }

  /* Alerts */
  .stAlert {
    background: var(--surface) !important;
    border: 1px solid var(--border) !important;
    border-radius: 12px !important;
    color: var(--text) !important;
    padding: 1rem 1.25rem !important;
    box-shadow: 0 2px 8px var(--shadow);
  }

  .stSuccess {
    border-left: 4px solid var(--success) !important;
    background: var(--success-light) !important;
  }

  .stError {
    border-left: 4px solid var(--danger) !important;
    background: var(--danger-light) !important;
  }

  .stWarning {
    border-left: 4px solid var(--warning) !important;
    background: var(--warning-light) !important;
  }

  .stInfo {
    border-left: 4px solid var(--info) !important;
    background: var(--info-light) !important;
  }

  /* Metrics */
  .stMetric {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 12px;
    padding: 1.25rem;
    box-shadow: 0 2px 8px var(--shadow);
  }

  .stMetric > div {
    color: var(--text) !important;
  }

  .stMetric label {
    color: var(--text-muted) !important;
    font-weight: 600 !important;
    font-size: 0.9rem !important;
    text-transform: uppercase;
    letter-spacing: 0.05em;
  }

  .stMetric [data-testid="stMetricValue"] {
    font-size: 2rem !important;
    font-weight: 700 !important;
    color: var(--text) !important;
    font-family: 'Source Sans Pro', sans-serif;
  }

  /* Download Buttons */
  .stDownloadButton > button {
    background: var(--surface) !important;
    color: var(--text) !important;
    border: 2px solid var(--border) !important;
    border-radius: 12px !important;
    font-weight: 600 !important;
    transition: all 0.2s ease !important;
    padding: 0.75rem 1.5rem !important;
  }

  .stDownloadButton > button:hover {
    background: var(--accent-light) !important;
    border-color: var(--accent) !important;
    transform: translateY(-2px) !important;
    box-shadow: 0 4px 15px var(--shadow-strong) !important;
  }

  /* Plotly Charts */
  .js-plotly-plot {
    background: var(--surface) !important;
    border-radius: 16px;
    overflow: hidden;
    border: 1px solid var(--border);
    box-shadow: 0 4px 20px var(--shadow);
  }

  /* Sidebar Content */
  .sidebar-title {
    font-weight: 800;
    font-size: 1.3rem;
    margin-bottom: 1rem;
    color: var(--accent);
    font-family: 'Source Sans Pro', sans-serif;
    background: var(--gradient);
    -webkit-background-clip: text;
    background-clip: text;
    color: transparent;
    padding-bottom: 0.5rem;
    border-bottom: 2px solid var(--border);
  }

  .sidebar-tip {
    color: var(--text-muted);
    font-size: 0.9rem;
    line-height: 1.5;
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 12px;
    padding: 1rem;
    margin: 0.5rem 0;
  }

  /* Expanders */
  .streamlit-expanderHeader {
    background: var(--panel) !important;
    border: 1px solid var(--border) !important;
    border-radius: 12px !important;
    color: var(--text) !important;
    font-weight: 600 !important;
  }

  .streamlit-expanderHeader:hover {
    background: var(--accent-light) !important;
    border-color: var(--accent) !important;
  }

  .streamlit-expanderContent {
    background: var(--surface) !important;
    border: 1px solid var(--border) !important;
    border-radius: 0 0 12px 12px !important;
    border-top: none !important;
  }

  /* Tooltips */
  [data-baseweb="tooltip"] {
    background: var(--surface) !important;
    color: var(--text) !important;
    border: 1px solid var(--border) !important;
    border-radius: 8px !important;
    box-shadow: 0 4px 20px var(--shadow) !important;
  }

  /* Scrollbars */
  ::-webkit-scrollbar {
    width: 10px;
    height: 10px;
  }

  ::-webkit-scrollbar-track {
    background: var(--panel);
    border-radius: 10px;
  }

  ::-webkit-scrollbar-thumb {
    background: var(--border-strong);
    border-radius: 10px;
    border: 2px solid var(--panel);
  }

  ::-webkit-scrollbar-thumb:hover {
    background: var(--accent);
  }

  /* Responsive Design */
  @media (max-width: 768px) {
    .hero {
      padding: 2rem 1.5rem;
      margin: 0.5rem 0 1.5rem 0;
    }
    
    .hero-content h1 {
      font-size: 2rem;
    }
    
    .block-container {
      padding-left: 1rem;
      padding-right: 1rem;
    }
    
    .section-title {
      font-size: 1.3rem;
      margin: 2rem 0 1rem 0;
    }

    .mode-button {
      padding: 1.5rem 1rem !important;
      height: 120px !important;
      font-size: 1rem !important;
    }
    
    .mode-icon {
      font-size: 2rem !important;
      margin-bottom: 0.75rem !important;
    }
  }

  /* Animations */
  @keyframes fadeIn {
    from { opacity: 0; transform: translateY(10px); }
    to { opacity: 1; transform: translateY(0); }
  }

  .animate-in {
    animation: fadeIn 0.5s ease-out;
  }

  /* Badges */
  .badge {
    display: inline-block;
    padding: 0.25rem 0.75rem;
    background: var(--accent-light);
    color: var(--accent);
    border-radius: 20px;
    font-size: 0.8rem;
    font-weight: 600;
    margin: 0.25rem;
  }

  /* Status Indicators */
  .status-indicator {
    display: inline-flex;
    align-items: center;
    gap: 0.5rem;
    padding: 0.25rem 0.75rem;
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 20px;
    font-size: 0.85rem;
    font-weight: 500;
  }

  .status-dot {
    width: 8px;
    height: 8px;
    border-radius: 50%;
    background: var(--success);
  }

  .status-dot.pending {
    background: var(--warning);
  }

  .status-dot.error {
    background: var(--danger);
  }

  /* Loading Spinner */
  .loading-spinner {
    display: inline-block;
    width: 20px;
    height: 20px;
    border: 3px solid var(--border);
    border-top-color: var(--accent);
    border-radius: 50%;
    animation: spin 1s linear infinite;
  }

  @keyframes spin {
    0% { transform: rotate(0deg); }
    100% { transform: rotate(360deg); }
  }

  /* Feature Cards */
  .feature-card {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 16px;
    padding: 1.5rem;
    text-align: center;
    transition: all 0.3s ease;
  }

  .feature-card:hover {
    transform: translateY(-5px);
    box-shadow: 0 12px 30px var(--shadow-strong);
    border-color: var(--accent);
  }

  .feature-icon {
    font-size: 2.5rem;
    margin-bottom: 1rem;
    color: var(--accent);
  }

  .feature-title {
    font-weight: 700;
    color: var(--text);
    margin-bottom: 0.5rem;
    font-family: 'Source Sans Pro', sans-serif;
  }

  .feature-desc {
    color: var(--text-muted);
    font-size: 0.9rem;
    line-height: 1.5;
  }
</style>
"""

# Apply light theme to entire app
st.markdown(get_light_theme_css(), unsafe_allow_html=True)

# =============================================================================
# DATA VALIDATION AND CLEANING FUNCTIONS
# =============================================================================

def validate_and_clean_count_matrix(count_matrix):
    """
    Validate and clean count matrix before analysis
    """
    # Create a copy to avoid modifying original
    cleaned_matrix = count_matrix.copy()
    
    # Check for duplicates
    duplicates = cleaned_matrix.index.duplicated().sum()
    if duplicates > 0:
        st.warning(f"Found {duplicates} duplicate gene identifiers. Making them unique...")
        
        # Make duplicates unique
        unique_index = []
        gene_counts = {}
        
        for gene in cleaned_matrix.index:
            if gene not in gene_counts:
                gene_counts[gene] = 1
                unique_index.append(gene)
            else:
                gene_counts[gene] += 1
                unique_index.append(f"{gene}_dup{gene_counts[gene]}")
        
        cleaned_matrix.index = unique_index
    
    # Check for missing values
    nan_count = cleaned_matrix.isna().sum().sum()
    if nan_count > 0:
        st.warning(f"Found {nan_count} NaN values. Replacing with zeros...")
        cleaned_matrix = cleaned_matrix.fillna(0)
    
    # Check for negative values (shouldn't exist in count data)
    negative_count = (cleaned_matrix < 0).sum().sum()
    if negative_count > 0:
        st.warning(f"Found {negative_count} negative values. Taking absolute values...")
        cleaned_matrix = cleaned_matrix.abs()
    
    # Check for all-zero rows
    zero_rows = (cleaned_matrix.sum(axis=1) == 0).sum()
    if zero_rows > 0:
        st.info(f"Found {zero_rows} genes with zero counts across all samples")
    
    return cleaned_matrix

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
    if 'degs_running' not in st.session_state:
        st.session_state.degs_running = False
    if 'current_pipeline_step' not in st.session_state:
        st.session_state.current_pipeline_step = 'selection'
    if 'degs_results' not in st.session_state:
        st.session_state.degs_results = None
    if 'analysis_results' not in st.session_state:
        st.session_state.analysis_results = {}
    if 'entrez_email_set' not in st.session_state:
        st.session_state.entrez_email_set = False
    if 'deseq2_run_completed' not in st.session_state:
        st.session_state.deseq2_run_completed = False

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
    batch_size = 1000
    
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
    
    # Check for zero variance genes
    all_samples = group1_samples + group2_samples
    variances = count_matrix[all_samples].var(axis=1)
    zero_var_genes = (variances == 0).sum()
    if zero_var_genes > 0:
        quality_issues.append(f"Found {zero_var_genes} genes with zero variance across all samples")
    
    # Check for negative values
    negative_count = (count_matrix[group1_samples + group2_samples] < 0).sum().sum()
    if negative_count > 0:
        quality_issues.append(f"Found {negative_count} negative values in the data")
    
    # Check for duplicate indices
    duplicate_genes = count_matrix.index.duplicated().sum()
    if duplicate_genes > 0:
        quality_issues.append(f"Found {duplicate_genes} duplicate gene identifiers")
    
    return quality_issues

# =============================================================================
# NORMALIZATION AND ANALYSIS FUNCTIONS
# =============================================================================

def normalize_to_cpm(count_matrix):
    """
    Normalize raw counts to CPM (Counts Per Million)
    """
    library_sizes = count_matrix.sum(axis=0)
    cpm_matrix = count_matrix.div(library_sizes, axis=1) * 1e6
    return cpm_matrix

@st.cache_resource
def compute_deseq2_results(_count_matrix, sample_group1, sample_group2, 
                          group1_name="Control", group2_name="Treatment"):
    """
    Compute DESeq2 results using cache_resource for DeseqDataSet objects
    Returns only the results DataFrame, not the DeseqDataSet object
    """
    try:
        # Prepare sample metadata
        samples = sample_group1 + sample_group2
        conditions = [group1_name] * len(sample_group1) + [group2_name] * len(sample_group2)
        
        # Create metadata DataFrame with samples as index
        metadata = pd.DataFrame({
            'condition': conditions
        }, index=samples)
        
        # Select the count data for the samples
        count_data = _count_matrix[samples]
        
        # IMPORTANT: Transpose the count matrix so that samples are rows and genes are columns
        # DESeq2 expects: rows = samples, columns = genes
        count_data_transposed = count_data.T
        
        # Create DESeqDataSet
        dds = DeseqDataSet(
            counts=count_data_transposed,
            metadata=metadata,
            design_factors='condition',
            ref_level=[group1_name]
        )
        
        # Run DESeq2
        dds.deseq2()
        
        # Get results
        stat_res = DeseqStats(dds, 
                             contrast=['condition', group2_name, group1_name])
        stat_res.summary()
        results_df = stat_res.results_df
        
        # Format results to match our expected format
        results = pd.DataFrame({
            'Gene': results_df.index,
            'logFC': results_df['log2FoldChange'],
            'p_value': results_df['pvalue'],
            'adj_p_value': results_df['padj'],
            'baseMean': results_df['baseMean'],
            'lfcSE': results_df['lfcSE'],
            'stat': results_df['stat']
        })
        
        # Handle NaN values
        results['p_value'] = results['p_value'].fillna(1.0)
        results['adj_p_value'] = results['adj_p_value'].fillna(1.0)
        
        # Add mean expression values from the original (non-transposed) data
        results['mean_group1'] = count_data[sample_group1].mean(axis=1)
        results['mean_group2'] = count_data[sample_group2].mean(axis=1)
        
        return results
        
    except ImportError:
        raise ImportError(
            "pyDESeq2 is not installed. Please install it with: "
            "pip install pydeseq2"
        )
    except Exception as e:
        raise Exception(f"DESeq2 analysis failed: {str(e)}")

def calculate_differential_expression_ttest(normalized_matrix, sample_group1, sample_group2):
    """
    Calculate differential expression using t-test on normalized data
    """
    # Apply log2 transformation for normalized data
    log_data = np.log2(normalized_matrix + 0.1)
    
    # Calculate means
    mean_group1 = log_data[sample_group1].mean(axis=1)
    mean_group2 = log_data[sample_group2].mean(axis=1)
    
    # Calculate logFC
    logFC = mean_group2 - mean_group1
    
    # Calculate p-values
    p_values = []
    for idx in range(len(log_data)):
        g1_data = log_data[sample_group1].iloc[idx].values
        g2_data = log_data[sample_group2].iloc[idx].values
        
        try:
            t_stat, p_val = ttest_ind(g1_data, g2_data, equal_var=False, nan_policy='omit')
            p_values.append(p_val if not np.isnan(p_val) else 1.0)
        except:
            p_values.append(1.0)
    
    # Create results DataFrame
    results = pd.DataFrame({
        'Gene': normalized_matrix.index,
        'logFC': logFC,
        'p_value': p_values,
        'mean_group1': np.exp2(mean_group1) - 0.1,
        'mean_group2': np.exp2(mean_group2) - 0.1
    })
    
    return results

@st.cache_data
def cached_calculate_ttest(_count_matrix, sample_group1, sample_group2, data_is_normalized):
    """
    Cached version of t-test analysis
    """
    if data_is_normalized:
        # Use t-test on already normalized data
        results = calculate_differential_expression_ttest(
            _count_matrix, sample_group1, sample_group2
        )
    else:
        # If data is not normalized, apply CPM normalization
        analysis_matrix = normalize_to_cpm(_count_matrix)
        results = calculate_differential_expression_ttest(
            analysis_matrix, sample_group1, sample_group2
        )
    
    return results

def run_differential_expression_analysis(count_matrix, sample_group1, sample_group2,
                                        data_is_normalized, analysis_method="t-test",
                                        use_caching=True):
    """
    Main function to run differential expression analysis with proper caching
    """
    if analysis_method == "DESeq2" and not data_is_normalized:
        if use_caching:
            # Use cached DESeq2 computation
            results = compute_deseq2_results(
                count_matrix, sample_group1, sample_group2
            )
        else:
            # Non-cached DESeq2
            results = compute_deseq2_results(
                count_matrix, sample_group1, sample_group2
            )
    else:
        # Use t-test (with caching if enabled)
        if use_caching:
            results = cached_calculate_ttest(
                count_matrix, sample_group1, sample_group2, data_is_normalized
            )
        else:
            if data_is_normalized:
                results = calculate_differential_expression_ttest(
                    count_matrix, sample_group1, sample_group2
                )
            else:
                normalized_matrix = normalize_to_cpm(count_matrix)
                results = calculate_differential_expression_ttest(
                    normalized_matrix, sample_group1, sample_group2
                )
    
    return results

def run_degs_analysis():
    """Run the DEGs analysis pipeline"""
    st.markdown("""
    <div class="hero">
        <div class="hero-content">
            <h1>🧬 Differential Expression Analyzer</h1>
            <p>Identify significantly expressed genes between sample groups using RNA-seq count data. 
            Upload your data and configure analysis parameters below.</p>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Sidebar content
    with st.sidebar:
        st.markdown('<div class="sidebar-title">📁 File Format Guide</div>', unsafe_allow_html=True)
        st.markdown('<div class="sidebar-tip">CSV/TSV file with genes as rows and samples as columns. First column should contain gene identifiers.</div>', unsafe_allow_html=True)
        
        st.markdown("---")
        
        st.markdown('<div class="sidebar-title">⚡ Performance Tips</div>', unsafe_allow_html=True)
        st.markdown('<div class="sidebar-tip">• Enable caching for faster re-analysis<br>• Use normalized data when available<br>• Limit samples for quicker processing</div>', unsafe_allow_html=True)
    
    # Main content in tabs
    tab1, tab2, tab3 = st.tabs(["📤 Upload Data", "⚙️ Configure", "📊 Results"])
    
    with tab1:
        st.markdown('<div class="section-title">Upload Your Data</div>', unsafe_allow_html=True)
        
        col1, col2 = st.columns([2, 1])
        
        with col1:
            uploaded_file = st.file_uploader(
                "Choose count matrix file",
                type=['csv', 'tsv'],
                help="Upload CSV or TSV file with genes as rows and samples as columns",
                label_visibility="collapsed"
            )
        
        with col2:
            data_is_normalized = st.checkbox(
                "Data is normalized",
                value=False,
                help="Check if your data is already normalized (FPKM, TPM, etc.)"
            )
        
        if uploaded_file is not None:
            try:
                if uploaded_file.name.endswith('.csv'):
                    count_matrix = pd.read_csv(uploaded_file, index_col=0)
                else:
                    count_matrix = pd.read_csv(uploaded_file, sep='\t', index_col=0)
                
                # Validate and clean the data
                count_matrix = validate_and_clean_count_matrix(count_matrix)
                
                st.success(f"✅ Data loaded: {count_matrix.shape[0]:,} genes, {count_matrix.shape[1]} samples")
                
                # Data preview - limit display rows
                with st.expander("📋 Data Preview"):
                    col_preview1, col_preview2 = st.columns(2)
                    with col_preview1:
                        display_rows = min(10, len(count_matrix))
                        st.dataframe(count_matrix.head(display_rows), width="stretch")
                    with col_preview2:
                        st.metric("Rows (Genes)", count_matrix.shape[0])
                        st.metric("Columns (Samples)", count_matrix.shape[1])
                        st.metric("Total Values", count_matrix.shape[0] * count_matrix.shape[1])
                
                # Show normalization warning if data is not normalized
                if not data_is_normalized:
                    st.warning("""
                    ⚠️ **Raw count data detected**
                    
                    For accurate differential expression analysis of raw RNA-seq counts, 
                    we recommend using the DESeq2 method which properly models count data.
                    """)
                
                st.session_state.count_matrix = count_matrix
                st.session_state.uploaded_file = uploaded_file.name
                st.session_state.data_is_normalized = data_is_normalized
                
            except Exception as e:
                st.error(f"❌ Error processing file: {str(e)}")
        else:
            st.info("👆 Upload a file to begin analysis")
    
    with tab2:
        st.markdown('<div class="section-title">Analysis Configuration</div>', unsafe_allow_html=True)
        
        # Check if count matrix exists in session state
        if 'count_matrix' in st.session_state and st.session_state.count_matrix is not None:
            count_matrix = st.session_state.count_matrix
            
            # Data quality check section
            st.subheader("🔍 Data Quality Check")
            
            duplicates = count_matrix.index.duplicated().sum()
            if duplicates > 0:
                st.error(f"❌ Found {duplicates} duplicate gene identifiers!")
            else:
                st.success("✅ All gene identifiers are unique")
            
            # Show other quality metrics
            col_qual1, col_qual2, col_qual3 = st.columns(3)
            
            with col_qual1:
                st.metric("Total Genes", len(count_matrix))
            
            with col_qual2:
                zeros = (count_matrix == 0).sum().sum()
                total = count_matrix.size
                zero_pct = (zeros / total) * 100 if total > 0 else 0
                st.metric("Zero Values", f"{zero_pct:.1f}%")
            
            with col_qual3:
                nans = count_matrix.isna().sum().sum()
                st.metric("NaN Values", nans)
            
            st.markdown("---")
            
            # Sample group configuration
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Group 1 - Control")
                group1_name = st.text_input("Group name", value="Control", key="group1_name")
                group1_samples = st.multiselect(
                    "Select samples",
                    options=count_matrix.columns.tolist(),
                    key="group1_samples"
                )
                if group1_samples:
                    st.markdown(f'<span class="badge">{len(group1_samples)} samples selected</span>', unsafe_allow_html=True)
            
            with col2:
                st.subheader("Group 2 - Treatment")
                group2_name = st.text_input("Group name", value="Treatment", key="group2_name")
                available_samples = [col for col in count_matrix.columns if col not in group1_samples]
                group2_samples = st.multiselect(
                    "Select samples",
                    options=available_samples,
                    key="group2_samples"
                )
                if group2_samples:
                    st.markdown(f'<span class="badge">{len(group2_samples)} samples selected</span>', unsafe_allow_html=True)
            
            # Analysis method selection
            st.subheader("Analysis Method")
            
            if st.session_state.data_is_normalized:
                st.info("📊 **Using t-test** - Your data is already normalized")
                analysis_method = "t-test"
            else:
                # Let user choose method for raw data
                analysis_method = st.radio(
                    "Select analysis method:",
                    options=["DESeq2 (recommended for raw counts)", "t-test (with CPM normalization)"],
                    index=0,
                    help="""
                    **DESeq2**: Full Negative Binomial model for raw counts (most accurate)
                    **t-test**: Simple t-test on CPM-normalized data (faster but less accurate)
                    """
                )
                
                # Extract method name
                if "DESeq2" in analysis_method:
                    analysis_method = "DESeq2"
                    st.success("✅ **DESeq2 selected** - Using Negative Binomial model for raw counts")
                else:
                    analysis_method = "t-test"
                    st.info("📊 **t-test selected** - Using CPM normalization + t-test")
            
            # Statistical parameters
            st.subheader("Statistical Parameters")
            col_param1, col_param2 = st.columns(2)
            
            with col_param1:
                logFC_threshold = st.slider(
                    "logFC Threshold",
                    min_value=0.0,
                    max_value=5.0,
                    value=1.0,
                    step=0.1,
                    help="Minimum absolute log2 fold change for significance"
                )
            
            with col_param2:
                p_value_threshold = st.select_slider(
                    "P-value Threshold",
                    options=[0.001, 0.01, 0.05, 0.1],
                    value=0.05,
                    help="Maximum p-value for significance"
                )
            
            # Performance settings
            st.subheader("Performance Settings")
            use_caching = st.checkbox("Enable caching", value=True)
            
            # Run analysis button
            st.markdown("---")
            if st.button("🚀 Run Differential Expression Analysis", type="primary", width="stretch"):
                if not group1_samples or not group2_samples:
                    st.warning("⚠️ Please select samples for both groups")
                else:
                    # Check for overlap
                    overlap = set(group1_samples) & set(group2_samples)
                    if overlap:
                        st.error(f"Sample overlap detected: {list(overlap)}")
                    else:
                        # Store parameters
                        st.session_state.analysis_params = {
                            'group1_name': group1_name,
                            'group2_name': group2_name,
                            'group1_samples': group1_samples,
                            'group2_samples': group2_samples,
                            'logFC_threshold': logFC_threshold,
                            'p_value_threshold': p_value_threshold,
                            'use_caching': use_caching,
                            'analysis_method': analysis_method
                        }
                        st.session_state.degs_running = True
                        st.rerun()
        else:
            st.info("📤 Please upload data in the 'Upload Data' tab first")
    
    with tab3:
        st.markdown('<div class="section-title">Analysis Results</div>', unsafe_allow_html=True)
        
        if ('analysis_params' in st.session_state and 
            'count_matrix' in st.session_state and 
            st.session_state.count_matrix is not None):
            
            params = st.session_state.analysis_params
            count_matrix = st.session_state.count_matrix
            
            # Run analysis if not already completed
            if st.session_state.degs_running and not st.session_state.degs_completed:
                # Setup progress
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                with st.spinner("🔬 Analyzing differential expression..."):
                    # Update progress
                    status_text.text("Initializing analysis...")
                    progress_bar.progress(10)
                    
                    try:
                        # Get analysis method
                        analysis_method = params['analysis_method']
                        
                        # Perform analysis
                        status_text.text(f"Running {analysis_method} analysis...")
                        progress_bar.progress(30)
                        
                        results = run_differential_expression_analysis(
                            count_matrix,
                            params['group1_samples'],
                            params['group2_samples'],
                            st.session_state.data_is_normalized,
                            analysis_method,
                            params['use_caching']
                        )
                        
                        status_text.text("Filtering significant genes...")
                        progress_bar.progress(70)
                        
                        # Apply filtering outside of cached function
                        upregulated, downregulated = filter_and_sort_degs(
                            results,
                            params['logFC_threshold'],
                            params['p_value_threshold']
                        )
                        
                        progress_bar.progress(100)
                        status_text.text("✅ Analysis complete!")
                        time.sleep(0.5)
                        
                        # Clear progress
                        progress_bar.empty()
                        status_text.empty()
                        
                        # Store results in session state
                        st.session_state.degs_results = {
                            'results': results,
                            'upregulated': upregulated,
                            'downregulated': downregulated,
                            'params': params,
                            'data_is_normalized': st.session_state.data_is_normalized,
                            'analysis_method': analysis_method
                        }
                        st.session_state.degs_completed = True
                        st.session_state.degs_running = False
                        st.session_state.deseq2_run_completed = True
                        
                    except ImportError as e:
                        st.error(f"""
                        ❌ **pyDESeq2 not installed**
                        
                        To use DESeq2 analysis, please install pyDESeq2:
                        
                        ```bash
                        pip install pydeseq2
                        ```
                        
                        Then restart the application.
                        """)
                        st.session_state.degs_running = False
                        return
                    except Exception as e:
                        st.error(f"❌ Analysis failed: {str(e)}")
                        st.session_state.degs_running = False
                        return
            
            # Display results if available
            if st.session_state.degs_completed and st.session_state.degs_results:
                results_data = st.session_state.degs_results
                
                # Summary metrics
                st.markdown("### 📈 Summary Statistics")
                col1, col2, col3, col4 = st.columns(4)
                
                with col1:
                    st.metric("Total Genes", len(results_data['results']))
                with col2:
                    total_degs = len(results_data['upregulated']) + len(results_data['downregulated'])
                    st.metric("Significant DEGs", total_degs)
                with col3:
                    st.metric("Upregulated", len(results_data['upregulated']))
                with col4:
                    st.metric("Downregulated", len(results_data['downregulated']))
                
                # Show analysis method info
                if results_data['analysis_method'] == "DESeq2":
                    st.success(f"✅ **Analysis Method:** DESeq2 (Negative Binomial model)")
                    if 'adj_p_value' in results_data['results'].columns:
                        st.info(f"🔬 **Multiple testing correction:** Benjamini-Hochberg FDR applied")
                else:
                    st.info(f"📊 **Analysis Method:** {results_data['analysis_method']}")
                
                # Gene tables - limit display rows
                st.markdown("### 🧬 Top Differentially Expressed Genes")
                
                tab_up, tab_down = st.tabs(["⬆️ Upregulated", "⬇️ Downregulated"])
                
                with tab_up:
                    if not results_data['upregulated'].empty:
                        # Include adjusted p-value if available
                        display_cols = ['Gene', 'logFC', 'p_value']
                        if 'adj_p_value' in results_data['upregulated'].columns:
                            display_cols.append('adj_p_value')
                        
                        # Limit display rows
                        display_rows = min(Config.MAX_DISPLAY_ROWS, len(results_data['upregulated']))
                        display_up = results_data['upregulated'].head(display_rows)[display_cols].round(4)
                        st.dataframe(display_up, width="stretch")
                        
                        if len(results_data['upregulated']) > Config.MAX_DISPLAY_ROWS:
                            st.info(f"Showing first {Config.MAX_DISPLAY_ROWS} of {len(results_data['upregulated'])} upregulated genes. Download for full results.")
                    else:
                        st.info("No upregulated genes found")
                
                with tab_down:
                    if not results_data['downregulated'].empty:
                        # Include adjusted p-value if available
                        display_cols = ['Gene', 'logFC', 'p_value']
                        if 'adj_p_value' in results_data['downregulated'].columns:
                            display_cols.append('adj_p_value')
                        
                        # Limit display rows
                        display_rows = min(Config.MAX_DISPLAY_ROWS, len(results_data['downregulated']))
                        display_down = results_data['downregulated'].head(display_rows)[display_cols].round(4)
                        st.dataframe(display_down, width="stretch")
                        
                        if len(results_data['downregulated']) > Config.MAX_DISPLAY_ROWS:
                            st.info(f"Showing first {Config.MAX_DISPLAY_ROWS} of {len(results_data['downregulated'])} downregulated genes. Download for full results.")
                    else:
                        st.info("No downregulated genes found")
                
                # Volcano plot with protection against log(0) errors
                st.markdown("### 🌋 Volcano Plot")
                try:
                    plot_data = results_data['results'].copy()
                    
                    # Use adjusted p-value if available, otherwise use raw p-value
                    if 'adj_p_value' in plot_data.columns:
                        pval_col = 'adj_p_value'
                        pval_label = "Adjusted p-value"
                    else:
                        pval_col = 'p_value'
                        pval_label = "p-value"
                    
                    # Clip p-values to avoid log(0) errors
                    min_pval = 1e-300  # Smallest safe value for log10
                    plot_data['pval_clipped'] = plot_data[pval_col].clip(lower=min_pval)
                    plot_data['-log10(p_value)'] = -np.log10(plot_data['pval_clipped'])
                    
                    plot_data['Significant'] = (np.abs(plot_data['logFC']) > params['logFC_threshold']) & (plot_data[pval_col] < params['p_value_threshold'])
                    
                    fig = px.scatter(
                        plot_data,
                        x='logFC',
                        y='-log10(p_value)',
                        color='Significant',
                        hover_data=['Gene'],
                        title=f"{params['group1_name']} vs {params['group2_name']}",
                        labels={'logFC': 'log2 Fold Change', '-log10(p_value)': f'-log10({pval_label})'},
                        color_discrete_map={True: '#ef4444', False: '#94a3b8'},
                        opacity=0.7
                    )
                    
                    # Add threshold lines
                    fig.add_vline(x=params['logFC_threshold'], line_dash="dash", line_color="#10b981")
                    fig.add_vline(x=-params['logFC_threshold'], line_dash="dash", line_color="#10b981")
                    fig.add_hline(y=-np.log10(params['p_value_threshold']), line_dash="dash", line_color="#10b981")
                    
                    fig.update_layout(
                        plot_bgcolor='white',
                        paper_bgcolor='white',
                        font_color='#1e293b',
                        title_font_color='#3b82f6'
                    )
                    
                    st.plotly_chart(fig, width="stretch")
                except Exception as e:
                    st.warning(f"Could not generate volcano plot: {e}")
                
                # Download section
                st.markdown("### 📥 Download Results")
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                
                col_dl1, col_dl2 = st.columns(2)
                
                with col_dl1:
                    # All DEGs
                    all_degs = pd.concat([results_data['upregulated'], results_data['downregulated']])
                    if not all_degs.empty:
                        csv_all = all_degs.to_csv(index=False)
                        st.download_button(
                            label="Download Significant DEGs",
                            data=csv_all,
                            file_name=f"significant_DEGs_{timestamp}.csv",
                            mime="text/csv",
                            width="stretch"
                        )
                
                with col_dl2:
                    # Complete results
                    csv_complete = results_data['results'].to_csv(index=False)
                    st.download_button(
                        label="Download Complete Results",
                        data=csv_complete,
                        file_name=f"complete_DE_results_{timestamp}.csv",
                        mime="text/csv",
                        width="stretch"
                    )
                
                # Continue to pathway analysis
                st.markdown("---")
                st.markdown("### 🔄 Continue Analysis")
                
                if st.button("🚀 Continue to Pathway & Drug Analysis", type="primary", width="stretch"):
                    # Get top genes for pathway analysis
                    top_up = results_data['upregulated'].head(10)['Gene'].tolist() if not results_data['upregulated'].empty else []
                    top_down = results_data['downregulated'].head(10)['Gene'].tolist() if not results_data['downregulated'].empty else []
                    pathway_genes = top_up + top_down
                    
                    st.session_state.pathway_genes = pathway_genes
                    st.session_state.analysis_mode = 'pathway_only'
                    st.session_state.current_pipeline_step = 'pathway'
                    st.rerun()
        
        elif 'count_matrix' not in st.session_state:
            st.info("📤 Please upload data in the 'Upload Data' tab first")
        else:
            st.info("⚙️ Configure and run the analysis in the 'Configure' tab")

# =============================================================================
# PATHWAY ANALYSIS FUNCTIONS (Second Pipeline)
# =============================================================================

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
                time.sleep(2 ** attempt)
                
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

# Retry decorator for API calls
def with_retry(func, max_retries=3, delay=1):
    """Decorator pattern for retrying API calls"""
    def wrapper(*args, **kwargs):
        for attempt in range(max_retries):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                if attempt == max_retries - 1:
                    raise e
                time.sleep(delay * (2 ** attempt))
        return None
    return wrapper

# Gene symbol validation
def validate_gene_symbol(gene: str) -> bool:
    """Validate gene symbol format"""
    if not gene or len(gene) < 2:
        return False
    # Remove version numbers if present (e.g., TP53.1 -> TP53)
    gene = re.sub(r'\..*', '', gene)
    # Check for valid characters
    if not re.match(r'^[A-Z0-9\-_]+$', gene.upper()):
        return False
    return True

def filter_valid_genes(gene_list: List[str]) -> Tuple[List[str], List[str]]:
    """Filter and validate gene symbols"""
    valid_genes = []
    invalid_genes = []
    
    for gene in gene_list:
        if validate_gene_symbol(gene):
            valid_genes.append(gene.upper())
        else:
            invalid_genes.append(gene)
    
    return valid_genes, invalid_genes

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
                if v and v not in seen and len(v) > 1:
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
            
            for key in ("gene.symbol", "symbol", "gene_symbol", "gene", "geneid"):
                if key in lower_map:
                    target_col = lower_map[key]
                    break
                    
            if target_col is None:
                target_col = df.columns[0]
                
            raw_genes = _clean_series_to_genes(df[target_col])
            valid_genes, invalid_genes = filter_valid_genes(raw_genes)
            if invalid_genes:
                logger.warning(f"Filtered out {len(invalid_genes)} invalid gene symbols")
            return valid_genes
            
    except Exception as e:
        logger.error(f"Error reading structured file: {e}")

    # Fallback to text parsing
    try:
        uploaded_file.seek(0)
        raw = uploaded_file.read()
        text = raw.decode("utf-8", errors="ignore") if isinstance(raw, (bytes, bytearray)) else str(raw)
        raw_genes = _clean_series_to_genes(pd.Series([text]))
        valid_genes, invalid_genes = filter_valid_genes(raw_genes)
        if invalid_genes:
            logger.warning(f"Filtered out {len(invalid_genes)} invalid gene symbols")
        return valid_genes
    except Exception as e:
        logger.error(f"Error reading as text: {e}")
        return []

# Enhanced Caching Functions
@with_retry
@st.cache_data(ttl=3600, show_spinner=False)
def kegg_get(path: str) -> str:
    """Enhanced KEGG API call with better error handling"""
    try:
        response = kegg_session.get(f"https://rest.kegg.jp{path}")
        return response.text if response else ""
    except Exception as e:
        logger.error(f"KEGG API error for {path}: {e}")
        return ""

@with_retry
@st.cache_data(ttl=3600, show_spinner=False)
def ncbi_esearch_gene_ids(gene_symbol: str, organism_entrez: str) -> list[str]:
    """Enhanced NCBI gene search"""
    try:
        handle = Entrez.esearch(
            db="gene",
            term=f"{gene_symbol}[Gene] AND {organism_entrez}[Organism]",
            retmode="xml",
            retmax=5
        )
        record = Entrez.read(handle)
        handle.close()
        return record.get("IdList", [])
    except Exception as e:
        logger.error(f"NCBI search error for {gene_symbol}: {e}")
        return []

@with_retry
@st.cache_data(ttl=3600, show_spinner=False)
def ncbi_esummary_description(gene_id: str) -> str:
    """Fixed: Enhanced NCBI gene description fetch with None check"""
    try:
        handle = Entrez.esummary(db="gene", id=gene_id, retmode="xml")
        raw_xml = handle.read()
        handle.close()
        root = ET.fromstring(raw_xml)
        docsum = root.find(".//DocumentSummary")
        if docsum is not None:
            description = docsum.findtext("Description", default="")
            return description.strip() if description else "Description unavailable"
        else:
            logger.warning(f"No DocumentSummary found for gene_id: {gene_id}")
            return "Description unavailable"
    except Exception as e:
        logger.error(f"NCBI summary error for {gene_id}: {e}")
        return "Description unavailable"

@with_retry
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

@with_retry
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

@with_retry
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
    
    if 'errors' in data and data['errors']:
        logger.error(f"OpenTargets API errors: {data['errors']}")
        return False
    
    return True

@with_retry
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

@with_retry
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

@with_retry
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

@with_retry
@st.cache_data(ttl=3600, show_spinner=False)
def ot_drugs_for_target(ensembl_id: str, size: int = 50) -> pd.DataFrame:
    """Fixed: Enhanced drug suggestions fetch with correct phase parsing"""
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
        target_data = ((data or {}).get("data", {})).get("target", {})
        if not target_data:
            return pd.DataFrame()
            
        known_drugs = target_data.get("knownDrugs", {})
        rows = known_drugs.get("rows", [])
        
        if not rows:
            return pd.DataFrame()

        out = []
        for r in rows:
            drug_obj = r.get("drug") or {}
            disease_obj = r.get("disease") or {}
            
            therapeutic_areas = disease_obj.get("therapeuticAreas", [])
            ta_names = "; ".join([ta.get("name", "") for ta in therapeutic_areas if ta.get("name")])
            
            # Parse phase value correctly
            phase = r.get("phase")
            try:
                # Handle "4" or "Phase 4" or "IV" formats
                if isinstance(phase, str):
                    if phase.isdigit():
                        phase_numeric = int(phase)
                    elif "Phase" in phase:
                        phase_numeric = int(phase.replace("Phase", "").strip())
                    elif phase.upper() == "IV":
                        phase_numeric = 4
                    elif phase.upper() == "III":
                        phase_numeric = 3
                    elif phase.upper() == "II":
                        phase_numeric = 2
                    elif phase.upper() == "I":
                        phase_numeric = 1
                    else:
                        phase_numeric = 0
                elif isinstance(phase, (int, float)):
                    phase_numeric = int(phase)
                else:
                    phase_numeric = 0
            except:
                phase_numeric = 0
            
            # Parse max_phase value
            max_phase = drug_obj.get("maximumClinicalTrialPhase")
            try:
                if isinstance(max_phase, str):
                    if max_phase.isdigit():
                        max_phase_numeric = int(max_phase)
                    else:
                        max_phase_numeric = phase_numeric  # Fallback to current phase
                elif isinstance(max_phase, (int, float)):
                    max_phase_numeric = int(max_phase)
                else:
                    max_phase_numeric = phase_numeric  # Fallback to current phase
            except:
                max_phase_numeric = phase_numeric  # Fallback to current phase
            
            out.append({
                "target": ensembl_id,
                "drug_id": drug_obj.get("id"),
                "drug_name": drug_obj.get("name"),
                "drug_type": drug_obj.get("drugType"),
                "phase": str(phase) if phase else "Unknown",
                "status": r.get("status"),
                "max_phase": str(max_phase) if max_phase else "Unknown",
                "moa": r.get("mechanismOfAction"),
                "disease_name": disease_obj.get("name"),
                "therapeutic_areas": ta_names,
                "phase_numeric": phase_numeric,
                "max_phase_numeric": max_phase_numeric,
            })
        
        return pd.DataFrame(out)
    
    except Exception as e:
        logger.error(f"Error fetching drugs for target {ensembl_id}: {e}")
        return pd.DataFrame()

# Enhanced Core Functions
@st.cache_data(ttl=3600, show_spinner=False)
def fetch_gene_metadata_and_kegg(gene_list: list[str], organism_entrez: str,
                                 kegg_org_prefix: str):
    """Enhanced metadata fetch with robust error handling"""
    results = []
    pathway_to_genes = defaultdict(set)
    
    total_genes = len(gene_list)
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    for i, gene in enumerate(gene_list, start=1):
        try:
            # Update progress
            progress = min(90, int((i / total_genes) * 90))
            progress_bar.progress(progress)
            status_text.text(f"Processing {i}/{total_genes}: {gene}")
            
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
            
            # Get description with error handling
            description = "Description unavailable"
            try:
                description = ncbi_esummary_description(gene_id)
            except Exception as desc_error:
                logger.warning(f"Could not get description for {gene}: {desc_error}")
            
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
            pids = []
            try:
                pids = kegg_gene_pathways(kegg_id)
            except Exception as pathway_error:
                logger.warning(f"Could not get pathways for {gene}: {pathway_error}")
                
            pathway_pairs = []
            
            for pid in pids:
                try:
                    name = kegg_pathway_name(pid) or "Unknown"
                    pathway_pairs.append(f"{pid.replace('path:', '')} - {name}")
                    pathway_to_genes[pid].add(gene)
                except Exception as name_error:
                    logger.warning(f"Could not get pathway name for {pid}: {name_error}")
                    pathway_pairs.append(f"{pid.replace('path:', '')} - Unknown")
                
            pathways_str = "; ".join(pathway_pairs) if pathway_pairs else None
            
            results.append({
                "Gene": gene, 
                "NCBI_ID": gene_id, 
                "Description": description, 
                "KEGG_Pathways": pathways_str,
                "Status": f"Found {len(pids)} pathways" if pids else "No pathways"
            })
            
            # Reduced sleep time
            time.sleep(0.05)
            
        except Exception as e:
            logger.error(f"Error processing {gene}: {e}")
            results.append({
                "Gene": gene, 
                "NCBI_ID": None, 
                "Description": f"Error: {str(e)[:100]}", 
                "KEGG_Pathways": None,
                "Status": "Error"
            })
    
    # Final progress update
    progress_bar.progress(100)
    status_text.text("Metadata fetch complete!")
    time.sleep(0.5)
    progress_bar.empty()
    status_text.empty()
    
    return pd.DataFrame(results), dict(pathway_to_genes)

@st.cache_data(ttl=3600, show_spinner=False)
def compute_enrichment_counts_only(pathway_to_genes: dict) -> pd.DataFrame:
    """Enhanced enrichment computation with better sorting"""
    if not pathway_to_genes:
        return pd.DataFrame(columns=["Pathway_ID", "Pathway_Name", "Count", "Gene_List"])
        
    rows = []
    for pid, genes in pathway_to_genes.items():
        try:
            pathway_name = kegg_pathway_name(pid) or "Unknown pathway"
            gene_list = sorted(list(genes))
            
            rows.append({
                "Pathway_ID": pid.replace("path:", ""),
                "Pathway_Name": pathway_name,
                "Count": len(genes),
                "Gene_List": ", ".join(gene_list[:10]) + ("..." if len(gene_list) > 10 else "")
            })
        except Exception as e:
            logger.error(f"Error processing pathway {pid}: {e}")
            continue
    
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["Count", "Pathway_Name"], ascending=[False, True]).reset_index(drop=True)
    
    return df

@st.cache_data(ttl=3600, show_spinner=False)
def build_gene_to_ot_target_map(genes: list[str], species: str = "Homo sapiens") -> dict:
    """Enhanced gene to target mapping with progress tracking"""
    g2t = {}
    total_genes = len(genes)
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    for i, g in enumerate(genes, start=1):
        progress = int((i / total_genes) * 100)
        progress_bar.progress(progress)
        status_text.text(f"Mapping genes to targets: {i}/{total_genes}")
        
        try:
            hit = ot_target_from_symbol(g, species)
            if hit:
                g2t[g] = hit
                logger.info(f"Mapped {g} to {hit.get('id')}")
            else:
                logger.warning(f"No target found for {g}")
            time.sleep(0.05)  # Reduced sleep time
        except Exception as e:
            logger.error(f"Error mapping {g}: {e}")
    
    progress_bar.progress(100)
    status_text.text("Gene mapping complete!")
    time.sleep(0.5)
    progress_bar.empty()
    status_text.empty()
    
    return g2t

@st.cache_data(ttl=3600, show_spinner=False)
def collect_disease_links(gene_to_target: dict) -> pd.DataFrame:
    """Enhanced disease collection with better error handling"""
    frames = []
    total_targets = len(gene_to_target)
    
    if total_targets == 0:
        return pd.DataFrame(columns=["gene", "gene_symbol", "target", "disease_id", "disease_name", "association_score", "therapeutic_areas"])
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    for i, (g, tgt) in enumerate(gene_to_target.items(), start=1):
        progress = int((i / total_targets) * 100)
        progress_bar.progress(progress)
        status_text.text(f"Collecting disease associations: {i}/{total_targets}")
        
        tid = tgt.get("id")
        if not tid:
            continue
            
        try:
            df = ot_diseases_for_target(tid)
            if not df.empty:
                df.insert(0, "gene", g)
                df.insert(1, "gene_symbol", tgt.get("approvedSymbol", g))
                frames.append(df)
            time.sleep(0.05)  # Reduced sleep time
        except Exception as e:
            logger.error(f"Error collecting diseases for {g}: {e}")
    
    progress_bar.progress(100)
    status_text.text("Disease collection complete!")
    time.sleep(0.5)
    progress_bar.empty()
    status_text.empty()
    
    if frames:
        combined = pd.concat(frames, ignore_index=True)
        combined = combined.drop_duplicates(subset=["gene", "disease_id"])
        combined = combined.sort_values("association_score", ascending=False)
        return combined
    
    return pd.DataFrame(columns=["gene", "gene_symbol", "target", "disease_id", "disease_name", "association_score", "therapeutic_areas"])

@st.cache_data(ttl=3600, show_spinner=False)
def collect_drug_suggestions(gene_to_target: dict) -> pd.DataFrame:
    """Enhanced drug collection with better filtering"""
    frames = []
    total_targets = len(gene_to_target)
    
    if total_targets == 0:
        return pd.DataFrame()
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    for i, (g, tgt) in enumerate(gene_to_target.items(), start=1):
        progress = int((i / total_targets) * 100)
        progress_bar.progress(progress)
        status_text.text(f"Collecting drug suggestions: {i}/{total_targets}")
        
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
            time.sleep(0.05)  # Reduced sleep time
        except Exception as e:
            logger.error(f"Error collecting drugs for {g}: {e}")
    
    progress_bar.progress(100)
    status_text.text("Drug collection complete!")
    time.sleep(0.5)
    progress_bar.empty()
    status_text.empty()
    
    if frames:
        combined = pd.concat(frames, ignore_index=True)
        combined = combined.drop_duplicates(subset=["gene", "drug_id", "disease_name"])
        logger.info(f"Total drugs collected: {len(combined)}")
        return combined
    
    logger.warning("No drugs collected from any target")
    return pd.DataFrame()

# Network Visualization Functions with size limits
def create_gene_disease_drug_network(df_diseases, df_drugs):
    """Create a gene-disease-drug network visualization with distinct colors"""
    
    # Check node count limits
    if len(df_diseases) > Config.MAX_NETWORK_NODES or len(df_drugs) > Config.MAX_NETWORK_NODES:
        logger.warning(f"Network too large for visualization. Skipping.")
        return None
    
    # Create network graph
    G = nx.Graph()
    
    # Track unique nodes
    gene_nodes = set()
    disease_nodes = set()
    drug_nodes = set()
    
    # Add gene-disease edges from disease associations
    if not df_diseases.empty:
        # Get top diseases by gene count
        top_diseases = (
            df_diseases.groupby('disease_name')
            .agg({'gene': lambda x: len(set(x))})
            .sort_values('gene', ascending=False)
            .head(10)  # Top 10 diseases
            .index.tolist()
        )
        
        for _, row in df_diseases.iterrows():
            if row['disease_name'] in top_diseases:
                gene_id = f"Gene: {row['gene']}"
                disease_id = f"Disease: {row['disease_name']}"
                
                # Add nodes
                G.add_node(gene_id, type='gene', size=10)
                G.add_node(disease_id, type='disease', size=15)
                
                # Add edge with weight based on association score
                weight = float(row['association_score']) if row['association_score'] else 0.5
                G.add_edge(gene_id, disease_id, weight=weight, type='gene-disease')
                
                gene_nodes.add(gene_id)
                disease_nodes.add(disease_id)
    
    # Add drug-disease edges from drug data
    if not df_drugs.empty and not df_diseases.empty:
        # Get top drugs
        top_drugs = df_drugs['drug_name'].value_counts().head(10).index.tolist()
        
        for _, row in df_drugs.iterrows():
            if row['drug_name'] in top_drugs and row['disease_name']:
                drug_id = f"Drug: {row['drug_name']}"
                disease_id = f"Disease: {row['disease_name']}"
                
                # Add drug node
                G.add_node(drug_id, type='drug', size=12)
                G.add_edge(drug_id, disease_id, weight=1.0, type='drug-disease')
                
                drug_nodes.add(drug_id)
    
    # Check if we have enough nodes
    if len(G.nodes()) < 3 or len(G.nodes()) > Config.MAX_NETWORK_NODES:
        return None
    
    # Create positions using spring layout
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
    
    # Prepare edge traces
    edge_traces = []
    
    # Gene-Disease edges (blue)
    gene_disease_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('type') == 'gene-disease']
    if gene_disease_edges:
        edge_x_gd = []
        edge_y_gd = []
        for u, v in gene_disease_edges:
            x0, y0 = pos[u]
            x1, y1 = pos[v]
            edge_x_gd.extend([x0, x1, None])
            edge_y_gd.extend([y0, y1, None])
        
        edge_trace_gd = go.Scatter(
            x=edge_x_gd, y=edge_y_gd,
            line=dict(width=1.5, color='rgba(59, 130, 246, 0.6)'),  # Blue
            hoverinfo='none',
            mode='lines',
            name='Gene-Disease'
        )
        edge_traces.append(edge_trace_gd)
    
    # Drug-Disease edges (green)
    drug_disease_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('type') == 'drug-disease']
    if drug_disease_edges:
        edge_x_dd = []
        edge_y_dd = []
        for u, v in drug_disease_edges:
            x0, y0 = pos[u]
            x1, y1 = pos[v]
            edge_x_dd.extend([x0, x1, None])
            edge_y_dd.extend([y0, y1, None])
        
        edge_trace_dd = go.Scatter(
            x=edge_x_dd, y=edge_y_dd,
            line=dict(width=1.5, color='rgba(16, 185, 129, 0.6)'),  # Green
            hoverinfo='none',
            mode='lines',
            name='Drug-Disease'
        )
        edge_traces.append(edge_trace_dd)
    
    # Prepare node traces
    node_traces = []
    
    # Gene nodes (blue)
    gene_node_x = [pos[node][0] for node in gene_nodes if node in pos]
    gene_node_y = [pos[node][1] for node in gene_nodes if node in pos]
    gene_node_text = [node.replace('Gene: ', '') for node in gene_nodes if node in pos]
    
    if gene_node_x:
        gene_trace = go.Scatter(
            x=gene_node_x, y=gene_node_y,
            mode='markers+text',
            hoverinfo='text',
            text=gene_node_text,
            textposition="top center",
            hovertext=[f"Gene: {text}" for text in gene_node_text],
            marker=dict(
                size=15,
                color='#3b82f6',  # Blue
                line=dict(width=2, color='white')
            ),
            name='Genes'
        )
        node_traces.append(gene_trace)
    
    # Disease nodes (red)
    disease_node_x = [pos[node][0] for node in disease_nodes if node in pos]
    disease_node_y = [pos[node][1] for node in disease_nodes if node in pos]
    disease_node_text = [node.replace('Disease: ', '') for node in disease_nodes if node in pos]
    
    if disease_node_x:
        disease_trace = go.Scatter(
            x=disease_node_x, y=disease_node_y,
            mode='markers+text',
            hoverinfo='text',
            text=disease_node_text,
            textposition="top center",
            hovertext=[f"Disease: {text}" for text in disease_node_text],
            marker=dict(
                size=20,
                color='#ef4444',  # Red
                line=dict(width=2, color='white'),
                symbol='diamond'
            ),
            name='Diseases'
        )
        node_traces.append(disease_trace)
    
    # Drug nodes (green)
    drug_node_x = [pos[node][0] for node in drug_nodes if node in pos]
    drug_node_y = [pos[node][1] for node in drug_nodes if node in pos]
    drug_node_text = [node.replace('Drug: ', '') for node in drug_nodes if node in pos]
    
    if drug_node_x:
        drug_trace = go.Scatter(
            x=drug_node_x, y=drug_node_y,
            mode='markers+text',
            hoverinfo='text',
            text=drug_node_text,
            textposition="top center",
            hovertext=[f"Drug: {text}" for text in drug_node_text],
            marker=dict(
                size=18,
                color='#10b981',  # Green
                line=dict(width=2, color='white'),
                symbol='square'
            ),
            name='Drugs'
        )
        node_traces.append(drug_trace)
    
    # Create figure
    fig = go.Figure(data=edge_traces + node_traces)
    
    fig.update_layout(
        title=dict(
            text='Gene-Disease-Drug Network',
            x=0.5,
            font=dict(size=20, color='#3b82f6')
        ),
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
            bgcolor='rgba(255, 255, 255, 0.9)',
            bordercolor='#e2e8f0',
            borderwidth=1
        ),
        hovermode='closest',
        margin=dict(b=20, l=5, r=5, t=50),
        annotations=[
            dict(
                text="Blue: Genes | Red: Diseases | Green: Drugs",
                showarrow=False,
                xref="paper", yref="paper",
                x=0.5, y=-0.05,
                xanchor='center',
                font=dict(color='#64748b', size=12)
            )
        ],
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=600
    )
    
    return fig

def create_gene_pathway_interaction_network(pathway_to_genes, kegg_org_prefix):
    """Create a gene-pathway interaction network visualization"""
    
    if not pathway_to_genes:
        return None
    
    # Get top pathways by gene count
    top_pathways = sorted(
        pathway_to_genes.items(), 
        key=lambda x: len(x[1]), 
        reverse=True
    )[:8]  # Top 8 pathways
    
    # Check total node count
    total_genes = sum(len(genes) for _, genes in top_pathways[:15])  # Limit genes per pathway
    if total_genes > Config.MAX_NETWORK_NODES:
        logger.warning(f"Network too large ({total_genes} nodes). Skipping.")
        return None
    
    # Create network graph
    G = nx.Graph()
    
    # Track nodes
    gene_nodes = set()
    pathway_nodes = set()
    
    # Add nodes and edges
    for pathway_id, genes in top_pathways:
        pathway_name = kegg_pathway_name(pathway_id) or pathway_id.replace("path:", "")
        # Shorten pathway name for display
        if len(pathway_name) > 30:
            pathway_name = pathway_name[:27] + "..."
        
        pathway_node = f"Pathway: {pathway_name}"
        
        # Add pathway node
        G.add_node(pathway_node, type='pathway', size=25)
        pathway_nodes.add(pathway_node)
        
        # Add gene nodes and edges
        for gene in list(genes)[:15]:  # Limit to 15 genes per pathway
            gene_node = f"Gene: {gene}"
            G.add_node(gene_node, type='gene', size=15)
            G.add_edge(gene_node, pathway_node, weight=1.0)
            
            gene_nodes.add(gene_node)
    
    # Check if we have enough nodes
    if len(G.nodes()) < 3 or len(G.nodes()) > Config.MAX_NETWORK_NODES:
        return None
    
    # Create positions using spring layout
    pos = nx.spring_layout(G, k=1.5, iterations=100, seed=42)
    
    # Prepare edge traces
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1.5, color='rgba(148, 163, 184, 0.6)'),  # Gray
        hoverinfo='none',
        mode='lines'
    )
    
    # Prepare node traces
    node_traces = []
    
    # Gene nodes (blue)
    gene_node_list = [node for node in gene_nodes if node in pos]
    if gene_node_list:
        gene_node_x = [pos[node][0] for node in gene_node_list]
        gene_node_y = [pos[node][1] for node in gene_node_list]
        gene_node_text = [node.replace('Gene: ', '') for node in gene_node_list]
        
        gene_trace = go.Scatter(
            x=gene_node_x, y=gene_node_y,
            mode='markers+text',
            hoverinfo='text',
            text=gene_node_text,
            textposition="middle center",
            hovertext=[f"Gene: {text}" for text in gene_node_text],
            marker=dict(
                size=15,
                color='#3b82f6',  # Blue
                line=dict(width=2, color='white'),
                symbol='circle'
            ),
            name='Genes'
        )
        node_traces.append(gene_trace)
    
    # Pathway nodes (green)
    pathway_node_list = [node for node in pathway_nodes if node in pos]
    if pathway_node_list:
        pathway_node_x = [pos[node][0] for node in pathway_node_list]
        pathway_node_y = [pos[node][1] for node in pathway_node_list]
        pathway_node_text = [node.replace('Pathway: ', '') for node in pathway_node_list]
        
        pathway_trace = go.Scatter(
            x=pathway_node_x, y=pathway_node_y,
            mode='markers+text',
            hoverinfo='text',
            text=pathway_node_text,
            textposition="middle center",
            hovertext=[f"Pathway: {text}" for text in pathway_node_text],
            marker=dict(
                size=25,
                color='#10b981',  # Green
                line=dict(width=2, color='white'),
                symbol='diamond'
            ),
            name='Pathways'
        )
        node_traces.append(pathway_trace)
    
    # Create figure
    fig = go.Figure(data=[edge_trace] + node_traces)
    
    fig.update_layout(
        title=dict(
            text='Gene-Pathway Interaction Network',
            x=0.5,
            font=dict(size=20, color='#3b82f6')
        ),
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
            bgcolor='rgba(255, 255, 255, 0.9)',
            bordercolor='#e2e8f0',
            borderwidth=1
        ),
        hovermode='closest',
        margin=dict(b=20, l=5, r=5, t=50),
        annotations=[
            dict(
                text="Blue: Genes | Green: Pathways",
                showarrow=False,
                xref="paper", yref="paper",
                x=0.5, y=-0.05,
                xanchor='center',
                font=dict(color='#64748b', size=12)
            )
        ],
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=600
    )
    
    return fig

# Light Plotly Theme
def apply_plotly_light_theme(fig: go.Figure) -> go.Figure:
    """Apply consistent light theme to Plotly figures"""
    fig.update_layout(
        paper_bgcolor='white',
        plot_bgcolor='white',
        font_color='#1e293b',
        font_family="Inter, sans-serif",
        title_font_size=18,
        title_font_color='#3b82f6',
        legend=dict(
            bgcolor='rgba(255, 255, 255, 0.9)',
            bordercolor='#e2e8f0',
            borderwidth=1,
            font=dict(color='#1e293b')
        ),
        margin=dict(l=50, r=50, t=50, b=50),
        hoverlabel=dict(
            bgcolor='white',
            font_color='#1e293b',
            bordercolor='#e2e8f0'
        )
    )
    
    # Update axes
    fig.update_xaxes(
        gridcolor='#f1f5f9',
        linecolor='#e2e8f0',
        linewidth=2,
        tickcolor='#94a3b8',
        title_font=dict(color='#64748b')
    )
    fig.update_yaxes(
        gridcolor='#f1f5f9',
        linecolor='#e2e8f0',
        linewidth=2,
        tickcolor='#94a3b8',
        title_font=dict(color='#64748b')
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
    sanitized = ''.join(c for c in gene.upper() if c.isalnum() or c in ['-', '_'])
    return sanitized[:50]

def validate_gene_list(genes: List[str]) -> Tuple[bool, str]:
    """Validate gene list before processing"""
    if not genes:
        return False, "Gene list is empty"
    
    if len(genes) > Config.MAX_GENES:
        return False, f"Too many genes (max {Config.MAX_GENES})"
    
    invalid_chars = set()
    for gene in genes:
        sanitized = sanitize_gene_symbol(gene)
        if sanitized != gene:
            invalid_chars.update(set(''.join(c for c in gene if not c.isalnum() and c not in ['-', '_'])))
    
    if invalid_chars:
        return False, f"Invalid characters found: {''.join(invalid_chars)}"
    
    return True, "Valid"

def set_entrez_email(email: str):
    """Set Entrez.email only once per session"""
    if not st.session_state.entrez_email_set:
        Entrez.email = email
        st.session_state.entrez_email_set = True
        logger.info("Entrez.email set for session")

def run_pathway_analysis(genes_from_input=None):
    """Run the pathway enrichment and drug discovery pipeline"""
    
    # Initialize session state
    if 'analysis_results' not in st.session_state:
        st.session_state.analysis_results = {}
    
    # Sidebar content
    with st.sidebar:
        st.markdown('<div class="sidebar-title">💊 Pathway & Drug Analyzer</div>', unsafe_allow_html=True)
        st.markdown('<div class="sidebar-tip">Complete gene-to-therapy pipeline for functional interpretation and drug discovery</div>', unsafe_allow_html=True)
        
        st.markdown("---")
        
        st.markdown("**📊 Analysis Pipeline**")
        st.markdown("""
        1. **Gene Mapping** → NCBI annotations
        2. **Pathway Analysis** → KEGG enrichment  
        3. **Disease Links** → OpenTargets associations
        4. **Drug Discovery** → Therapeutic compounds
        5. **Visualization** → Interactive networks
        """)
        
        st.markdown("---")
        
        st.markdown("**⚡ Performance Tips**")
        st.markdown('<div class="sidebar-tip">• Keep gene lists under 100 for optimal speed<br>• Results cached for 1 hour<br>• Built-in rate limiting prevents API timeouts</div>', unsafe_allow_html=True)
    
    # Main content
    st.markdown("""
    <div class="hero">
        <div class="hero-content">
            <h1>💊 Pathway & Drug Analyzer</h1>
            <p>Discover enriched pathways, disease associations, and potential therapeutic compounds for your gene list. 
            Perfect for interpreting DEG results or analyzing custom gene sets.</p>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # If genes are provided from DEGs analysis
    if genes_from_input:
        st.success(f"🎯 Using {len(genes_from_input)} genes from DEGs analysis")
        run_pathway_analysis_with_genes(genes_from_input)
        return
    
    # Input Section
    st.markdown('<div class="section-title">Input Configuration</div>', unsafe_allow_html=True)
    
    # Organism and email
    col1, col2 = st.columns(2)
    
    with col1:
        email = st.text_input(
            "📧 NCBI Email Address",
            placeholder="your.email@example.com",
            help="Required for NCBI API access"
        )
    
    with col2:
        organisms = {
            "🧑 Human (Homo sapiens)": {"entrez": "Homo sapiens", "kegg": "hsa"},
            "🐭 Mouse (Mus musculus)": {"entrez": "Mus musculus", "kegg": "mmu"},
            "🐀 Rat (Rattus norvegicus)": {"entrez": "Rattus norvegicus", "kegg": "rno"},
        }
        org_label = st.selectbox("🔬 Organism", list(organisms.keys()), index=0)
        organism_entrez = organisms[org_label]["entrez"]
        kegg_org_prefix = organisms[org_label]["kegg"]
    
    # Gene input methods
    st.markdown("### 🧬 Gene Input Methods")
    
    tab_upload, tab_manual = st.tabs(["📁 Upload File", "✏️ Manual Input"])
    
    genes_from_user = []
    
    with tab_upload:
        uploaded = st.file_uploader(
            "Upload gene list file",
            type=["csv", "tsv", "txt", "xlsx"],
            help="Supports CSV, TSV, XLSX, and text files",
            label_visibility="collapsed"
        )
        
        if uploaded is not None:
            genes_from_user = load_genes_from_any(uploaded)
            if genes_from_user:
                st.success(f"✅ Loaded {len(genes_from_user)} genes from file")
    
    with tab_manual:
        manual_input = st.text_area(
            "Enter gene symbols (one per line or comma-separated)",
            placeholder="TP53\nBRCA1\nEGFR\nMYC\nPIK3CA",
            height=150,
            label_visibility="collapsed"
        )
        
        if manual_input.strip():
            genes_from_user = (
                pd.Series([manual_input])
                .str.replace(r"[,;|\t\n ]+", "\n", regex=True)
                .str.split("\n").explode().str.strip().str.upper()
            )
            genes_from_user = [g for g in genes_from_user.tolist() if g and len(g) > 1][:Config.MAX_GENES]
            if genes_from_user:
                st.success(f"✅ Parsed {len(genes_from_user)} genes from input")
    
    # Drug filtering options
    st.markdown("### ⚙️ Analysis Filters")
    
    col_filter1, col_filter2 = st.columns(2)
    
    with col_filter1:
        opt_only_phase4 = st.checkbox(
            "Show only approved drugs (Phase 4+)",
            value=True,
            help="Filter to show only drugs that have completed clinical trials"
        )
    
    with col_filter2:
        show_investigational = st.checkbox(
            "Include investigational compounds",
            value=False,
            help="Include drugs in earlier clinical trial phases"
        )
    
    # Validate and run
    st.markdown("---")
    
    if genes_from_user:
        # Validate gene list
        valid_genes, invalid_genes = filter_valid_genes(genes_from_user)
        
        if invalid_genes:
            st.warning(f"⚠️ Filtered out {len(invalid_genes)} invalid gene symbols")
            if len(valid_genes) == 0:
                st.error("No valid gene symbols found")
                return
        
        # Show summary
        col_sum1, col_sum2, col_sum3 = st.columns(3)
        
        with col_sum1:
            st.metric("Valid Genes", len(valid_genes))
        
        with col_sum2:
            st.metric("Invalid Genes", len(invalid_genes))
        
        with col_sum3:
            if email and validate_email(email):
                st.metric("Email", "✅ Valid")
            else:
                st.metric("Email", "❌ Required")
        
        # Run button
        col_btn1, col_btn2 = st.columns([1, 2])
        
        with col_btn1:
            run_analysis = st.button(
                "🚀 Start Analysis", 
                type="primary", 
                disabled=(not valid_genes or not email),
                width="stretch"
            )
        
        with col_btn2:
            if not email:
                st.warning("Please provide your NCBI email address")
            elif not valid_genes:
                st.info("Please provide valid gene symbols")
    
    else:
        st.info("👆 Please provide gene symbols using one of the input methods above")
        run_analysis = False
    
    if run_analysis and email:
        set_entrez_email(email)
        run_pathway_analysis_with_genes(valid_genes, organism_entrez, kegg_org_prefix, 
                                       opt_only_phase4, show_investigational)

def run_pathway_analysis_with_genes(genes_from_input, organism_entrez=None, 
                                   kegg_org_prefix=None, opt_only_phase4=False, 
                                   show_investigational=False):
    """Run the actual pathway analysis with provided genes"""
    
    if organism_entrez is None:
        organism_entrez = "Homo sapiens"
    if kegg_org_prefix is None:
        kegg_org_prefix = "hsa"
    
    # Overall progress tracker
    overall_progress = st.progress(0)
    overall_status = st.empty()
    
    # Create tabs for different analysis sections
    tab_meta, tab_path, tab_disease, tab_drug, tab_viz = st.tabs([
        "📇 Gene Metadata", 
        "📊 Pathway Enrichment", 
        "🧬 Disease Links", 
        "💊 Drug Discovery", 
        "🌐 Network View"
    ])
    
    # Gene Metadata
    with tab_meta:
        st.markdown('<div class="section-title">Gene Annotations & Metadata</div>', unsafe_allow_html=True)
        
        if 'metadata' not in st.session_state.analysis_results:
            overall_status.text("🔄 Step 1/5: Fetching gene metadata and pathways...")
            
            with st.spinner("🔍 Fetching gene metadata and pathways..."):
                df_meta, pathway_to_genes = fetch_gene_metadata_and_kegg(
                    genes_from_input, organism_entrez, kegg_org_prefix
                )
            
            st.session_state.analysis_results['metadata'] = df_meta
            st.session_state.analysis_results['pathways'] = pathway_to_genes
            overall_progress.progress(20)
        else:
            df_meta = st.session_state.analysis_results['metadata']
            pathway_to_genes = st.session_state.analysis_results['pathways']
        
        if not df_meta.empty:
            # Summary cards
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
            
            # Data table - limit display rows
            st.markdown("### Gene Details")
            display_df = df_meta.copy()
            display_df.insert(0, "#", range(1, len(display_df) + 1))
            
            display_rows = min(Config.MAX_DISPLAY_ROWS, len(display_df))
            st.dataframe(
                display_df.head(display_rows), 
                width="stretch", 
                hide_index=True,
                column_config={
                    "Description": st.column_config.TextColumn(width="large"),
                    "KEGG_Pathways": st.column_config.TextColumn(width="large"),
                    "Status": st.column_config.TextColumn(width="small")
                }
            )
            
            if len(df_meta) > Config.MAX_DISPLAY_ROWS:
                st.info(f"Showing first {Config.MAX_DISPLAY_ROWS} of {len(df_meta)} genes. Download for full results.")
            
            # Download
            csv_data = df_meta.to_csv(index=False).encode("utf-8")
            st.download_button(
                "📥 Download Gene Metadata",
                data=csv_data,
                file_name="gene_metadata.csv",
                mime="text/csv",
                width="stretch"
            )
        else:
            st.error("❌ No gene metadata could be retrieved")
    
    # Pathway Enrichment
    with tab_path:
        st.markdown('<div class="section-title">Pathway Enrichment Analysis</div>', unsafe_allow_html=True)
        
        if 'pathways' in st.session_state.analysis_results:
            pathway_to_genes = st.session_state.analysis_results['pathways']
            
            overall_status.text("🔄 Step 2/5: Computing pathway enrichment...")
            
            with st.spinner("📈 Computing pathway enrichment..."):
                df_enrich = compute_enrichment_counts_only(pathway_to_genes)
            
            overall_progress.progress(40)
            
            if not df_enrich.empty:
                # Summary
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.metric("Enriched Pathways", len(df_enrich))
                with col2:
                    max_genes = df_enrich['Count'].max() if not df_enrich.empty else 0
                    st.metric("Max Genes/Pathway", max_genes)
                with col3:
                    avg_genes = df_enrich['Count'].mean() if not df_enrich.empty else 0
                    st.metric("Avg Genes/Pathway", f"{avg_genes:.1f}")
                
                # Data table - limit display rows
                st.markdown("### Enriched Pathways")
                display_enrich = df_enrich[['Pathway_ID', 'Pathway_Name', 'Count', 'Gene_List']].copy()
                display_enrich.insert(0, "#", range(1, len(display_enrich) + 1))
                
                display_rows = min(Config.MAX_DISPLAY_ROWS, len(display_enrich))
                st.dataframe(
                    display_enrich.head(display_rows), 
                    width="stretch", 
                    hide_index=True,
                    column_config={
                        "Pathway_Name": st.column_config.TextColumn(width="large"),
                        "Gene_List": st.column_config.TextColumn(width="large"),
                        "Count": st.column_config.NumberColumn(width="small")
                    }
                )
                
                if len(df_enrich) > Config.MAX_DISPLAY_ROWS:
                    st.info(f"Showing first {Config.MAX_DISPLAY_ROWS} of {len(df_enrich)} pathways. Download for full results.")
                
                # Visualization
                if len(df_enrich) > 0:
                    st.markdown("### 📊 Top Pathways Visualization")
                    
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
                    fig = apply_plotly_light_theme(fig)
                    
                    st.plotly_chart(fig, width="stretch")
                
                # Download
                download_enrich = df_enrich[['Pathway_ID', 'Pathway_Name', 'Count', 'Gene_List']].copy()
                csv_data = download_enrich.to_csv(index=False).encode("utf-8")
                st.download_button(
                    "📥 Download Pathway Results",
                    data=csv_data,
                    file_name="pathway_enrichment.csv",
                    mime="text/csv",
                    width="stretch"
                )
            else:
                st.info("ℹ️ No pathway enrichment found")
    
    # Disease Associations
    with tab_disease:
        st.markdown('<div class="section-title">Disease Associations</div>', unsafe_allow_html=True)
        
        if 'diseases' not in st.session_state.analysis_results:
            overall_status.text("🔄 Step 3/5: Mapping genes and fetching disease associations...")
            
            with st.spinner("🔍 Mapping genes to targets and fetching disease associations..."):
                gene_to_target = build_gene_to_ot_target_map(genes_from_input)
                df_diseases = collect_disease_links(gene_to_target)
                
            st.session_state.analysis_results['gene_to_target'] = gene_to_target
            st.session_state.analysis_results['diseases'] = df_diseases
            overall_progress.progress(60)
        else:
            gene_to_target = st.session_state.analysis_results['gene_to_target']
            df_diseases = st.session_state.analysis_results['diseases']
        
        if not df_diseases.empty:
            # Summary
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
            
            # Display - limit rows
            display_diseases = disease_summary.copy()
            display_diseases.insert(0, "#", range(1, len(display_diseases) + 1))
            
            display_rows = min(Config.MAX_DISPLAY_ROWS, len(display_diseases))
            st.dataframe(
                display_diseases.head(display_rows), 
                width="stretch", 
                hide_index=True,
                column_config={
                    "Disease_Name": st.column_config.TextColumn(width="large"),
                    "Therapeutic_Areas": st.column_config.TextColumn(width="medium"),
                    "Max_Score": st.column_config.NumberColumn(format="%.3f"),
                    "Avg_Score": st.column_config.NumberColumn(format="%.3f")
                }
            )
            
            if len(disease_summary) > Config.MAX_DISPLAY_ROWS:
                st.info(f"Showing first {Config.MAX_DISPLAY_ROWS} of {len(disease_summary)} diseases. Download for full results.")
            
            # Visualization
            if len(disease_summary) > 0:
                st.markdown("### 🎯 Disease Association Network")
                
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
                    },
                    color="Gene_Count",
                    color_continuous_scale="viridis"
                )
                
                fig = apply_plotly_light_theme(fig)
                st.plotly_chart(fig, width="stretch")
            
            # Download
            col_dl1, col_dl2 = st.columns(2)
            
            with col_dl1:
                csv_detailed = df_diseases.to_csv(index=False).encode("utf-8")
                st.download_button(
                    "📥 Download Detailed Associations",
                    data=csv_detailed,
                    file_name="disease_associations.csv",
                    mime="text/csv",
                    width="stretch"
                )
            
            with col_dl2:
                csv_summary = disease_summary.to_csv(index=False).encode("utf-8")
                st.download_button(
                    "📥 Download Disease Summary",
                    data=csv_summary,
                    file_name="disease_summary.csv",
                    mime="text/csv",
                    width="stretch"
                )
        else:
            st.info("ℹ️ No disease associations found")
    
    # Drug Suggestions
    with tab_drug:
        st.markdown('<div class="section-title">Therapeutic Drug Discovery</div>', unsafe_allow_html=True)
        
        if 'drugs' not in st.session_state.analysis_results:
            if 'gene_to_target' in st.session_state.analysis_results:
                gene_to_target = st.session_state.analysis_results['gene_to_target']
                
                overall_status.text("🔄 Step 4/5: Collecting drug suggestions...")
                
                with st.spinner("💊 Collecting drug suggestions..."):
                    df_drugs = collect_drug_suggestions(gene_to_target)
                
                st.session_state.analysis_results['drugs'] = df_drugs
                overall_progress.progress(80)
            else:
                st.error("❌ Complete disease association step first")
                df_drugs = pd.DataFrame()
        else:
            df_drugs = st.session_state.analysis_results['drugs']
        
        if not df_drugs.empty:
            # Apply filters
            filtered_drugs = df_drugs.copy()
            
            st.markdown(f"**Found {len(filtered_drugs)} total drug entries**")
            
            if opt_only_phase4:
                # Show drugs that have reached phase 4 or higher (approved drugs)
                filtered_drugs = filtered_drugs[
                    (filtered_drugs['max_phase_numeric'] >= 4)
                ]
                st.info(f"Showing only approved drugs (Phase 4+): {len(filtered_drugs)} drugs")
            
            if not show_investigational:
                filtered_drugs = filtered_drugs[
                    filtered_drugs['status'] != 'Investigational'
                ]
            
            if not filtered_drugs.empty:
                # Calculate metrics
                valid_max_phases = filtered_drugs[filtered_drugs['max_phase_numeric'] > 0]['max_phase_numeric']
                avg_phase = valid_max_phases.mean() if len(valid_max_phases) > 0 else 0
                
                # Summary
                col1, col2, col3, col4 = st.columns(4)
                
                with col1:
                    unique_drugs = filtered_drugs['drug_name'].nunique()
                    st.metric("Unique Drugs", unique_drugs)
                with col2:
                    unique_targets = filtered_drugs['gene'].nunique()
                    st.metric("Targeted Genes", unique_targets)
                with col3:
                    st.metric("Avg Phase", f"{avg_phase:.1f}")
                with col4:
                    approved_drugs = len(filtered_drugs[
                        (filtered_drugs['phase_numeric'] >= 4) | 
                        (filtered_drugs['max_phase_numeric'] >= 4)
                    ])
                    st.metric("Approved Drugs", approved_drugs)
                
                # Data table - limit display rows
                display_columns = [
                    'gene', 'gene_symbol', 'drug_id', 'drug_name', 'drug_type', 
                    'phase', 'status', 'max_phase', 'moa', 'disease_name', 'therapeutic_areas'
                ]
                
                available_columns = [col for col in display_columns if col in filtered_drugs.columns]
                display_drugs = filtered_drugs[available_columns].copy()
                display_drugs.insert(0, "#", range(1, len(display_drugs) + 1))
                
                display_rows = min(Config.MAX_DISPLAY_ROWS, len(display_drugs))
                st.dataframe(
                    display_drugs.head(display_rows), 
                    width="stretch",
                    column_config={
                        "#": st.column_config.NumberColumn(width="small"),
                        "gene": st.column_config.TextColumn("Gene", width="small"),
                        "gene_symbol": st.column_config.TextColumn("Symbol", width="small"),
                        "drug_name": st.column_config.TextColumn("Drug Name", width="medium"),
                        "drug_type": st.column_config.TextColumn("Type", width="small"),
                        "phase": st.column_config.TextColumn("Phase", width="small"),
                        "status": st.column_config.TextColumn("Status", width="small"),
                        "max_phase": st.column_config.TextColumn("Max Phase", width="small"),
                        "moa": st.column_config.TextColumn("Mechanism", width="large"),
                        "disease_name": st.column_config.TextColumn("Disease", width="medium"),
                        "therapeutic_areas": st.column_config.TextColumn("Areas", width="medium"),
                    }
                )
                
                if len(filtered_drugs) > Config.MAX_DISPLAY_ROWS:
                    st.info(f"Showing first {Config.MAX_DISPLAY_ROWS} of {len(filtered_drugs)} drugs. Download for full results.")
                
                # Visualization
                if len(filtered_drugs) > 0:
                    st.markdown("### 📈 Drug Development Phase Distribution")
                    
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
                    
                    fig = apply_plotly_light_theme(fig)
                    st.plotly_chart(fig, width="stretch")
                
                # Download
                download_columns = [col for col in display_columns if col in filtered_drugs.columns]
                download_data = filtered_drugs[download_columns].to_csv(index=False).encode("utf-8")
                st.download_button(
                    "📥 Download Drug Suggestions",
                    data=download_data,
                    file_name="drug_suggestions.csv",
                    mime="text/csv",
                    width="stretch"
                )
            else:
                st.info("ℹ️ No drugs match the filter criteria")
        else:
            st.info("ℹ️ No drug suggestions available")
    
    # Network Visualization
    with tab_viz:
        st.markdown('<div class="section-title">🌐 Interactive Network Visualization</div>', unsafe_allow_html=True)
        
        # Check if we have the required data
        has_diseases = 'diseases' in st.session_state.analysis_results and not st.session_state.analysis_results['diseases'].empty
        has_drugs = 'drugs' in st.session_state.analysis_results and not st.session_state.analysis_results['drugs'].empty
        has_pathways = 'pathways' in st.session_state.analysis_results and st.session_state.analysis_results['pathways']
        
        overall_status.text("🔄 Step 5/5: Generating visualizations...")
        overall_progress.progress(100)
        overall_status.text("✅ Analysis Complete!")
        time.sleep(1)
        overall_progress.empty()
        overall_status.empty()
        
        if has_diseases or has_drugs or has_pathways:
            col_left, col_right = st.columns(2)
            
            # Gene-Disease-Drug Network
            with col_left:
                if has_diseases:
                    st.markdown("### 🔗 Gene-Disease-Drug Network")
                    st.markdown("*Visualizing relationships between genes, diseases, and potential drugs*")
                    
                    try:
                        df_diseases = st.session_state.analysis_results['diseases']
                        df_drugs = st.session_state.analysis_results.get('drugs', pd.DataFrame())
                        
                        # Check if network size is manageable
                        total_nodes_estimate = df_diseases['gene'].nunique() + df_diseases['disease_name'].nunique()
                        if df_drugs is not None and not df_drugs.empty:
                            total_nodes_estimate += df_drugs['drug_name'].nunique()
                        
                        if total_nodes_estimate > Config.MAX_NETWORK_NODES:
                            st.warning(f"⚠️ Network too large ({total_nodes_estimate} estimated nodes). Showing summary instead.")
                            # Show summary instead
                            with st.expander("📊 Network Statistics"):
                                col_stat1, col_stat2, col_stat3 = st.columns(3)
                                with col_stat1:
                                    unique_genes = df_diseases['gene'].nunique()
                                    st.metric("Unique Genes", unique_genes)
                                with col_stat2:
                                    unique_diseases = df_diseases['disease_name'].nunique()
                                    st.metric("Diseases", unique_diseases)
                                with col_stat3:
                                    unique_drugs = df_drugs['drug_name'].nunique() if not df_drugs.empty else 0
                                    st.metric("Drugs", unique_drugs)
                        else:
                            fig_network = create_gene_disease_drug_network(df_diseases, df_drugs)
                            
                            if fig_network:
                                fig_network = apply_plotly_light_theme(fig_network)
                                st.plotly_chart(fig_network, width="stretch")
                                
                                # Show network statistics
                                with st.expander("📊 Network Statistics"):
                                    col_stat1, col_stat2, col_stat3 = st.columns(3)
                                    with col_stat1:
                                        unique_genes = df_diseases['gene'].nunique()
                                        st.metric("Unique Genes", unique_genes)
                                    with col_stat2:
                                        unique_diseases = df_diseases['disease_name'].nunique()
                                        st.metric("Diseases", unique_diseases)
                                    with col_stat3:
                                        unique_drugs = df_drugs['drug_name'].nunique() if not df_drugs.empty else 0
                                        st.metric("Drugs", unique_drugs)
                            else:
                                st.info("Not enough data to create network visualization. Try with more genes.")
                                
                    except Exception as e:
                        st.error(f"Error creating gene-disease-drug network: {e}")
            
            # Gene-Pathway Interaction Network
            with col_right:
                if has_pathways:
                    st.markdown("### 🛤️ Gene-Pathway Interaction Network")
                    st.markdown("*Visualizing connections between genes and biological pathways*")
                    
                    try:
                        pathway_to_genes = st.session_state.analysis_results['pathways']
                        
                        # Check if network size is manageable
                        total_genes = sum(len(genes) for genes in pathway_to_genes.values())
                        if total_genes > Config.MAX_NETWORK_NODES:
                            st.warning(f"⚠️ Network too large ({total_genes} genes). Showing summary instead.")
                            # Show summary instead
                            with st.expander("📊 Pathway Statistics"):
                                col_stat1, col_stat2 = st.columns(2)
                                with col_stat1:
                                    total_pathways = len(pathway_to_genes)
                                    st.metric("Total Pathways", total_pathways)
                                with col_stat2:
                                    avg_genes = total_genes / max(1, total_pathways)
                                    st.metric("Avg Genes/Pathway", f"{avg_genes:.1f}")
                        else:
                            fig_pathway = create_gene_pathway_interaction_network(pathway_to_genes, kegg_org_prefix)
                            
                            if fig_pathway:
                                fig_pathway = apply_plotly_light_theme(fig_pathway)
                                st.plotly_chart(fig_pathway, width="stretch")
                                
                                # Show pathway statistics
                                with st.expander("📊 Pathway Statistics"):
                                    col_stat1, col_stat2 = st.columns(2)
                                    with col_stat1:
                                        total_pathways = len(pathway_to_genes)
                                        st.metric("Total Pathways", total_pathways)
                                    with col_stat2:
                                        avg_genes = total_genes / max(1, total_pathways)
                                        st.metric("Avg Genes/Pathway", f"{avg_genes:.1f}")
                                    
                                    # Show top pathways
                                    pathway_counts = []
                                    for pid, genes in pathway_to_genes.items():
                                        pathway_name = kegg_pathway_name(pid) or pid.replace("path:", "")
                                        pathway_counts.append({
                                            'Pathway': pathway_name[:50] + "..." if len(pathway_name) > 50 else pathway_name,
                                            'Gene Count': len(genes)
                                        })
                                    
                                    if pathway_counts:
                                        df_pathway_counts = pd.DataFrame(pathway_counts)
                                        df_pathway_counts = df_pathway_counts.sort_values('Gene Count', ascending=False).head(5)
                                        st.dataframe(df_pathway_counts, width="stretch", hide_index=True)
                            else:
                                st.info("Not enough pathway data to create network visualization.")
                                
                    except Exception as e:
                        st.error(f"Error creating gene-pathway network: {e}")
            
            # Additional insights
            st.markdown("---")
            st.markdown("### 📋 Network Insights")
            
            if has_diseases and has_drugs:
                col_insight1, col_insight2 = st.columns(2)
                
                with col_insight1:
                    st.markdown("**🧬 Key Gene Connections**")
                    if not df_diseases.empty:
                        # Find genes with most disease associations
                        gene_disease_counts = df_diseases.groupby('gene').agg({
                            'disease_name': 'nunique',
                            'association_score': 'mean'
                        }).sort_values('disease_name', ascending=False).head(5)
                        
                        for gene, (count, avg_score) in gene_disease_counts.iterrows():
                            st.markdown(f"• **{gene}**: {int(count)} diseases (avg score: {avg_score:.3f})")
                
                with col_insight2:
                    st.markdown("**💊 Top Therapeutic Targets**")
                    if not df_drugs.empty:
                        # Find genes with most drugs
                        gene_drug_counts = df_drugs.groupby('gene').agg({
                            'drug_name': 'nunique',
                            'phase': lambda x: x.max() if not x.empty else 'N/A'
                        }).sort_values('drug_name', ascending=False).head(5)
                        
                        for gene, (count, phase) in gene_drug_counts.iterrows():
                            st.markdown(f"• **{gene}**: {int(count)} drugs (max phase: {phase})")
        
        else:
            st.info("🔄 Complete previous analysis steps to generate network visualizations.")
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: #64748b; font-size: 0.9rem; padding: 1rem;'>
        <p><strong>Gene2Therapy</strong> • Modern light theme interface with optimized performance</p>
        <p>Data sources: NCBI, KEGG, OpenTargets • Results cached for 1 hour</p>
        <p>⚠️ Always validate findings with primary literature and clinical sources</p>
    </div>
    """, unsafe_allow_html=True)

# =============================================================================
# MAIN APP INTEGRATION
# =============================================================================

def main():
    # Initialize session state
    initialize_session_state()
    
    # Main title with gradient
    st.markdown("""
    <div style="text-align: center; padding: 2rem 0;">
        <h1 style="
            font-size: 3.5rem;
            font-weight: 800;
            background: linear-gradient(135deg, #3b82f6 0%, #10b981 100%);
            -webkit-background-clip: text;
            background-clip: text;
            color: transparent;
            margin-bottom: 0.5rem;
            font-family: 'Source Sans Pro', sans-serif;
        ">
            🧬 Gene2Therapy
        </h1>
        <p style="
            color: #64748b;
            font-size: 1.2rem;
            max-width: 800px;
            margin: 0 auto 2rem auto;
            line-height: 1.6;
        ">
            Integrated pipeline for differential expression analysis, pathway enrichment, 
            disease association, and drug discovery
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Mode selection with feature cards
    st.markdown('<div class="section-title">Select Analysis Mode</div>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="feature-card">
            <div class="mode-icon">🧬</div>
            <div class="feature-title">DEGs Analysis</div>
            <div class="feature-desc">
                Identify differentially expressed genes from RNA-seq data.
                Upload count matrices, configure statistical parameters, and
                visualize results with interactive plots.
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        if st.button("Start DEGs Analysis", key="degs_btn", width="stretch"):
            st.session_state.analysis_mode = 'degs_only'
            st.session_state.current_pipeline_step = 'degs'
            st.session_state.show_pathway_analysis = False
            st.session_state.degs_completed = False
            st.session_state.degs_running = False
            st.rerun()
    
    with col2:
        st.markdown("""
        <div class="feature-card">
            <div class="mode-icon">💊</div>
            <div class="feature-title">Pathway & Drug Analysis</div>
            <div class="feature-desc">
                Discover enriched pathways, disease associations, and
                therapeutic compounds. Use DEG results or upload custom
                gene lists for comprehensive analysis.
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        if st.button("Start Pathway Analysis", key="pathway_btn", width="stretch"):
            st.session_state.analysis_mode = 'pathway_only'
            st.session_state.current_pipeline_step = 'pathway'
            st.session_state.show_pathway_analysis = True
            st.session_state.degs_completed = False
            st.session_state.degs_running = False
            st.rerun()
    
    # Display appropriate pipeline
    st.markdown("---")
    
    if st.session_state.analysis_mode == 'degs_only':
        st.markdown('<div style="text-align: center; margin: 2rem 0;"><h2>🧬 Differential Expression Analysis</h2></div>', unsafe_allow_html=True)
        run_degs_analysis()
    
    elif st.session_state.analysis_mode == 'pathway_only':
        st.markdown('<div style="text-align: center; margin: 2rem 0;"><h2>💊 Pathway & Drug Analysis</h2></div>', unsafe_allow_html=True)
        run_pathway_analysis()
    
    else:
        # Welcome and feature overview
        st.markdown("""
        <div style="text-align: center; padding: 3rem 0;">
            <h2 style="color: #3b82f6; margin-bottom: 2rem;">Welcome to Gene2Therapy</h2>
            <p style="color: #64748b; font-size: 1.1rem; max-width: 800px; margin: 0 auto 3rem auto;">
                A comprehensive bioinformatics platform for gene expression analysis and therapeutic discovery.
                Choose one of the analysis modes above to get started.
            </p>
        </div>
        """, unsafe_allow_html=True)
        
        # Feature highlights
        col_feat1, col_feat2, col_feat3 = st.columns(3)
        
        with col_feat1:
            st.markdown("""
            <div class="feature-card">
                <div style="font-size: 2rem; margin-bottom: 1rem;">⚡</div>
                <div class="feature-title">Fast Processing</div>
                <div class="feature-desc">
                    Optimized algorithms and caching for rapid analysis of large datasets
                </div>
            </div>
            """, unsafe_allow_html=True)
        
        with col_feat2:
            st.markdown("""
            <div class="feature-card">
                <div style="font-size: 2rem; margin-bottom: 1rem;">🔗</div>
                <div class="feature-title">Integrated Databases</div>
                <div class="feature-desc">
                    Access NCBI, KEGG, and OpenTargets databases for comprehensive analysis
                </div>
            </div>
            """, unsafe_allow_html=True)
        
        with col_feat3:
            st.markdown("""
            <div class="feature-card">
                <div style="font-size: 2rem; margin-bottom: 1rem;">📊</div>
                <div class="feature-title">Interactive Visualization</div>
                <div class="feature-desc">
                    Beautiful, interactive plots and networks for data exploration
                </div>
            </div>
            """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()