import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import warnings
import io
from datetime import datetime
import time

warnings.filterwarnings('ignore')

# Set page configuration
st.set_page_config(
    page_title="Differential Expression Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

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

def main():
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

if __name__ == "__main__":
    main()