import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import ttest_ind

def benjamini_hochberg(pvals):
    """Apply Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    sorted_idx = np.argsort(pvals)
    sorted_pvals = np.array(pvals)[sorted_idx]
    
    # BH correction
    adjusted = np.zeros(n)
    for i, p in enumerate(sorted_pvals):
        adjusted[i] = p * n / (i + 1)
    
    # Ensure monotonicity
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])
    
    # Clip to [0, 1]
    adjusted = np.clip(adjusted, 0, 1)
    
    # Reorder to original indices
    result = np.zeros(n)
    for i, idx in enumerate(sorted_idx):
        result[idx] = adjusted[i]
    
    return result.tolist()

def process_bulk_rna_seq(file_path, is_normalized):
    """
    Process Bulk RNA-seq data.
    Expected format: Genes (rows) x Samples (columns).
    First column should be gene names (index).
    """
    try:
        # Load data
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path, index_col=0)
        else:
            df = pd.read_csv(file_path, sep='\t', index_col=0)
        
        # Clean data - remove any non-numeric columns
        df = df.select_dtypes(include=[np.number])
        
        # Remove rows with all zeros
        df = df.loc[(df != 0).any(axis=1)]
        
        # Remove rows with NaN
        df = df.dropna()
        
        if df.shape[0] == 0:
            raise ValueError("No valid genes found after filtering.")
        
        if df.shape[1] < 2:
            raise ValueError("Need at least 2 samples for DEG analysis.")
        
        # Split samples: assume first half is Control, second half is Case
        n_samples = df.shape[1]
        mid = n_samples // 2
        
        if mid == 0:
            mid = 1
        
        control_cols = df.columns[:mid].tolist()
        case_cols = df.columns[mid:].tolist()
        
        controls = df[control_cols]
        cases = df[case_cols]
        
        if is_normalized:
            # For normalized data: compute log2FC and t-test
            case_mean = cases.mean(axis=1)
            ctrl_mean = controls.mean(axis=1)
            
            # Check if data appears to be log2 transformed
            max_val = df.max().max()
            
            if max_val > 50:
                # Data is likely linear scale, need to log2 transform
                # Add pseudocount to avoid log(0)
                case_log = np.log2(cases + 1)
                ctrl_log = np.log2(controls + 1)
                log2fc = case_log.mean(axis=1) - ctrl_log.mean(axis=1)
            else:
                # Data is likely already log2 transformed
                log2fc = case_mean - ctrl_mean
            
            # T-test for each gene
            pvals = []
            for gene in df.index:
                try:
                    stat, p = ttest_ind(cases.loc[gene].values, controls.loc[gene].values)
                    pvals.append(p if not np.isnan(p) else 1.0)
                except:
                    pvals.append(1.0)
            
            # BH correction
            padj = benjamini_hochberg(pvals)
            
            results_df = pd.DataFrame({
                'gene': df.index,
                'log2FoldChange': log2fc.values,
                'pvalue': pvals,
                'padj': padj
            })
            
        else:
            # For raw counts: use pyDESeq2
            try:
                from pydeseq2.dds import DeseqDataSet
                from pydeseq2.ds import DeseqStats
                
                # Create metadata
                sample_names = df.columns.tolist()
                conditions = ['control'] * mid + ['case'] * (n_samples - mid)
                
                metadata = pd.DataFrame({
                    'condition': conditions
                }, index=sample_names)
                
                # Transpose: pyDESeq2 expects samples as rows, genes as columns
                counts_df = df.T.astype(int)
                
                # Create DeseqDataSet
                dds = DeseqDataSet(
                    counts=counts_df,
                    metadata=metadata,
                    design_factors="condition",
                    ref_level=["condition", "control"]
                )
                
                # Fit model
                dds.deseq2()
                
                # Get results
                stat_res = DeseqStats(dds, contrast=["condition", "case", "control"])
                stat_res.summary()
                
                results = stat_res.results_df
                
                results_df = pd.DataFrame({
                    'gene': results.index,
                    'log2FoldChange': results['log2FoldChange'].values,
                    'pvalue': results['pvalue'].values,
                    'padj': results['padj'].fillna(1.0).values
                })
                
            except ImportError:
                # Fallback if pyDESeq2 not installed: use simple method
                # Log2 transform counts
                case_log = np.log2(cases + 1)
                ctrl_log = np.log2(controls + 1)
                log2fc = case_log.mean(axis=1) - ctrl_log.mean(axis=1)
                
                # T-test
                pvals = []
                for gene in df.index:
                    try:
                        case_vals = np.log2(cases.loc[gene].values + 1)
                        ctrl_vals = np.log2(controls.loc[gene].values + 1)
                        stat, p = ttest_ind(case_vals, ctrl_vals)
                        pvals.append(p if not np.isnan(p) else 1.0)
                    except:
                        pvals.append(1.0)
                
                padj = benjamini_hochberg(pvals)
                
                results_df = pd.DataFrame({
                    'gene': df.index,
                    'log2FoldChange': log2fc.values,
                    'pvalue': pvals,
                    'padj': padj
                })
        
        # Sort by adjusted p-value
        results_df = results_df.sort_values('padj')
        
        # Replace NaN with defaults
        results_df['log2FoldChange'] = results_df['log2FoldChange'].fillna(0)
        results_df['pvalue'] = results_df['pvalue'].fillna(1)
        results_df['padj'] = results_df['padj'].fillna(1)
        
        # Convert to list of dicts for JSON serialization
        serialized_results = results_df.to_dict(orient='records')
        
        return serialized_results

    except Exception as e:
        import traceback
        error_msg = f"Processing failed: {str(e)}\n{traceback.format_exc()}"
        raise ValueError(error_msg)


def process_scrna_seq(file_path):
    """
    Process scRNA-seq data using Scanpy.
    """
    try:
        import scanpy as sc
        
        # Read data - support multiple formats
        if file_path.endswith('.h5ad'):
            adata = sc.read_h5ad(file_path)
        elif file_path.endswith('.csv'):
            adata = sc.read_csv(file_path)
        elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
            adata = sc.read_csv(file_path, delimiter='\t')
        else:
            # Try reading as CSV first
            try:
                adata = sc.read_csv(file_path)
            except:
                adata = sc.read_csv(file_path, delimiter='\t')
        
        # Basic QC filtering
        if adata.n_obs > 100:
            sc.pp.filter_cells(adata, min_genes=200)
        if adata.n_vars > 100:
            sc.pp.filter_genes(adata, min_cells=3)
        
        # Normalize
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Find variable genes
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        
        # PCA
        sc.tl.pca(adata, svd_solver='arpack')
        
        # Neighbors and clustering
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=min(40, adata.n_vars - 1))
        sc.tl.leiden(adata)
        
        # Find marker genes
        sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
        
        # Extract results
        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        
        output = []
        for group in groups:
            genes = result['names'][group]
            scores = result['scores'][group]
            pvals = result['pvals_adj'][group]
            logfc = result['logfoldchanges'][group] if 'logfoldchanges' in result else [0] * len(genes)
            
            for i, (g, s, p) in enumerate(zip(genes[:20], scores[:20], pvals[:20])):
                lfc = logfc[i] if i < len(logfc) else 0
                output.append({
                    'gene': str(g),
                    'cluster': str(group),
                    'log2FoldChange': float(lfc) if not np.isnan(lfc) else 0,
                    'pvalue': float(p) if not np.isnan(p) else 1.0,
                    'padj': float(p) if not np.isnan(p) else 1.0,
                    'score': float(s) if not np.isnan(s) else 0
                })
        
        return output

    except Exception as e:
        import traceback
        error_msg = f"scRNA processing failed: {str(e)}\n{traceback.format_exc()}"
        raise ValueError(error_msg)
