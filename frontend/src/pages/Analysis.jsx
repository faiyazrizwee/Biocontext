import React, { useState, useCallback } from 'react';
import { useDropzone } from 'react-dropzone';
import Plot from 'react-plotly.js';
import { degApi } from '../api';
import './Analysis.css';

const Analysis = () => {
    const [file, setFile] = useState(null);
    const [jobType, setJobType] = useState('BULK');
    const [isNormalized, setIsNormalized] = useState(true);
    const [job, setJob] = useState(null);
    const [results, setResults] = useState(null);
    const [loading, setLoading] = useState(false);
    const [progress, setProgress] = useState({ step: 0, status: '', message: '' });
    const [error, setError] = useState(null);
    const [activeTab, setActiveTab] = useState('volcano');

    const onDrop = useCallback((acceptedFiles) => {
        if (acceptedFiles.length > 0) {
            setFile(acceptedFiles[0]);
            setError(null);
        }
    }, []);

    const { getRootProps, getInputProps, isDragActive } = useDropzone({
        onDrop,
        accept: {
            'text/csv': ['.csv'],
            'text/tab-separated-values': ['.tsv', '.txt']
        },
        maxFiles: 1
    });

    const handleSubmit = async (e) => {
        e.preventDefault();
        if (!file) {
            setError('Please upload a file');
            return;
        }

        setLoading(true);
        setError(null);
        setProgress({ step: 1, status: 'Uploading', message: 'Uploading data file...' });

        try {
            const formData = new FormData();
            formData.append('data_file', file);
            formData.append('job_type', jobType);
            formData.append('is_normalized', isNormalized);

            const jobData = await degApi.createJob(formData);
            setJob(jobData);
            setProgress({ step: 2, status: 'Processing', message: 'Running DEG analysis...' });
            pollResults(jobData.id);
        } catch (err) {
            setError(err.response?.data?.error || 'Failed to create analysis job');
            setLoading(false);
        }
    };

    const pollResults = async (jobId) => {
        const maxAttempts = 60;
        let attempts = 0;

        const poll = async () => {
            try {
                const jobData = await degApi.getJob(jobId);
                setJob(jobData);

                if (jobData.status === 'COMPLETED') {
                    setProgress({ step: 3, status: 'Loading', message: 'Fetching results...' });
                    const resultsData = await degApi.getResults(jobId);
                    setResults(resultsData);
                    setLoading(false);
                    setProgress({ step: 0, status: '', message: '' });
                } else if (jobData.status === 'FAILED') {
                    setError(jobData.error_message || 'Analysis failed');
                    setLoading(false);
                    setProgress({ step: 0, status: '', message: '' });
                } else if (attempts < maxAttempts) {
                    attempts++;
                    setTimeout(poll, 2000);
                } else {
                    setError('Analysis timed out');
                    setLoading(false);
                }
            } catch (err) {
                setError('Failed to fetch results');
                setLoading(false);
            }
        };

        poll();
    };

    const getVolcanoData = () => {
        if (!results || !Array.isArray(results)) return [];

        const upregulated = results.filter(r => {
            const fc = r.log2FoldChange || 0;
            const pval = r.padj || r.pvalue || 1;
            return pval < 0.05 && fc > 1;
        });

        const downregulated = results.filter(r => {
            const fc = r.log2FoldChange || 0;
            const pval = r.padj || r.pvalue || 1;
            return pval < 0.05 && fc < -1;
        });

        const nonsig = results.filter(r => {
            const fc = r.log2FoldChange || 0;
            const pval = r.padj || r.pvalue || 1;
            return !(pval < 0.05 && Math.abs(fc) > 1);
        });

        return [
            {
                x: nonsig.map(r => r.log2FoldChange || 0),
                y: nonsig.map(r => -Math.log10(r.padj || r.pvalue || 1)),
                text: nonsig.map(r => r.gene || ''),
                mode: 'markers',
                type: 'scatter',
                name: 'Not Significant',
                marker: {
                    color: '#94a3b8',
                    size: 8,
                    opacity: 0.6,
                    line: { width: 1, color: '#64748b' }
                },
                hovertemplate: '<b>%{text}</b><br>log2FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
            },
            {
                x: upregulated.map(r => r.log2FoldChange || 0),
                y: upregulated.map(r => -Math.log10(r.padj || r.pvalue || 1)),
                text: upregulated.map(r => r.gene || ''),
                mode: 'markers',
                type: 'scatter',
                name: 'Upregulated',
                marker: {
                    color: '#10b981',
                    size: 10,
                    opacity: 0.85,
                    symbol: 'circle',
                    line: { width: 1.5, color: '#059669' }
                },
                hovertemplate: '<b>%{text}</b><br>log2FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
            },
            {
                x: downregulated.map(r => r.log2FoldChange || 0),
                y: downregulated.map(r => -Math.log10(r.padj || r.pvalue || 1)),
                text: downregulated.map(r => r.gene || ''),
                mode: 'markers',
                type: 'scatter',
                name: 'Downregulated',
                marker: {
                    color: '#ef4444',
                    size: 10,
                    opacity: 0.85,
                    symbol: 'circle',
                    line: { width: 1.5, color: '#dc2626' }
                },
                hovertemplate: '<b>%{text}</b><br>log2FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
            }
        ];
    };

    const volcanoLayout = {
        title: {
            text: '<b>Volcano Plot</b>',
            font: { size: 20, color: '#134e4a', family: 'Inter, sans-serif' }
        },
        xaxis: {
            title: { text: 'log₂ Fold Change', font: { size: 14, color: '#3f6b68' } },
            zeroline: true,
            zerolinecolor: '#a7f3d0',
            zerolinewidth: 2,
            gridcolor: '#e0f2f1',
            tickfont: { color: '#3f6b68' }
        },
        yaxis: {
            title: { text: '-log₁₀(adjusted p-value)', font: { size: 14, color: '#3f6b68' } },
            gridcolor: '#e0f2f1',
            tickfont: { color: '#3f6b68' }
        },
        paper_bgcolor: '#ffffff',
        plot_bgcolor: '#f8fffe',
        font: { family: 'Inter, sans-serif' },
        hovermode: 'closest',
        legend: {
            orientation: 'h',
            yanchor: 'bottom',
            y: 1.02,
            xanchor: 'center',
            x: 0.5,
            bgcolor: 'rgba(255,255,255,0.9)',
            bordercolor: '#a7f3d0',
            borderwidth: 1
        },
        shapes: [
            { type: 'line', x0: -1, x1: -1, y0: 0, y1: 1, yref: 'paper', line: { dash: 'dot', color: '#0ea5e9', width: 2 } },
            { type: 'line', x0: 1, x1: 1, y0: 0, y1: 1, yref: 'paper', line: { dash: 'dot', color: '#0ea5e9', width: 2 } },
            { type: 'line', x0: 0, x1: 1, xref: 'paper', y0: -Math.log10(0.05), y1: -Math.log10(0.05), line: { dash: 'dot', color: '#f59e0b', width: 2 } }
        ],
        annotations: [
            { x: 1.05, y: 0.5, xref: 'paper', yref: 'paper', text: 'FC > 2', showarrow: false, font: { size: 10, color: '#0ea5e9' } },
            { x: -0.05, y: 0.5, xref: 'paper', yref: 'paper', text: 'FC < -2', showarrow: false, font: { size: 10, color: '#0ea5e9' } }
        ],
        margin: { t: 80, r: 60, b: 60, l: 70 }
    };

    // Get heatmap data for a specific direction (up, down, or all)
    const getHeatmapData = (direction = 'all') => {
        if (!results || !Array.isArray(results)) return [];

        let filteredGenes = [...results];

        if (direction === 'up') {
            filteredGenes = results.filter(r => (r.log2FoldChange || 0) > 0 && (r.padj || r.pvalue || 1) < 0.05);
        } else if (direction === 'down') {
            filteredGenes = results.filter(r => (r.log2FoldChange || 0) < 0 && (r.padj || r.pvalue || 1) < 0.05);
        }

        // Sort by absolute log2FC and take top 25
        const sortedGenes = filteredGenes
            .sort((a, b) => Math.abs(b.log2FoldChange || 0) - Math.abs(a.log2FoldChange || 0))
            .slice(0, 25);

        if (sortedGenes.length === 0) return [];

        const z = sortedGenes.map(r => [r.log2FoldChange || 0]);
        const y = sortedGenes.map(r => r.gene || '');

        const colorscale = direction === 'up'
            ? [[0, '#a7f3d0'], [1, '#059669']]  // Green gradient for upregulated
            : direction === 'down'
                ? [[0, '#dc2626'], [1, '#fecaca']]  // Red gradient for downregulated
                : [[0, '#3b82f6'], [0.5, '#fefefe'], [1, '#ef4444']];  // Blue-White-Red for all

        return [{
            z: z,
            y: y,
            x: ['log₂FC'],
            type: 'heatmap',
            colorscale: colorscale,
            zmid: direction === 'all' ? 0 : undefined,
            showscale: true,
            colorbar: {
                title: {
                    text: 'log₂FC',
                    font: { size: 12, color: '#3f6b68' }
                },
                tickfont: { color: '#3f6b68' },
                thickness: 15,
                len: 0.9
            },
            hovertemplate: '<b>%{y}</b><br>log₂FC: %{z:.3f}<extra></extra>'
        }];
    };

    const getHeatmapLayout = (title) => ({
        title: {
            text: `<b>${title}</b>`,
            font: { size: 16, color: '#134e4a', family: 'Inter, sans-serif' }
        },
        paper_bgcolor: '#ffffff',
        plot_bgcolor: '#ffffff',
        font: { family: 'Inter, sans-serif' },
        xaxis: {
            tickfont: { size: 12, color: '#3f6b68' }
        },
        yaxis: {
            tickfont: { size: 9, color: '#3f6b68' },
            automargin: true
        },
        margin: { t: 50, r: 80, b: 30, l: 120 }
    });

    const getTopGenes = (direction) => {
        if (!results || !Array.isArray(results)) return [];

        return results
            .filter(r => direction === 'up' ? (r.log2FoldChange || 0) > 0 : (r.log2FoldChange || 0) < 0)
            .sort((a, b) => Math.abs(b.log2FoldChange || 0) - Math.abs(a.log2FoldChange || 0))
            .slice(0, 10);
    };

    const downloadCSV = (data, filename) => {
        if (!data.length) return;
        const headers = Object.keys(data[0]).join(',');
        const rows = data.map(row => Object.values(row).join(','));
        const csv = [headers, ...rows].join('\n');
        const blob = new Blob([csv], { type: 'text/csv' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        a.click();
        URL.revokeObjectURL(url);
    };

    const plotConfig = {
        responsive: true,
        displayModeBar: true,
        modeBarButtonsToRemove: ['lasso2d', 'select2d'],
        toImageButtonOptions: {
            format: 'svg',
            filename: 'deg_plot',
            height: 800,
            width: 1200,
            scale: 2
        },
        displaylogo: false
    };

    return (
        <div className="analysis-page">
            <div className="container">
                <div className="page-header">
                    <h1 className="page-title">
                        <span className="title-icon">🧬</span>
                        DEG Analysis
                    </h1>
                    <p className="page-description">
                        Upload your gene expression data for differential expression analysis with publication-quality visualizations
                    </p>
                </div>

                {!results ? (
                    <div className="analysis-form-container">
                        <form onSubmit={handleSubmit} className="analysis-form">
                            <div className="form-section">
                                <h3><span className="step-badge">1</span> Upload Data</h3>
                                <div {...getRootProps()} className={`dropzone ${isDragActive ? 'active' : ''} ${file ? 'has-file' : ''}`}>
                                    <input {...getInputProps()} />
                                    {file ? (
                                        <div className="file-info">
                                            <span className="file-icon">✅</span>
                                            <div>
                                                <span className="file-name">{file.name}</span>
                                                <span className="file-size">{(file.size / 1024).toFixed(1)} KB</span>
                                            </div>
                                        </div>
                                    ) : (
                                        <div className="dropzone-content">
                                            <span className="upload-icon">📂</span>
                                            <h4>Drag & drop your file here</h4>
                                            <p>or click to browse</p>
                                            <span className="file-types">Supports: CSV, TSV, TXT</span>
                                        </div>
                                    )}
                                </div>
                            </div>

                            <div className="form-section">
                                <h3><span className="step-badge">2</span> Configure Analysis</h3>
                                <div className="form-options">
                                    <div className="form-group">
                                        <label>Data Type</label>
                                        <div className="radio-group">
                                            <label className={`radio-option ${jobType === 'BULK' ? 'selected' : ''}`}>
                                                <input type="radio" value="BULK" checked={jobType === 'BULK'} onChange={(e) => setJobType(e.target.value)} />
                                                <div className="radio-content">
                                                    <span className="radio-icon">🔬</span>
                                                    <div>
                                                        <strong>Bulk RNA-seq</strong>
                                                        <small>Standard DEG analysis</small>
                                                    </div>
                                                </div>
                                            </label>
                                            <label className={`radio-option ${jobType === 'SCRNA' ? 'selected' : ''}`}>
                                                <input type="radio" value="SCRNA" checked={jobType === 'SCRNA'} onChange={(e) => setJobType(e.target.value)} />
                                                <div className="radio-content">
                                                    <span className="radio-icon">🧫</span>
                                                    <div>
                                                        <strong>Single-cell RNA-seq</strong>
                                                        <small>Scanpy clustering</small>
                                                    </div>
                                                </div>
                                            </label>
                                        </div>
                                    </div>

                                    <div className="form-group">
                                        <label>Data State</label>
                                        <div className="radio-group">
                                            <label className={`radio-option ${isNormalized ? 'selected' : ''}`}>
                                                <input type="radio" checked={isNormalized} onChange={() => setIsNormalized(true)} />
                                                <div className="radio-content">
                                                    <span className="radio-icon">📊</span>
                                                    <div>
                                                        <strong>Normalized</strong>
                                                        <small>Log2 or TPM/FPKM</small>
                                                    </div>
                                                </div>
                                            </label>
                                            <label className={`radio-option ${!isNormalized ? 'selected' : ''}`}>
                                                <input type="radio" checked={!isNormalized} onChange={() => setIsNormalized(false)} />
                                                <div className="radio-content">
                                                    <span className="radio-icon">🔢</span>
                                                    <div>
                                                        <strong>Raw Counts</strong>
                                                        <small>Uses pyDESeq2</small>
                                                    </div>
                                                </div>
                                            </label>
                                        </div>
                                    </div>
                                </div>
                            </div>

                            {error && <div className="error-message">⚠️ {error}</div>}

                            {/* Progress Bar */}
                            {loading && (
                                <div className="progress-container">
                                    <div className="progress-info">
                                        <span className="progress-step">Step {progress.step} of 3</span>
                                        <span className="progress-current">{progress.status}</span>
                                    </div>
                                    <div className="progress-bar">
                                        <div
                                            className="progress-fill"
                                            style={{ width: `${(progress.step / 3) * 100}%` }}
                                        ></div>
                                    </div>
                                    <div className="progress-message">{progress.message}</div>
                                </div>
                            )}

                            <button type="submit" className="btn btn-primary btn-lg submit-btn" disabled={loading}>
                                {loading ? (
                                    <>
                                        <span className="spinner"></span>
                                        <span>Analyzing...</span>
                                    </>
                                ) : (
                                    <>
                                        <span>🚀</span>
                                        <span>Run Analysis</span>
                                    </>
                                )}
                            </button>
                        </form>
                    </div>
                ) : (
                    <div className="results-container">
                        <div className="results-header">
                            <h2>Analysis Results</h2>
                            <div className="results-stats">
                                <span className="stat-badge up">↑ {results.filter(r => (r.log2FoldChange || 0) > 1 && (r.padj || r.pvalue || 1) < 0.05).length} Up</span>
                                <span className="stat-badge down">↓ {results.filter(r => (r.log2FoldChange || 0) < -1 && (r.padj || r.pvalue || 1) < 0.05).length} Down</span>
                                <span className="stat-badge total">Total: {results.length} genes</span>
                            </div>
                        </div>

                        <div className="results-tabs">
                            <button className={`tab-btn ${activeTab === 'volcano' ? 'active' : ''}`} onClick={() => setActiveTab('volcano')}>
                                <span>🌋</span> Volcano Plot
                            </button>
                            <button className={`tab-btn ${activeTab === 'heatmap' ? 'active' : ''}`} onClick={() => setActiveTab('heatmap')}>
                                <span>🔥</span> Heatmap
                            </button>
                            <button className={`tab-btn ${activeTab === 'table' ? 'active' : ''}`} onClick={() => setActiveTab('table')}>
                                <span>📋</span> Results Table
                            </button>
                        </div>

                        <div className="tab-content">
                            {activeTab === 'volcano' && (
                                <div className="plot-container">
                                    <Plot
                                        data={getVolcanoData()}
                                        layout={volcanoLayout}
                                        config={plotConfig}
                                        style={{ width: '100%', height: '600px' }}
                                    />
                                </div>
                            )}

                            {activeTab === 'heatmap' && (
                                <div className="heatmaps-grid">
                                    <div className="plot-container heatmap-up">
                                        <h4 className="heatmap-title up">🔺 Upregulated Genes</h4>
                                        <Plot
                                            data={getHeatmapData('up')}
                                            layout={getHeatmapLayout('Top 25 Upregulated')}
                                            config={plotConfig}
                                            style={{ width: '100%', height: '600px' }}
                                        />
                                    </div>
                                    <div className="plot-container heatmap-down">
                                        <h4 className="heatmap-title down">🔻 Downregulated Genes</h4>
                                        <Plot
                                            data={getHeatmapData('down')}
                                            layout={getHeatmapLayout('Top 25 Downregulated')}
                                            config={plotConfig}
                                            style={{ width: '100%', height: '600px' }}
                                        />
                                    </div>
                                </div>
                            )}

                            {activeTab === 'table' && (
                                <div className="tables-container">
                                    <div className="table-section">
                                        <div className="table-header">
                                            <h3>🔺 Top 10 Upregulated Genes</h3>
                                            <button className="btn btn-secondary" onClick={() => downloadCSV(getTopGenes('up'), 'upregulated_genes.csv')}>
                                                📥 Download CSV
                                            </button>
                                        </div>
                                        <div className="table-wrapper">
                                            <table className="results-table">
                                                <thead>
                                                    <tr>
                                                        <th>#</th>
                                                        <th>Gene</th>
                                                        <th>log₂FC</th>
                                                        <th>p-value</th>
                                                        <th>adj. p-value</th>
                                                    </tr>
                                                </thead>
                                                <tbody>
                                                    {getTopGenes('up').map((gene, i) => (
                                                        <tr key={i}>
                                                            <td>{i + 1}</td>
                                                            <td className="gene-name">{gene.gene}</td>
                                                            <td className="positive">+{(gene.log2FoldChange || 0).toFixed(3)}</td>
                                                            <td>{(gene.pvalue || 0).toExponential(2)}</td>
                                                            <td>{(gene.padj || gene.pvalue || 0).toExponential(2)}</td>
                                                        </tr>
                                                    ))}
                                                </tbody>
                                            </table>
                                        </div>
                                    </div>

                                    <div className="table-section">
                                        <div className="table-header">
                                            <h3>🔻 Top 10 Downregulated Genes</h3>
                                            <button className="btn btn-secondary" onClick={() => downloadCSV(getTopGenes('down'), 'downregulated_genes.csv')}>
                                                📥 Download CSV
                                            </button>
                                        </div>
                                        <div className="table-wrapper">
                                            <table className="results-table">
                                                <thead>
                                                    <tr>
                                                        <th>#</th>
                                                        <th>Gene</th>
                                                        <th>log₂FC</th>
                                                        <th>p-value</th>
                                                        <th>adj. p-value</th>
                                                    </tr>
                                                </thead>
                                                <tbody>
                                                    {getTopGenes('down').map((gene, i) => (
                                                        <tr key={i}>
                                                            <td>{i + 1}</td>
                                                            <td className="gene-name">{gene.gene}</td>
                                                            <td className="negative">{(gene.log2FoldChange || 0).toFixed(3)}</td>
                                                            <td>{(gene.pvalue || 0).toExponential(2)}</td>
                                                            <td>{(gene.padj || gene.pvalue || 0).toExponential(2)}</td>
                                                        </tr>
                                                    ))}
                                                </tbody>
                                            </table>
                                        </div>
                                    </div>

                                    <div className="full-export">
                                        <button className="btn btn-primary btn-lg" onClick={() => downloadCSV(results, 'all_deg_results.csv')}>
                                            📦 Download Full Results (CSV)
                                        </button>
                                    </div>
                                </div>
                            )}
                        </div>

                        <div className="results-actions">
                            <button className="btn btn-secondary" onClick={() => { setResults(null); setFile(null); setJob(null); }}>
                                ↩️ New Analysis
                            </button>
                        </div>
                    </div>
                )}
            </div>
        </div>
    );
};

export default Analysis;
