import React, { useState, useCallback } from 'react';
import { useDropzone } from 'react-dropzone';
import * as XLSX from 'xlsx';
import { annotationApi } from '../api';
import './Annotation.css';

const Annotation = () => {
    const [geneInput, setGeneInput] = useState('');
    const [uploadedFile, setUploadedFile] = useState(null);
    const [loading, setLoading] = useState({});
    const [progress, setProgress] = useState({ step: 0, total: 5, current: '' });
    const [results, setResults] = useState({
        geneMap: null,
        kegg: null,
        go: null,
        diseases: null,
        drugs: null
    });
    const [error, setError] = useState(null);
    const [activeSection, setActiveSection] = useState(null);

    // Parse genes from text input
    const parseGenes = () => {
        return geneInput
            .split(/[\n,\s]+/)
            .map(g => g.trim().toUpperCase())
            .filter(g => g.length > 0);
    };

    // Handle file upload
    const onDrop = useCallback(async (acceptedFiles) => {
        if (acceptedFiles.length === 0) return;

        const file = acceptedFiles[0];
        setUploadedFile(file);
        setError(null);

        try {
            const genes = await parseFile(file);
            if (genes.length > 0) {
                setGeneInput(genes.join('\n'));
            } else {
                setError('No valid gene symbols found in file');
            }
        } catch (err) {
            setError(`Error reading file: ${err.message}`);
        }
    }, []);

    // Parse different file formats
    const parseFile = (file) => {
        return new Promise((resolve, reject) => {
            const fileName = file.name.toLowerCase();
            const reader = new FileReader();

            if (fileName.endsWith('.xlsx') || fileName.endsWith('.xls')) {
                // Excel file
                reader.onload = (e) => {
                    try {
                        const data = new Uint8Array(e.target.result);
                        const workbook = XLSX.read(data, { type: 'array' });
                        const firstSheet = workbook.Sheets[workbook.SheetNames[0]];
                        const jsonData = XLSX.utils.sheet_to_json(firstSheet, { header: 1 });

                        // Extract genes from first column (skip header if present)
                        const genes = [];
                        jsonData.forEach((row, idx) => {
                            if (row[0]) {
                                const gene = String(row[0]).trim().toUpperCase();
                                // Skip common header names
                                if (idx === 0 && ['GENE', 'GENES', 'SYMBOL', 'GENE_SYMBOL', 'GENE_NAME', 'ID'].includes(gene)) {
                                    return;
                                }
                                if (gene.length > 0 && gene.length < 20) {
                                    genes.push(gene);
                                }
                            }
                        });
                        resolve(genes);
                    } catch (err) {
                        reject(err);
                    }
                };
                reader.readAsArrayBuffer(file);
            } else {
                // Text files (txt, csv, tsv)
                reader.onload = (e) => {
                    try {
                        const text = e.target.result;
                        let genes = [];

                        if (fileName.endsWith('.csv')) {
                            // CSV: split by comma or newline
                            genes = text.split(/[,\n\r]+/);
                        } else if (fileName.endsWith('.tsv')) {
                            // TSV: split by tab or newline
                            genes = text.split(/[\t\n\r]+/);
                        } else {
                            // TXT: split by any whitespace, comma, or newline
                            genes = text.split(/[\s,\n\r\t]+/);
                        }

                        // Clean and filter genes
                        genes = genes
                            .map(g => g.trim().toUpperCase())
                            .filter(g => {
                                // Skip common headers
                                if (['GENE', 'GENES', 'SYMBOL', 'GENE_SYMBOL', 'GENE_NAME', 'ID', ''].includes(g)) {
                                    return false;
                                }
                                return g.length > 0 && g.length < 20;
                            });

                        resolve(genes);
                    } catch (err) {
                        reject(err);
                    }
                };
                reader.readAsText(file);
            }
        });
    };

    const { getRootProps, getInputProps, isDragActive } = useDropzone({
        onDrop,
        accept: {
            'text/plain': ['.txt'],
            'text/csv': ['.csv'],
            'text/tab-separated-values': ['.tsv'],
            'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet': ['.xlsx'],
            'application/vnd.ms-excel': ['.xls']
        },
        maxFiles: 1
    });

    const clearFile = () => {
        setUploadedFile(null);
        setGeneInput('');
    };

    const runAnalysis = async (type) => {
        const genes = parseGenes();
        if (genes.length === 0) {
            setError('Please enter gene symbols');
            return;
        }

        setLoading(prev => ({ ...prev, [type]: true }));
        setError(null);

        try {
            let data;
            switch (type) {
                case 'geneMap':
                    data = await annotationApi.mapGenes(genes);
                    setResults(prev => ({ ...prev, geneMap: data.results }));
                    break;
                case 'kegg':
                    data = await annotationApi.getKEGG(genes);
                    setResults(prev => ({ ...prev, kegg: data.pathways }));
                    break;
                case 'go':
                    data = await annotationApi.getGO(genes);
                    setResults(prev => ({ ...prev, go: data.go_terms }));
                    break;
                case 'diseases':
                    data = await annotationApi.getDiseases(genes);
                    setResults(prev => ({ ...prev, diseases: data.associations }));
                    break;
                case 'drugs':
                    data = await annotationApi.getDrugs(genes);
                    setResults(prev => ({ ...prev, drugs: data.drugs }));
                    break;
            }
            setActiveSection(type);
        } catch (err) {
            setError(err.response?.data?.error || `Failed to fetch ${type} data`);
        } finally {
            setLoading(prev => ({ ...prev, [type]: false }));
        }
    };

    const runAll = async () => {
        const genes = parseGenes();
        if (genes.length === 0) {
            setError('Please enter gene symbols');
            return;
        }

        setLoading({ all: true });
        setError(null);
        setProgress({ step: 0, total: 5, current: 'Starting...' });

        const steps = [
            { name: 'NCBI Gene Mapping', key: 'geneMap', fn: () => annotationApi.mapGenes(genes) },
            { name: 'KEGG Pathways', key: 'kegg', fn: () => annotationApi.getKEGG(genes) },
            { name: 'Gene Ontology', key: 'go', fn: () => annotationApi.getGO(genes) },
            { name: 'Disease Associations', key: 'diseases', fn: () => annotationApi.getDiseases(genes) },
            { name: 'Drug Interactions', key: 'drugs', fn: () => annotationApi.getDrugs(genes) }
        ];

        const newResults = { ...results };

        try {
            for (let i = 0; i < steps.length; i++) {
                setProgress({ step: i + 1, total: steps.length, current: steps[i].name });

                try {
                    const data = await steps[i].fn();
                    switch (steps[i].key) {
                        case 'geneMap': newResults.geneMap = data.results; break;
                        case 'kegg': newResults.kegg = data.pathways; break;
                        case 'go': newResults.go = data.go_terms; break;
                        case 'diseases': newResults.diseases = data.associations; break;
                        case 'drugs': newResults.drugs = data.drugs; break;
                    }
                } catch (stepErr) {
                    console.error(`Error in ${steps[i].name}:`, stepErr);
                    // Continue with other steps even if one fails
                }
            }

            setResults(newResults);
            setActiveSection('geneMap');
        } catch (err) {
            setError('Failed to complete all analyses');
        } finally {
            setLoading({ all: false });
            setProgress({ step: 0, total: 5, current: '' });
        }
    };

    return (
        <div className="annotation-page">
            <div className="container">
                <div className="page-header">
                    <h1 className="page-title">
                        <span className="title-icon">🧬</span>
                        Functional Annotation
                    </h1>
                    <p className="page-description">
                        Map genes to pathways, diseases, and drug targets using multiple databases
                    </p>
                </div>

                <div className="annotation-layout">
                    {/* Input Panel */}
                    <div className="input-panel card">
                        <h3>📁 Upload Gene List</h3>
                        <div
                            {...getRootProps()}
                            className={`file-dropzone ${isDragActive ? 'active' : ''} ${uploadedFile ? 'has-file' : ''}`}
                        >
                            <input {...getInputProps()} />
                            {uploadedFile ? (
                                <div className="file-info">
                                    <span className="file-icon">✅</span>
                                    <div className="file-details">
                                        <span className="file-name">{uploadedFile.name}</span>
                                        <span className="file-size">{(uploadedFile.size / 1024).toFixed(1)} KB</span>
                                    </div>
                                    <button className="clear-btn" onClick={(e) => { e.stopPropagation(); clearFile(); }}>✕</button>
                                </div>
                            ) : (
                                <div className="dropzone-content">
                                    <span className="upload-icon">📤</span>
                                    <p><strong>Drop file here</strong> or click to browse</p>
                                    <span className="file-formats">TXT, CSV, TSV, XLSX</span>
                                </div>
                            )}
                        </div>

                        <div className="divider">
                            <span>or enter manually</span>
                        </div>

                        <h3>✏️ Gene List</h3>
                        <textarea
                            className="gene-input"
                            placeholder="Enter gene symbols (one per line or comma-separated)&#10;&#10;Example:&#10;TP53&#10;BRCA1&#10;EGFR&#10;MYC"
                            value={geneInput}
                            onChange={(e) => setGeneInput(e.target.value)}
                        />
                        <div className="gene-count">
                            <span className="count-badge">{parseGenes().length}</span> genes detected
                        </div>

                        {error && <div className="error-message">⚠️ {error}</div>}

                        {/* Progress Bar */}
                        {loading.all && (
                            <div className="progress-container">
                                <div className="progress-info">
                                    <span className="progress-step">Step {progress.step} of {progress.total}</span>
                                    <span className="progress-current">{progress.current}</span>
                                </div>
                                <div className="progress-bar">
                                    <div
                                        className="progress-fill"
                                        style={{ width: `${(progress.step / progress.total) * 100}%` }}
                                    ></div>
                                </div>
                            </div>
                        )}

                        <div className="analysis-buttons">
                            <button
                                className="btn btn-primary btn-lg"
                                onClick={runAll}
                                disabled={loading.all}
                            >
                                {loading.all ? '⏳ Running All...' : '🚀 Run All Analyses'}
                            </button>

                            <div className="individual-buttons">
                                <button
                                    className="btn btn-secondary"
                                    onClick={() => runAnalysis('geneMap')}
                                    disabled={loading.geneMap}
                                >
                                    🧬 NCBI
                                </button>
                                <button
                                    className="btn btn-secondary"
                                    onClick={() => runAnalysis('kegg')}
                                    disabled={loading.kegg}
                                >
                                    🛤️ KEGG
                                </button>
                                <button
                                    className="btn btn-secondary"
                                    onClick={() => runAnalysis('go')}
                                    disabled={loading.go}
                                >
                                    🔬 GO
                                </button>
                                <button
                                    className="btn btn-secondary"
                                    onClick={() => runAnalysis('diseases')}
                                    disabled={loading.diseases}
                                >
                                    🏥 Diseases
                                </button>
                                <button
                                    className="btn btn-secondary"
                                    onClick={() => runAnalysis('drugs')}
                                    disabled={loading.drugs}
                                >
                                    💊 Drugs
                                </button>
                            </div>
                        </div>
                    </div>

                    {/* Results Panel */}
                    <div className="results-panel">
                        {/* Section Tabs */}
                        {(results.geneMap || results.kegg || results.go || results.diseases || results.drugs) && (
                            <div className="section-tabs">
                                {results.geneMap && (
                                    <button
                                        className={`section-tab ${activeSection === 'geneMap' ? 'active' : ''}`}
                                        onClick={() => setActiveSection('geneMap')}
                                    >
                                        🧬 NCBI ({results.geneMap.length})
                                    </button>
                                )}
                                {results.kegg && (
                                    <button
                                        className={`section-tab ${activeSection === 'kegg' ? 'active' : ''}`}
                                        onClick={() => setActiveSection('kegg')}
                                    >
                                        🛤️ KEGG ({results.kegg.length})
                                    </button>
                                )}
                                {results.go && (
                                    <button
                                        className={`section-tab ${activeSection === 'go' ? 'active' : ''}`}
                                        onClick={() => setActiveSection('go')}
                                    >
                                        🔬 GO ({results.go.length})
                                    </button>
                                )}
                                {results.diseases && (
                                    <button
                                        className={`section-tab ${activeSection === 'diseases' ? 'active' : ''}`}
                                        onClick={() => setActiveSection('diseases')}
                                    >
                                        🏥 Diseases ({results.diseases.length})
                                    </button>
                                )}
                                {results.drugs && (
                                    <button
                                        className={`section-tab ${activeSection === 'drugs' ? 'active' : ''}`}
                                        onClick={() => setActiveSection('drugs')}
                                    >
                                        💊 Drugs ({results.drugs.length})
                                    </button>
                                )}
                            </div>
                        )}

                        {/* Results Content */}
                        <div className="results-content">
                            {activeSection === 'geneMap' && results.geneMap && (
                                <div className="result-section card">
                                    <h3>🧬 NCBI Gene Mapping</h3>
                                    <div className="table-wrapper">
                                        <table className="results-table">
                                            <thead>
                                                <tr>
                                                    <th>Symbol</th>
                                                    <th>NCBI Gene IDs</th>
                                                    <th>Status</th>
                                                </tr>
                                            </thead>
                                            <tbody>
                                                {results.geneMap.map((item, i) => (
                                                    <tr key={i}>
                                                        <td className="gene-name">{item.symbol}</td>
                                                        <td>{item.ncbi_gene_ids?.join(', ') || '-'}</td>
                                                        <td className={item.status === 'found' ? 'positive' : 'negative'}>
                                                            {item.status}
                                                        </td>
                                                    </tr>
                                                ))}
                                            </tbody>
                                        </table>
                                    </div>
                                </div>
                            )}

                            {activeSection === 'kegg' && results.kegg && (
                                <div className="result-section card">
                                    <h3>🛤️ KEGG Pathways</h3>
                                    <div className="table-wrapper">
                                        <table className="results-table">
                                            <thead>
                                                <tr>
                                                    <th>Pathway</th>
                                                    <th>Gene Count</th>
                                                    <th>Genes</th>
                                                </tr>
                                            </thead>
                                            <tbody>
                                                {results.kegg.map((item, i) => (
                                                    <tr key={i}>
                                                        <td>{item.pathway}</td>
                                                        <td className="count-cell">{item.count}</td>
                                                        <td className="genes-cell">{item.genes?.join(', ')}</td>
                                                    </tr>
                                                ))}
                                            </tbody>
                                        </table>
                                    </div>
                                </div>
                            )}

                            {activeSection === 'go' && results.go && (
                                <div className="result-section card">
                                    <h3>🔬 Gene Ontology</h3>
                                    <div className="table-wrapper">
                                        <table className="results-table">
                                            <thead>
                                                <tr>
                                                    <th>Gene</th>
                                                    <th>GO ID</th>
                                                    <th>Term</th>
                                                </tr>
                                            </thead>
                                            <tbody>
                                                {results.go.map((item, i) => (
                                                    <tr key={i}>
                                                        <td className="gene-name">{item.gene}</td>
                                                        <td className="go-id">{item.go_id}</td>
                                                        <td>{item.term}</td>
                                                    </tr>
                                                ))}
                                            </tbody>
                                        </table>
                                    </div>
                                </div>
                            )}

                            {activeSection === 'diseases' && results.diseases && (
                                <div className="result-section card">
                                    <h3>🏥 Disease Associations</h3>
                                    <div className="table-wrapper">
                                        <table className="results-table">
                                            <thead>
                                                <tr>
                                                    <th>Gene</th>
                                                    <th>Disease</th>
                                                    <th>Score</th>
                                                </tr>
                                            </thead>
                                            <tbody>
                                                {results.diseases.map((item, i) => (
                                                    <tr key={i}>
                                                        <td className="gene-name">{item.gene}</td>
                                                        <td>{item.disease}</td>
                                                        <td className="score-cell">{item.score?.toFixed(3) || '-'}</td>
                                                    </tr>
                                                ))}
                                            </tbody>
                                        </table>
                                    </div>
                                </div>
                            )}

                            {activeSection === 'drugs' && results.drugs && (
                                <div className="result-section card">
                                    <h3>💊 Drug Interactions</h3>
                                    <div className="table-wrapper">
                                        <table className="results-table">
                                            <thead>
                                                <tr>
                                                    <th>Gene</th>
                                                    <th>Drug</th>
                                                    <th>Interaction Type</th>
                                                </tr>
                                            </thead>
                                            <tbody>
                                                {results.drugs.map((item, i) => (
                                                    <tr key={i}>
                                                        <td className="gene-name">{item.gene}</td>
                                                        <td className="drug-name">{item.drug_name}</td>
                                                        <td>{item.interaction_types?.join(', ') || '-'}</td>
                                                    </tr>
                                                ))}
                                            </tbody>
                                        </table>
                                    </div>
                                </div>
                            )}

                            {!activeSection && (
                                <div className="empty-state">
                                    <span className="empty-icon">🔍</span>
                                    <h3>Ready to Analyze</h3>
                                    <p>Upload a file or enter gene symbols, then run an analysis to see results</p>
                                </div>
                            )}
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default Annotation;
