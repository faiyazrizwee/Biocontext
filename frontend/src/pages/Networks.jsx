import React, { useState, useEffect, useRef, useCallback } from 'react';
import cytoscape from 'cytoscape';
import { useDropzone } from 'react-dropzone';
import * as XLSX from 'xlsx';
import { networksApi, annotationApi } from '../api';
import './Networks.css';

const Networks = () => {
    const [geneInput, setGeneInput] = useState('');
    const [uploadedFile, setUploadedFile] = useState(null);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState(null);
    const [networkData, setNetworkData] = useState(null);
    const [networkType, setNetworkType] = useState('therapy');
    const cyRef = useRef(null);
    const containerRef = useRef(null);

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

    const buildNetwork = async () => {
        const genes = parseGenes();
        if (genes.length === 0) {
            setError('Please enter gene symbols');
            return;
        }

        setLoading(true);
        setError(null);

        try {
            // Fetch annotation data
            const [diseasesRes, drugsRes] = await Promise.all([
                annotationApi.getDiseases(genes),
                annotationApi.getDrugs(genes)
            ]);

            // Build network
            const data = await networksApi.getTherapyNetwork({
                genes,
                diseases: diseasesRes.associations || [],
                drugs: drugsRes.drugs || []
            });

            setNetworkData(data.elements);
        } catch (err) {
            setError('Failed to build network');
        } finally {
            setLoading(false);
        }
    };

    useEffect(() => {
        if (networkData && containerRef.current) {
            // Destroy previous instance
            if (cyRef.current) {
                cyRef.current.destroy();
            }

            // Create Cytoscape instance
            cyRef.current = cytoscape({
                container: containerRef.current,
                elements: [...networkData.nodes, ...networkData.edges],
                style: [
                    {
                        selector: 'node',
                        style: {
                            'label': 'data(label)',
                            'text-valign': 'center',
                            'text-halign': 'center',
                            'font-size': '10px',
                            'color': '#fff',
                            'text-outline-color': '#000',
                            'text-outline-width': 1,
                            'width': 40,
                            'height': 40
                        }
                    },
                    {
                        selector: 'node[type="gene"]',
                        style: {
                            'background-color': '#0066cc',
                            'shape': 'ellipse'
                        }
                    },
                    {
                        selector: 'node[type="disease"]',
                        style: {
                            'background-color': '#e0245e',
                            'shape': 'round-rectangle'
                        }
                    },
                    {
                        selector: 'node[type="drug"]',
                        style: {
                            'background-color': '#17bf63',
                            'shape': 'diamond'
                        }
                    },
                    {
                        selector: 'node[type="pathway"]',
                        style: {
                            'background-color': '#ffad1f',
                            'shape': 'hexagon'
                        }
                    },
                    {
                        selector: 'edge',
                        style: {
                            'width': 2,
                            'line-color': '#38444d',
                            'curve-style': 'bezier',
                            'target-arrow-shape': 'triangle',
                            'target-arrow-color': '#38444d'
                        }
                    }
                ],
                layout: {
                    name: 'cose',
                    animate: true,
                    animationDuration: 500,
                    fit: true,
                    padding: 50,
                    nodeRepulsion: 8000,
                    idealEdgeLength: 100
                }
            });

            // Add hover tooltips
            cyRef.current.on('mouseover', 'node', (e) => {
                e.target.style('border-width', 3);
                e.target.style('border-color', '#fff');
            });

            cyRef.current.on('mouseout', 'node', (e) => {
                e.target.style('border-width', 0);
            });
        }

        return () => {
            if (cyRef.current) {
                cyRef.current.destroy();
            }
        };
    }, [networkData]);

    const exportNetwork = (format) => {
        if (!cyRef.current) return;

        if (format === 'png') {
            const png = cyRef.current.png({ scale: 2, bg: '#15202b' });
            const link = document.createElement('a');
            link.href = png;
            link.download = 'network.png';
            link.click();
        } else if (format === 'json') {
            const json = JSON.stringify(networkData, null, 2);
            const blob = new Blob([json], { type: 'application/json' });
            const url = URL.createObjectURL(blob);
            const link = document.createElement('a');
            link.href = url;
            link.download = 'network_data.json';
            link.click();
            URL.revokeObjectURL(url);
        }
    };

    const resetLayout = () => {
        if (cyRef.current) {
            cyRef.current.layout({
                name: 'cose',
                animate: true,
                animationDuration: 500
            }).run();
        }
    };

    const fitNetwork = () => {
        if (cyRef.current) {
            cyRef.current.fit(50);
        }
    };

    return (
        <div className="networks-page">
            <div className="container">
                <div className="page-header">
                    <h1 className="page-title">
                        <span className="title-icon">🕸️</span>
                        Network Visualization
                    </h1>
                    <p className="page-description">
                        Explore gene-disease-drug relationships through interactive network graphs
                    </p>
                </div>

                <div className="networks-layout">
                    {/* Controls Panel */}
                    <div className="controls-panel card">
                        <h3>Build Network</h3>

                        <div {...getRootProps()} className={`file-dropzone ${isDragActive ? 'active' : ''} ${uploadedFile ? 'has-file' : ''}`}>
                            <input {...getInputProps()} />
                            <div className="dropzone-content">
                                {uploadedFile ? (
                                    <div className="file-info">
                                        <span className="file-icon">📄</span>
                                        <div className="file-details">
                                            <span className="file-name">{uploadedFile.name}</span>
                                            <span className="file-size">{(uploadedFile.size / 1024).toFixed(1)} KB</span>
                                        </div>
                                        <button className="clear-btn" onClick={(e) => { e.stopPropagation(); clearFile(); }}>✕</button>
                                    </div>
                                ) : (
                                    <>
                                        <span className="upload-icon">📁</span>
                                        <p>{isDragActive ? 'Drop file here' : 'Click or drag file to upload'}</p>
                                        <span className="file-formats">TXT, CSV, TSV, XLSX</span>
                                    </>
                                )}
                            </div>
                        </div>

                        <div className="divider">
                            <span>OR</span>
                        </div>

                        <textarea
                            className="gene-input"
                            placeholder="Enter gene symbols&#10;(one per line or comma-separated)"
                            value={geneInput}
                            onChange={(e) => setGeneInput(e.target.value)}
                        />
                        <div className="gene-count">
                            <span className="count-badge">{parseGenes().length}</span> genes detected
                        </div>

                        {error && <div className="error-message">{error}</div>}

                        <button
                            className="btn btn-primary build-btn"
                            onClick={buildNetwork}
                            disabled={loading}
                        >
                            {loading ? 'Building...' : 'Build Network'}
                        </button>

                        {networkData && (
                            <>
                                <div className="legend">
                                    <h4>Legend</h4>
                                    <div className="legend-item">
                                        <span className="legend-shape gene"></span>
                                        <span>Gene</span>
                                    </div>
                                    <div className="legend-item">
                                        <span className="legend-shape disease"></span>
                                        <span>Disease</span>
                                    </div>
                                    <div className="legend-item">
                                        <span className="legend-shape drug"></span>
                                        <span>Drug</span>
                                    </div>
                                </div>

                                <div className="network-controls">
                                    <h4>Controls</h4>
                                    <button className="btn btn-secondary" onClick={resetLayout}>
                                        Reset Layout
                                    </button>
                                    <button className="btn btn-secondary" onClick={fitNetwork}>
                                        Fit to View
                                    </button>
                                </div>

                                <div className="export-controls">
                                    <h4>Export</h4>
                                    <button className="btn btn-secondary" onClick={() => exportNetwork('png')}>
                                        Download PNG
                                    </button>
                                    <button className="btn btn-secondary" onClick={() => exportNetwork('json')}>
                                        Download JSON
                                    </button>
                                </div>
                            </>
                        )}
                    </div>

                    {/* Network Container */}
                    <div className="network-container card">
                        {networkData ? (
                            <div ref={containerRef} className="cytoscape-container"></div>
                        ) : (
                            <div className="empty-state">
                                <span className="empty-icon">🕸️</span>
                                <p>Enter genes and build a network to visualize relationships</p>
                            </div>
                        )}
                    </div>
                </div>

                {/* Stats */}
                {networkData && (
                    <div className="network-stats">
                        <div className="stat-card card">
                            <div className="stat-value">{networkData.nodes?.length || 0}</div>
                            <div className="stat-label">Nodes</div>
                        </div>
                        <div className="stat-card card">
                            <div className="stat-value">{networkData.edges?.length || 0}</div>
                            <div className="stat-label">Edges</div>
                        </div>
                        <div className="stat-card card">
                            <div className="stat-value">
                                {networkData.nodes?.filter(n => n.data.type === 'gene').length || 0}
                            </div>
                            <div className="stat-label">Genes</div>
                        </div>
                        <div className="stat-card card">
                            <div className="stat-value">
                                {networkData.nodes?.filter(n => n.data.type === 'disease').length || 0}
                            </div>
                            <div className="stat-label">Diseases</div>
                        </div>
                        <div className="stat-card card">
                            <div className="stat-value">
                                {networkData.nodes?.filter(n => n.data.type === 'drug').length || 0}
                            </div>
                            <div className="stat-label">Drugs</div>
                        </div>
                    </div>
                )}
            </div>
        </div>
    );
};

export default Networks;
