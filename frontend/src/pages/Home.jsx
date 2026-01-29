import React from 'react';
import { Link } from 'react-router-dom';
import './Home.css';

const Home = () => {
    return (
        <div className="home">
            {/* Floating Particles */}
            <div className="particles-bg">
                {[...Array(30)].map((_, i) => (
                    <div
                        key={i}
                        className="particle"
                        style={{
                            left: `${Math.random() * 100}%`,
                            top: `${Math.random() * 100}%`,
                            animationDelay: `${Math.random() * 20}s`,
                            animationDuration: `${15 + Math.random() * 10}s`
                        }}
                    />
                ))}
            </div>

            {/* Hero Section */}
            <section className="hero">
                <div className="hero-container">
                    <div className="hero-content">
                        <div className="hero-badge">
                            <span className="badge-icon">🧬</span>
                            <span>Next-Gen Bioinformatics Platform</span>
                        </div>

                        <h1 className="hero-title">
                            <span className="title-line">From Genes to</span>
                            <span className="title-line gradient-text">Therapeutic Insights</span>
                        </h1>

                        <p className="hero-description">
                            A powerful platform for differential gene expression analysis,
                            functional annotation, and drug-disease network discovery.
                            Transform your transcriptomics data into actionable insights.
                        </p>

                        <div className="hero-actions">
                            <Link to="/analysis" className="btn btn-primary btn-lg glow">
                                <span>Start Analysis</span>
                                <span className="btn-arrow">→</span>
                            </Link>
                            <Link to="/annotation" className="btn btn-secondary btn-lg">
                                <span>Explore Features</span>
                            </Link>
                        </div>

                        <div className="hero-stats">
                            <div className="stat">
                                <span className="stat-value">pyDESeq2</span>
                                <span className="stat-label">DEG Analysis</span>
                            </div>
                            <div className="stat-divider"></div>
                            <div className="stat">
                                <span className="stat-value">KEGG + GO</span>
                                <span className="stat-label">Enrichment</span>
                            </div>
                            <div className="stat-divider"></div>
                            <div className="stat">
                                <span className="stat-value">Open Targets</span>
                                <span className="stat-label">Disease Links</span>
                            </div>
                        </div>
                    </div>

                    <div className="hero-visual">
                        <div className="dna-orb">
                            <div className="orb-glow"></div>
                            <div className="dna-helix">
                                {[...Array(20)].map((_, i) => (
                                    <div key={i} className="dna-pair" style={{ '--i': i }}>
                                        <div className="dna-nucleotide left"></div>
                                        <div className="dna-backbone"></div>
                                        <div className="dna-nucleotide right"></div>
                                    </div>
                                ))}
                            </div>
                        </div>
                    </div>
                </div>
            </section>

            {/* Features Section */}
            <section className="features section">
                <div className="container">
                    <div className="section-header">
                        <span className="section-label">What We Offer</span>
                        <h2 className="section-title">
                            Comprehensive Analysis <span className="gradient-text">Pipeline</span>
                        </h2>
                        <p className="section-description">
                            End-to-end solution for transcriptomics research and drug discovery
                        </p>
                    </div>

                    <div className="features-grid stagger">
                        <div className="feature-card card">
                            <div className="feature-icon-wrapper">
                                <span className="feature-icon">📊</span>
                            </div>
                            <h3>DEG Analysis</h3>
                            <p>Compute differentially expressed genes with pyDESeq2. Support for bulk and single-cell RNA-seq data.</p>
                        </div>

                        <div className="feature-card card">
                            <div className="feature-icon-wrapper">
                                <span className="feature-icon">🛤️</span>
                            </div>
                            <h3>Pathway Enrichment</h3>
                            <p>KEGG pathway analysis with pathway names and Gene Ontology annotation via UniProt.</p>
                        </div>

                        <div className="feature-card card">
                            <div className="feature-icon-wrapper">
                                <span className="feature-icon">🏥</span>
                            </div>
                            <h3>Disease Mapping</h3>
                            <p>Connect genes to diseases using Open Targets Platform with association scores.</p>
                        </div>

                        <div className="feature-card card">
                            <div className="feature-icon-wrapper">
                                <span className="feature-icon">💊</span>
                            </div>
                            <h3>Drug Discovery</h3>
                            <p>Find drug-gene interactions via DGIdb and Open Targets for therapeutic insights.</p>
                        </div>

                        <div className="feature-card card">
                            <div className="feature-icon-wrapper">
                                <span className="feature-icon">🕸️</span>
                            </div>
                            <h3>Network Visualization</h3>
                            <p>Interactive gene-disease-drug networks with Cytoscape.js for hypothesis generation.</p>
                        </div>

                        <div className="feature-card card">
                            <div className="feature-icon-wrapper">
                                <span className="feature-icon">📈</span>
                            </div>
                            <h3>Publication Plots</h3>
                            <p>High-quality volcano plots, heatmaps, and export to SVG for publications.</p>
                        </div>
                    </div>
                </div>
            </section>

            {/* CTA Section */}
            <section className="cta-section section">
                <div className="container">
                    <div className="cta-card card">
                        <div className="cta-content">
                            <h2>Ready to accelerate your research?</h2>
                            <p>Start analyzing your gene expression data today</p>
                        </div>
                        <Link to="/analysis" className="btn btn-primary btn-lg">
                            Launch Analysis
                        </Link>
                    </div>
                </div>
            </section>
        </div>
    );
};

export default Home;
