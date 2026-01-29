import React from 'react';
import './Footer.css';

const Footer = () => {
    return (
        <footer className="footer">
            <div className="footer-container">
                <div className="footer-brand">
                    <div className="footer-logo">
                        <span className="logo-icon">🧬</span>
                        <span className="logo-text">
                            <span className="logo-primary">Gene</span>
                            <span className="logo-accent">2Therapy</span>
                        </span>
                    </div>
                    <p className="footer-tagline">
                        From differential expression to therapeutic insights
                    </p>
                </div>

                <div className="footer-links">
                    <div className="footer-section">
                        <h4>Platform</h4>
                        <ul>
                            <li><a href="/analysis">DEG Analysis</a></li>
                            <li><a href="/annotation">Annotation</a></li>
                            <li><a href="/networks">Networks</a></li>
                        </ul>
                    </div>

                    <div className="footer-section">
                        <h4>Resources</h4>
                        <ul>
                            <li><a href="https://www.kegg.jp/" target="_blank" rel="noopener noreferrer">KEGG</a></li>
                            <li><a href="https://platform.opentargets.org/" target="_blank" rel="noopener noreferrer">Open Targets</a></li>
                            <li><a href="https://dgidb.org/" target="_blank" rel="noopener noreferrer">DGIdb</a></li>
                        </ul>
                    </div>
                </div>
            </div>

            <div className="footer-bottom">
                <div className="footer-container">
                    <p>&copy; 2025 Gene2Therapy. Built for Research.</p>
                    <p className="footer-tech">Django · React · Plotly · Cytoscape</p>
                </div>
            </div>
        </footer>
    );
};

export default Footer;
