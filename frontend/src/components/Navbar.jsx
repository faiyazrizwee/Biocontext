import React, { useState, useEffect } from 'react';
import { Link, useLocation } from 'react-router-dom';
import './Navbar.css';

const Navbar = () => {
    const [scrolled, setScrolled] = useState(false);
    const [menuOpen, setMenuOpen] = useState(false);
    const location = useLocation();

    useEffect(() => {
        const handleScroll = () => {
            setScrolled(window.scrollY > 20);
        };
        window.addEventListener('scroll', handleScroll);
        return () => window.removeEventListener('scroll', handleScroll);
    }, []);

    const navLinks = [
        { path: '/', label: 'Home' },
        { path: '/analysis', label: 'Analysis' },
        { path: '/annotation', label: 'Annotation' },
        { path: '/networks', label: 'Networks' }
    ];

    return (
        <nav className={`navbar ${scrolled ? 'scrolled' : ''}`}>
            <div className="navbar-container">
                <Link to="/" className="logo">
                    <span className="logo-icon">🧬</span>
                    <span className="logo-text">
                        <span className="logo-primary">Gene</span>
                        <span className="logo-accent">2Therapy</span>
                    </span>
                </Link>

                <div className={`nav-menu ${menuOpen ? 'open' : ''}`}>
                    {navLinks.map(link => (
                        <Link
                            key={link.path}
                            to={link.path}
                            className={`nav-link ${location.pathname === link.path ? 'active' : ''}`}
                            onClick={() => setMenuOpen(false)}
                        >
                            <span className="nav-link-text">{link.label}</span>
                            <span className="nav-link-indicator"></span>
                        </Link>
                    ))}
                </div>

                <div className="nav-actions">
                    <Link to="/analysis" className="btn btn-primary nav-cta">
                        Get Started
                    </Link>
                </div>

                <button
                    className={`menu-toggle ${menuOpen ? 'open' : ''}`}
                    onClick={() => setMenuOpen(!menuOpen)}
                    aria-label="Toggle menu"
                >
                    <span></span>
                    <span></span>
                    <span></span>
                </button>
            </div>
        </nav>
    );
};

export default Navbar;
