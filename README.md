# Gene2Therapy

A modern, scalable bioinformatics web application for differential gene expression (DEG) analysis and therapeutic target discovery.

## 🧬 Features

### Module 1: DEG Analysis
- **Bulk RNA-seq** support with statistical DEG computation
- **Single-cell RNA-seq** support with Scanpy integration
- Support for both **normalized** and **raw count** data (pyDESeq2)
- Interactive **Volcano plots** and **Heatmaps** (Plotly.js)
- Download results as CSV and high-resolution images (PNG/SVG)

### Module 2: Functional Annotation
- **NCBI Gene** mapping
- **KEGG** pathway enrichment
- **Gene Ontology** (GO) enrichment
- **Disease associations** (Open Targets)
- **Drug-gene interactions** (DGIdb)

### Module 3: Network Visualization
- Interactive **Gene → Disease → Drug** networks
- Built with **Cytoscape.js**
- Zoom, pan, and hover interactions
- Export to PNG and JSON

## 🏗️ Architecture

```
Gene2Therapy/
├── backend/           # Django project settings
├── deg_analysis/      # DEG analysis API
├── annotation/        # External API integrations
├── networks/          # Network construction API
├── frontend/          # React (Vite) application
└── requirements.txt   # Python dependencies
```

## 🚀 Getting Started

### Prerequisites
- Python 3.10+
- Node.js 18+ & npm

### Backend Setup

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\Scripts\activate   # Windows

# Install dependencies
pip install -r requirements.txt

# Run migrations
python manage.py migrate

# Start server
python manage.py runserver
```

Backend runs at: `http://127.0.0.1:8000`

### Frontend Setup

```bash
cd frontend

# Install dependencies
npm install

# Start development server
npm run dev
```

Frontend runs at: `http://localhost:5173`

## 📡 API Endpoints

### DEG Analysis
| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/api/deg/jobs/` | Create analysis job |
| GET | `/api/deg/jobs/{id}/` | Get job status |
| GET | `/api/deg/jobs/{id}/results/` | Get results |

### Annotation
| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/api/annotation/gene-map/` | NCBI gene mapping |
| POST | `/api/annotation/kegg/` | KEGG pathway enrichment |
| POST | `/api/annotation/go/` | GO enrichment |
| POST | `/api/annotation/diseases/` | Disease associations |
| POST | `/api/annotation/drugs/` | Drug-gene interactions |

### Networks
| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/api/networks/gene-pathway-disease/` | Gene-Pathway-Disease network |
| POST | `/api/networks/gene-drug/` | Gene-Drug network |
| POST | `/api/networks/therapy/` | Full therapy network |

## 🛠️ Tech Stack

**Backend:**
- Django 4.2+
- Django REST Framework
- pandas, numpy, scipy
- scanpy (scRNA-seq)
- pyDESeq2 (raw counts)

**Frontend:**
- React 18
- Vite
- Plotly.js (visualizations)
- Cytoscape.js (networks)
- Axios (API client)

## 📊 Data Formats

### Input
- CSV or TSV files
- Bulk RNA-seq: Genes (rows) × Samples (columns)
- scRNA-seq: Gene × Cell matrix

### Output
- DEG tables (CSV)
- Volcano plots (PNG/SVG)
- Heatmaps (PNG/SVG)
- Network graphs (PNG/JSON)

## 📝 License

This project is intended for academic and research use.

## 🤝 Contributing

Contributions are welcome! Please open an issue or submit a pull request.
