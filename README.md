# Gene2Therapy - Integrated Gene Analysis

Gene2Therapy is a powerful, web-based bioinformatics tool designed to streamline the analysis of gene expression data and its integration with clinical and therapeutic insights.

## Features
- **Differential Gene Expression (DGE)**: Analyze RNA-seq data using t-tests or the industry-standard DESeq2 pipeline.
- **Pathway Analysis**: Identify biological pathways associated with differentially expressed genes.
- **NCBI Integration**: Fetch gene information and research data directly from NCBI biological databases.
- **Interactive Visualizations**: High-quality plots including Volcano plots, Heatmaps, and Network structures built with Plotly and NetworkX.
- **Modern UI**: A clean, intuitive interface designed for researchers.

## Installation
1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd <project-directory>
   ```
2. Install dependencies:
   ```bash
   pip install streamlit pandas numpy scipy pydeseq2 plotly networkx biopython requests
   ```
3. Run the application:
   ```bash
   streamlit run main.py
   ```
