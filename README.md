
# Gene2Therapy – BioContext Web Application

**Gene2Therapy** is a web-based bioinformatics application designed to analyze gene lists, enrich pathways, associate diseases, and suggest therapeutic drug repurposing. It integrates multiple data sources, such as **NCBI**, **KEGG**, and **OpenTargets**, to facilitate data-driven insights for gene-based research.

The app is optimized for **dark mode** and features an intuitive **Streamlit** interface, making it easy for researchers to perform complex bioinformatics analyses without requiring coding skills.

---

## 🚀 Features
- **Gene Annotation**: Fetches gene metadata from **NCBI**.
- **KEGG Pathway Enrichment**: Pathway analysis with **KEGG**.
- **Disease Associations**: Links genes to diseases using **OpenTargets**.
- **Drug Repurposing**: Suggests therapeutic drugs targeting specific genes.
- **Interactive Visualizations**: Displays pathway enrichment, gene-disease-drug networks, and more.
- **Rate-Limited API Calls**: Optimized with error handling and caching for efficiency.

---

## 📦 Installation

To get started with **Gene2Therapy**, follow these installation steps:

### 1. Clone the Repository

```bash
git clone https://github.com/faiyazrizwee/Gene2Therapy.git
cd Gene2Therapy
```

### 2. Set up the Conda Environment

```bash
conda env create -f environment.yml
conda activate gene2therapy
```

### 3. Install the Required Packages

```bash
pip install -r requirements.txt
```

### 4. Run the Application

```bash
streamlit run main.py
```

Open the application in your web browser at:
[http://localhost:8501](http://localhost:8501)

---

## 📂 File Structure
```
Gene2Therapy/
│── main.py              # Streamlit app with bioinformatics pipeline
│── requirements.txt     # Conda environment dependencies
│── environment.yml      # Conda environment configuration
│── assets/              # Application assets (e.g., logos, images)
│── README.md            # Project documentation
```

---

## 📧 Configuration
- The app requires a valid **NCBI Entrez email** for API access.
- You can configure the **organism** of interest from the dropdown, selecting from:
  - Homo sapiens (Human)
  - Mus musculus (Mouse)
  - Rattus norvegicus (Rat)

---

## 🧬 How It Works

### Gene Metadata and KEGG Pathways
1. **Upload Gene List**: Upload a CSV, TSV, Excel, or plain text file containing gene symbols, or manually input gene symbols.
2. **NCBI Gene Search**: The app searches **NCBI** for gene descriptions and associated KEGG pathways.
3. **KEGG Pathway Enrichment**: The app performs KEGG pathway enrichment analysis using the Hypergeometric test.

### Disease Associations and Drug Repurposing
4. **Disease Association**: Using **OpenTargets**, genes are mapped to their associated diseases with therapeutic scores.
5. **Drug Suggestions**: Based on the mapped genes and diseases, the app suggests drugs that are known to target those genes, showing their clinical trial phases and mechanisms of action.

### Visualizations
- **KEGG Pathway Enrichment**: Displays enriched pathways based on gene counts.
- **Gene-Disease-Drug Networks**: Visualizes the relationships between genes, diseases, and drugs in a Sankey diagram or network graph.

---

## 📊 Example Input and Output

### Example Input (gene_list.txt)
```
TP53
BRCA1
EGFR
MYC
```

### Example Output (gene_metadata_with_kegg.csv)
```
Gene,NCBI_ID,Description,KEGG_Pathways
TP53,7157,tumor protein p53,path:hsa04115 - p53 signaling pathway
BRCA1,672,BRCA1 DNA repair associated,path:hsa04137 - Mitophagy; path:hsa04140 - Autophagy
...
```

---

## 💡 Tips for Best Performance
- Keep gene lists under 100 genes for faster processing.
- The app supports **CSV**, **TSV**, **XLSX**, and plain **text** files.
- **Rate limiting** is implemented to avoid API timeouts; results are cached for 1 hour.

---

## 🤝 Contributing
We welcome contributions to improve **Gene2Therapy**. Please fork the repository, make changes, and submit a pull request. For major changes, open an issue first to discuss what you would like to change.

---

## 📜 License
MIT License © 2025 Md Faiyaz Rizwee
