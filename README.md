# 🧬 BioContext

BioContext is a lightweight bioinformatics tool for **gene annotation and pathway enrichment analysis**.  
It integrates **NCBI Entrez** and **KEGG REST API** to fetch gene metadata, descriptions, and pathway information.  
A simple **Streamlit web app** is also included for easy use without coding.

---

## 🚀 Features
- Load **gene lists** from `.csv` or `.txt`
- Retrieve:
  - **NCBI Gene ID**
  - **Gene Description**
  - **KEGG Pathways** (mapped with IDs and names)
- Perform **Pathway Enrichment Analysis** using the Hypergeometric test
- Download results as `.csv`
- Interactive **Streamlit Web Interface**

---

## 📦 Installation

Clone the repository:

```bash
git clone https://github.com/faiyazrizwee/Biocontext.git
cd Biocontext
```

Create and activate a conda environment:

```bash
conda env create -f environment.yml
conda activate biocontext
```

---

## 🖥️ Usage

### 1. Command-line mode
Run with a gene list:

```bash
python main.py gene_list.csv
```

Results will be saved in `data/gene_metadata_with_kegg.csv`.

---

### 2. Web interface (Streamlit)
Launch the Streamlit app:

```bash
streamlit run streamlit_app.py
```

Then open your browser at:  
👉 http://localhost:8501  

Upload your **gene list file**, select an organism, and download results.

---

## 📂 File structure
```
Biocontext/
│── main.py              # Command-line tool
│── streamlit_app.py     # Streamlit web interface
│── environment.yml      # Conda environment
│── data/                # Example input/output data
│   └── gene_list.txt
│── README.md            # Project documentation
```

---

## 🧪 Example Input (gene_list.txt)
```
TP53
BRCA1
EGFR
MYOSIN
```

## 📊 Example Output
**gene_metadata_with_kegg.csv**
```
Gene,NCBI_ID,Description,KEGG_Pathways
TP53,7157,tumor protein p53,path:hsa04115 - p53 signaling pathway
BRCA1,672,BRCA1 DNA repair associated,path:hsa04137 - Mitophagy; path:hsa04140 - Autophagy
...
```

---

## 🤝 Contributing
Pull requests are welcome. For major changes, please open an issue first.

---

## 📜 License
MIT License © 2025 Md Faiyaz Rizwee  
